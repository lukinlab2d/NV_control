"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
"""
from  qcodes.actions import Task as qctask
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot
import numpy as np
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS

import nidaqmx, time
from nidaqmx.constants import *
from qcodes.instrument.base import Instrument
from qcodes.instrument.parameter import Parameter
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
from PIL import Image
from PlotPulse import *


class ODMR_CW(Instrument):

    def __init__(self, name='ODMRObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        self.ifPlotPulse=ifPlotPulse

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserParam =       {'delay_time': 2, 'channel':3}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.AFGParam =         {'delay_time': 2, 'channel':1}
        self.MWswitchParam =    {'delay_time': 2, 'channel':2}
        global laserChannel; laserChannel = self.LaserParam['channel']

        settings_extra = {'clock_speed': self.clock_speed, 'Laser': self.LaserParam, 'Counter': self.CounterParam, 
                        'AFG': self.AFGParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'min_pulse_dur': int(5*1e3/self.clock_speed)}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        # List of frequencies and power
        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = self.settings['uwPower']

        # Pulse lengths
        num_loops                = self.settings['num_loops'];                wait_btwn_sig_ref          = self.settings['wait_btwn_sig_ref']
        AFG_delay_in_ns          = self.settings['AFG_delay_in_ns'];          AFG_duration_in_ns         = self.settings['AFG_duration_in_ns'];                             when_AFG_ends         = AFG_delay_in_ns + AFG_duration_in_ns      
        read_signal_delay_in_ns  = AFG_delay_in_ns;                           read_signal_duration_in_ns = AFG_duration_in_ns + self.settings['AFG_off_to_read_signal_off'];when_read_signal_ends = read_signal_delay_in_ns + read_signal_duration_in_ns 
        read_ref_delay_in_ns     = when_read_signal_ends + wait_btwn_sig_ref; read_ref_duration_in_ns    = read_signal_duration_in_ns;                                      when_read_ref_ends    = read_ref_delay_in_ns + read_ref_duration_in_ns
        laser_delay_in_ns        = self.settings['laser_delay_in_ns'];        laser_duration_in_ns       = when_read_ref_ends - laser_delay_in_ns + 100
        
        if read_signal_duration_in_ns != read_ref_duration_in_ns:
            raise Exception("Duration of reading signal and reference must be the same")    

        # Make pulse sequence (per each freq)
        pulse_sequence = []
        pulse_sequence += [spc.Pulse('Laser',  laser_delay_in_ns,  duration=int(laser_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('AFG', AFG_delay_in_ns, duration=int(AFG_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('MWswitch', AFG_delay_in_ns, duration=int(AFG_duration_in_ns))]
        pulse_sequence += [spc.Pulse('Counter', read_signal_delay_in_ns,   duration=int(read_signal_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter', read_ref_delay_in_ns,  duration=int(read_ref_duration_in_ns))] # times are in ns
        self.pulse_sequence = pulse_sequence
        
        # SRS object
        self.srs = SRS()
        self.srs.set_freq(3e9) #Hz
        self.srs.set_RFAmplitude(uwPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        num_reads = int(num_loops * num_reads_per_iter) # the results look like this: sig-ref-sig-ref-...-sig-ref for num_loops times

        # Pulse width counter. Timebase = signal; gate = PB signal
        self.ctrtask = nidaqmx.Task()
        pulseWidthChan = self.ctrtask.ci_channels.add_ci_pulse_width_chan( # define the pulse width counter
            counter = "cDAQ1Mod1/ctr0",
            name_to_assign_to_channel = "",
            min_val = 0,
            max_val = int(1e8),
            units = TimeUnits.TICKS,
            starting_edge = Edge.RISING,
            )
        self.ctrtask.timing.cfg_implicit_timing(
            sample_mode = AcquisitionType.CONTINUOUS,
            samps_per_chan = int(1.2*num_reads) # x2 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        self.add_parameter(
            name = "sig",
            parameter_class = Signal,
            num_reads = num_reads,
            num_loops = num_loops,
            num_reads_per_iter = num_reads_per_iter,
            read_duration_in_ns = read_signal_duration_in_ns,
            pulse_sequence = pulse_sequence,
            settings = self.settings
        )
        self.add_parameter(
            name = "ref",
            parameter_class = Reference,
        )
        self.add_parameter(
            name = "sigOverRef",
            parameter_class = SigOverRef,
        )

        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        global ctrtask; ctrtask = self.ctrtask
        global srs; srs = self.srs
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigOverRef = self.sigOverRef

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(self.freqsArray[0], self.freqsArray[-1], num=len(self.freqsArray)),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,sigOverRef).then(qctask(sig.close))

        data = loop.get_data_set(name='ODMR_CW')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.ODMRObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.ODMRObject_ref, name='ref')
        # plot = QtPlot(
        #     data.ODMRObject_sigOverRef, # this is implemented as a Parameter
        #     figsize = (1200, 600),
        #     interval = 1,
        #     name = 'sig/ref'
        #     )

        loop.with_bg_task(plot.update)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        # img = Image.open(dataPlotFile)
        # img.show()
        
        if self.ifPlotPulse:
            pulsePlotFilename = data.location + "/pulsePlot.png"
            plotPulseObject = PlotPulse(measurementObject=self, plotFilename=pulsePlotFilename, ifShown=True)
            plotPulseObject.makePulsePlot()
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ODMRObject_sig_set.dat'

    
class Signal(Parameter):
    def __init__(self, num_reads: int, num_reads_per_iter: int, num_loops: int,
                 read_duration_in_ns: int, name='sig', pulse_sequence=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.num_reads = num_reads
        self.num_reads_per_iter = num_reads_per_iter
        self.num_loops = num_loops
        self.read_duration_in_ns = read_duration_in_ns
        self.loopCounter = 0
        self.pulse_sequence = pulse_sequence
        self.settings = settings

    def get_raw(self):
        self.loopCounter += 1
        ctrtask.start()

        # Pulse Blaster object inittialized here to avoid long delay between initialization and first run
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        self.pb.program_pb(self.pulse_sequence, num_loops=self.num_loops)
        pb = self.pb

        pb.start_pulse_seq()
        pb.wait()
        pb.stop_pulse_seq()
        pb.close()

        xLineData = np.array(ctrtask.read(self.num_reads, timeout=10))

        ctrtask.stop()
        
        rate = xLineData/(self.read_duration_in_ns/1e9)/1e3
        sig = rate[::self.num_reads_per_iter]
        ref = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)
        global sig_avg_over_ref_avg; sig_avg_over_ref_avg = sig_avg/ref_avg
        return sig_avg
    
    def set_raw(self, value):
        srs.set_freq(value)
        print("Loop " + str(self.loopCounter))
        print("Set SRS freq to " + str(np.round(value/1e6,3)) + ' MHz') # set the SRS frequency in reality

    def close(self):
        ctrtask.close()
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        pb.turn_on_infinite(channel=laserChannel)

class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg

class SigOverRef(Parameter):
    def __init__(self, name='sigOverRef',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_avg_over_ref_avg