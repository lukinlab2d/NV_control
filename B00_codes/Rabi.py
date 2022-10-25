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
from copy import deepcopy
from qcodes.actions import Task as qctask
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

class Rabi(Instrument):

    def __init__(self, name='RabiObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserParam =       {'delay_time': 2, 'channel':3}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.AFGParam =         {'delay_time': 2, 'channel':1}
        self.MWswitchParam =    {'delay_time': 2, 'channel':2}
        settings_extra = {'clock_speed': self.clock_speed, 'Laser': self.LaserParam, 'Counter': self.CounterParam, 
                        'AFG': self.AFGParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.tausArray = np.linspace(start, stop, num_sweep_points)
        self.uwPower = self.settings['uwPower']; self.uwFreq = self.settings['uwFreq']

        self.add_parameter(
            name = "sig",
            parameter_class  = Signal,
            settings = self.settings,
            measurementObject = self
        )
        self.add_parameter(
            name = "ref",
            parameter_class = Reference,
        )
        self.savedPulseSequencePlots = {}

        # SRS object
        self.srs = SRS()
        self.srs.set_freq(self.uwFreq) #Hz
        self.srs.set_RFAmplitude(self.uwPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(self.tausArray[0], self.tausArray[-1], num=len(self.tausArray)),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig, ref,
                                            qctask(sig.plotPulseSequences),
                                            )

        data = loop.get_data_set(name='Rabi')
        data.add_metadata(self.settings)
        
        plot = QtPlot(
            data.RabiObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.RabiObject_ref, name='ref')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        img = Image.open(dataPlotFile)
        img.show()
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.RabiObject = measurementObject
        self.loopCounter = 0
        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.tausArray = np.linspace(start, stop, num_sweep_points)

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))

        self.ctrtask.stop(); self.ctrtask.close()

        rate = xLineData/(self.read_duration_in_ns/1e9)/1e3
        sig = rate[::self.num_reads_per_iter]
        ref = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)
        return sig_avg

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster

        print("Loop " + str(self.loopCounter))
        num_loops                 = self.settings['num_loops'];               
        laser_init_delay_in_ns    = self.settings['laser_init_delay_in_ns'];    laser_init_duration_in_ns  = self.settings['laser_init_duration_in_ns']; when_init_end = laser_init_delay_in_ns+laser_init_duration_in_ns
        AFG_delay_after_init_in_ns= self.settings['AFG_delay_after_init_in_ns'];AFG_duration_in_ns         = tau_ns; 
        AFG_delay_in_ns = when_init_end + AFG_delay_after_init_in_ns
        
        when_pulse_end = AFG_delay_in_ns + AFG_duration_in_ns
        laser_read_delay_after_pulse_in_ns   = self.settings['laser_read_delay_after_pulse_in_ns'];  laser_read_duration_in_ns  = self.settings['laser_read_duration_in_ns']
        read_signal_delay_after_pulse_in_ns  = self.settings['read_signal_delay_after_pulse_in_ns']; read_signal_duration_in_ns = self.settings['read_signal_duration_in_ns']
        read_ref_delay_after_pulse_in_ns     = self.settings['read_ref_delay_after_pulse_in_ns'];    read_ref_duration_in_ns    = self.settings['read_ref_duration_in_ns']
        
        laser_read_delay_in_ns = when_pulse_end + laser_read_delay_after_pulse_in_ns
        read_signal_delay_in_ns = when_pulse_end + read_signal_delay_after_pulse_in_ns
        read_ref_delay_in_ns = when_pulse_end + read_ref_delay_after_pulse_in_ns

        self.read_duration_in_ns = read_signal_duration_in_ns
        if read_signal_duration_in_ns != read_ref_duration_in_ns:
            raise Exception("Duration of reading signal and reference must be the same")

        pulseSequence = []
        pulseSequence += [spc.Pulse('Laser',   laser_init_delay_in_ns, duration=int(laser_init_duration_in_ns))] # times are in ns
        pulseSequence += [spc.Pulse('AFG',     AFG_delay_in_ns,        duration=int(AFG_duration_in_ns))] # times are in ns
        pulseSequence += [spc.Pulse('MWswitch',AFG_delay_in_ns,        duration=int(AFG_duration_in_ns))]
        pulseSequence += [spc.Pulse('Laser',   laser_read_delay_in_ns, duration=int(laser_read_duration_in_ns))] # times are in ns
        pulseSequence += [spc.Pulse('Counter', read_signal_delay_in_ns,duration=int(read_signal_duration_in_ns))] # times are in ns
        pulseSequence += [spc.Pulse('Counter', read_ref_delay_in_ns,   duration=int(read_ref_duration_in_ns))] # times are in ns
        self.pulseSequence = pulseSequence
        
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        self.pb.program_pb(pulseSequence, num_loops=num_loops)

        num_reads_per_iter = 0
        for pulse in pulseSequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        self.num_reads = int(num_loops * num_reads_per_iter)
        self.num_reads_per_iter = num_reads_per_iter

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
            samps_per_chan = 2*self.num_reads # x2 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        print('Set tau to ' + str(tau_ns) + " ns")
    
    def plotPulseSequences(self):
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulseSequence, ifShown=True, ifSave=False)
                fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.RabiObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1


class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg