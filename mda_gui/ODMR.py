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
import datetime as dt
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS

import nidaqmx
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


from qcodes.dataset import do1d, do2d, dond, LinSweep, LogSweep, ArraySweep
from qcodes.utils.dataset.doNd import plot
from qcodes.dataset.sqlite.database import initialise_or_create_database_at
from qcodes.dataset.experiment_container import load_or_create_experiment
from qcodes.tests.instrument_mocks import DummyInstrument, DummyInstrumentWithMeasurement
from qcodes.dataset.measurements import Measurement
from qcodes.dataset.plotting import plot_dataset

def ns2cycles(time_in_ns, samp_rate=1e7):
        return int(time_in_ns/1e9*samp_rate)
    

class ODMR(Instrument):

    def __init__(self, name='ODMRObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        self.ifPlotPulse=ifPlotPulse

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserParam = {'delay_time': 2, 'channel':3}
        self.CounterParam = {'delay_time': 2, 'channel':4}
        self.AFGParam = {'delay_time': 2, 'channel':1}
        self.MWswitchParam = {'delay_time': 2, 'channel':2}

        settings_extra = {'clock_speed': self.clock_speed, 'Laser': self.LaserParam, 'Counter': self.CounterParam, 
                        'AFG': self.AFGParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'min_pulse_dur': int(5*1e3/self.clock_speed)}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.freqsArray = np.linspace(start, stop, num_sweep_points)

        num_loops                = self.settings['num_loops'];               #samp_rate_DAQ              = self.settings['samp_rate_DAQ']
        laser_init_delay_in_ns   = self.settings['laser_init_delay_in_ns'];  laser_init_duration_in_ns  = self.settings['laser_init_duration_in_ns']
        laser_read_delay_in_ns   = self.settings['laser_read_delay_in_ns'];  laser_read_duration_in_ns  = self.settings['laser_read_duration_in_ns']
        AFG_delay_in_ns          = self.settings['AFG_delay_in_ns'];         AFG_duration_in_ns         = self.settings['AFG_duration_in_ns']
        read_signal_delay_in_ns  = self.settings['read_signal_delay_in_ns']; read_signal_duration_in_ns = self.settings['read_signal_duration_in_ns']
        read_ref_delay_in_ns     = self.settings['read_ref_delay_in_ns'];    read_ref_duration_in_ns    = self.settings['read_ref_duration_in_ns']
        
        if read_signal_duration_in_ns != read_ref_duration_in_ns:
            raise Exception("Duration of reading signal and reference must be the same") 

        print()
        self.now = dt.datetime.now(); self.date = self.now.strftime("%Y-%m-%d")
        print(self.now)   
    
        pulse_sequence = []
        pulse_sequence += [spc.Pulse('Laser',  laser_init_delay_in_ns,  duration=int(laser_init_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('Laser',  laser_read_delay_in_ns,  duration=int(laser_read_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('AFG', AFG_delay_in_ns, duration=int(AFG_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('MWswitch', AFG_delay_in_ns, duration=int(AFG_duration_in_ns))]
        pulse_sequence += [spc.Pulse('Counter', read_signal_delay_in_ns,   duration=int(read_signal_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter', read_ref_delay_in_ns,  duration=int(read_ref_duration_in_ns))] # times are in ns
        self.pulse_sequence = pulse_sequence
        
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings)
        self.pb.program_pb(pulse_sequence, num_loops=num_loops)

        self.srs = SRS()
        global srs; srs = self.srs
        self.srs.set_freq(3e9) #Hz
        self.srs.set_RFAmplitude(-30) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        num_reads = int(num_loops * num_reads_per_iter)


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
            samps_per_chan = 2*num_reads # x2 to make sure buffer doesn't overflow
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
        )

        self.add_parameter(
            name = "ref",
            parameter_class = Reference,
        )

        global pb; pb = self.pb
        global ctrtask; ctrtask = self.ctrtask
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref
        loop = Loop(
            sig.sweep(self.freqsArray[0], self.freqsArray[-1], num=len(self.freqsArray)),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref).then(qctask(sig.close))

        data = loop.get_data_set(name='ODMR')
        data.add_metadata(self.settings)
        
        plot = QtPlot(
            data.ODMRObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.ODMRObject_ref, name='ref')

        loop.with_bg_task(plot.update)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        img = Image.open(dataPlotFile)
        img.show()
        
        if self.ifPlotPulse:
            pulsePlotFilename = data.location + "/pulsePlot.png"
            plotPulseObject = PlotPulse(measurementObject=self, plotFilename=pulsePlotFilename, ifShown=True)
            plotPulseObject.makePulsePlot()

    
class Signal(Parameter):
    def __init__(self, num_reads: int, num_reads_per_iter: int, num_loops: int,
                 read_duration_in_ns: int, name='sig',**kwargs):
        super().__init__(name, **kwargs)
        self.num_reads = num_reads
        self.num_reads_per_iter = num_reads_per_iter
        self.num_loops = num_loops
        self.read_duration_in_ns = read_duration_in_ns
        self.loopCounter = 0

    def get_raw(self):
        self.loopCounter += 1
        ctrtask.start()

        pb.start_pulse_seq()
        pb.wait()
        pb.stop_pulse_seq_without_closing()

        xLineData = np.array(ctrtask.read(self.num_reads, timeout=0.5))

        ctrtask.stop()
        
        rate = xLineData/(self.read_duration_in_ns/1e9)/1e3
        sig = rate[::self.num_reads_per_iter]
        ref = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)
        return sig_avg
    
    def set_raw(self, value):
        srs.set_freq(value)
        print("Loop " + str(self.loopCounter))
        print("Set SRS freq to " + str(np.round(value/1e6,3)) + ' MHz') # set the SRS frequency in reality

    def close(self):
        ctrtask.close()
        pb.stop_pulse_seq()

class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg