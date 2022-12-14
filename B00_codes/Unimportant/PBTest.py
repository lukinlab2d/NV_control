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

from operator import contains
from typing import Dict, Union
from collections import namedtuple
import numpy as np
import datetime as dt
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS

import nidaqmx, datetime, ctypes, os, warnings, time, itertools
from nidaqmx.constants import *
from qcodes.instrument.base import Instrument
# from qcodes.instrument.parameter import ParameterWithSetpoints, Parameter
# from qcodes.utils.validators import *
from pylabcontrol.core.read_write_functions import get_config_value
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
import matplotlib as mpl
import matplotlib.pyplot as plt


from qcodes.dataset import do1d, do2d, dond, LinSweep, LogSweep, ArraySweep #TogetherSweep
from qcodes.utils.dataset.doNd import plot
from qcodes.dataset.sqlite.database import initialise_or_create_database_at
from qcodes.dataset.experiment_container import load_or_create_experiment
from qcodes.tests.instrument_mocks import DummyInstrument, DummyInstrumentWithMeasurement
from qcodes.dataset.measurements import Measurement
from qcodes.dataset.plotting import plot_dataset

def ns2cycles(time_in_ns, samp_rate=1e7):
        return int(time_in_ns/1e9*samp_rate)
    
class Trace:
    def __init__(self, arr, name):
        self.name = name
        self.arr = arr

        if "aser" in name: 
            self.vert_offset = 3; self.color = 'C2'
        elif "AFG" in name or "MWswitch" in name: 
            self.vert_offset = 1.5; self.color = 'C0'
        elif "ounter" in name: 
            self.vert_offset = 0; self.color = 'k'

        self.arr = [x + self.vert_offset for x in arr]
        self.length = len(arr)

####################################################################################################################
if __name__ == '__main__':
    print()
    now = dt.datetime.now(); date = now.strftime("%Y-%m-%d")
    print(now)
    
    for i in range(1000):
        # dataFolder = os.path.join(os.getcwd(), 'data', date)
        # if not os.path.isdir(dataFolder): os.mkdir(dataFolder)
        # dataFilePath = os.path.join(dataFolder, 'test.db')
        # print(dataFilePath)
        # initialise_or_create_database_at(dataFilePath)

        # clock speed is in MHz - is 'status' needed in the dictionary?
        clock_speed = 500 # MHz
        LaserParam = {'delay_time': 2, 'channel':1}
        CounterParam = {'delay_time': 2, 'channel':2}

        num_loops = int(1e6)
        laser_delay_in_ns  = 0.2e3;  laser_duration_in_ns   = 20
        read_delay_in_ns   = 0.2e3;  read_duration_in_ns    = 20

        uwFreq = 400e6; uwPower = 0
        srs = SRS()
        srs.set_freq(uwFreq) #Hz
        srs.set_RFAmplitude(uwPower) #dBm
        srs.enableIQmodulation()
        srs.enable_RFOutput()

        settings = {'clock_speed': clock_speed, 'Laser': LaserParam, 'Counter': CounterParam, 
                    'PB_type': 'USB', 'min_pulse_dur': int(1*1e3/clock_speed)}

        pb = spc.B00PulseBlaster("SpinCorePB", settings=settings, verbose=False)
        # print("Is pb an Instrument? " + str(isinstance(pb, Instrument)))
        
        pulse_sequence = []
        pulse_sequence += [spc.Pulse('Laser',  laser_delay_in_ns,  duration=int(laser_duration_in_ns))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter', read_delay_in_ns,  duration=int(read_duration_in_ns))] # times are in ns
        # total_duration_in_ns = (read_delay_in_ns + read_duration_in_ns + 100)*num_loops - read_delay_in_ns
        pb.program_pb(pulse_sequence, num_loops=num_loops)

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        num_reads = int(num_loops * num_reads_per_iter)

        # total_duration_in_DAQ_cycles = int(total_duration_in_ns / 1e9 * samp_rate)
        # laser_duration_in_DAQ_cycles = int(laser_duration_in_ns / 1e9 * samp_rate)
        # read_duration_in_DAQ_cycles = int(read_duration_in_ns / 1e9 * samp_rate)

        # # Pulse width counter. Timebase = signal; gate = PB signal
        # ctrtask = nidaqmx.Task()
        # pulseWidthChan = ctrtask.ci_channels.add_ci_pulse_width_chan( # define the pulse width counter
        #     counter = "cDAQ1Mod1/ctr0",
        #     name_to_assign_to_channel = "",
        #     min_val = 0,
        #     max_val = int(1e8),
        #     units = TimeUnits.TICKS,
        #     starting_edge = Edge.RISING,
        #     )
        # ctrtask.timing.cfg_implicit_timing(
        #     sample_mode = AcquisitionType.CONTINUOUS,
        #     samps_per_chan = 2*num_reads # x2 to make sure buffer doesn't overflow
        #     )
        # pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        # pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string
        
        # ctrtask.start()

        pb.start_pulse_seq()
        pb.wait()
        pb.stop_pulse_seq()
        pb.close()

        # xLineData = np.array(ctrtask.read(num_reads, timeout=0.5))

        # ctrtask.stop(); ctrtask.close()
        
        # rate = xLineData/(read_duration_in_ns/1e9)/1e3
        # signal = rate[::num_reads_per_iter]
        # ref = rate[1::num_reads_per_iter]

        # # Plot results
        # x_axis = np.array(range(0, len(xLineData)))
    
        # fig, ax = plt.subplots()
        # ax.plot(x_axis[::num_reads_per_iter], signal,'r', label='Sig')
        # ax.plot(x_axis[1::num_reads_per_iter], ref,'k', label='Ref')
        # ax.set_xlabel('Iterations')
        # ax.set_ylabel('Count rate (kc/s)')
        # rate_avg = np.average(xLineData)/(read_duration_in_ns/1e9)/1e3
        # print('average rate = ' + str(rate_avg) + ' kc/s')
        # ax.legend(loc='best')

        # Plot pulses
        # samp_rate = 1e7
        # fig2, ax = plt.subplots()
        # PulseTrace = {}
        # maxLength = np.max([pulse.start_time + pulse.duration for pulse in pulse_sequence])
        
        # for pulse in pulse_sequence:
        #     channel_id = pulse.channel_id
        #     trace_array = np.concatenate((np.zeros(ns2cycles(pulse.start_time)),
        #                             np.ones(ns2cycles(pulse.duration)),
        #                             np.zeros(ns2cycles(maxLength - pulse.start_time - pulse.duration))))
        #     if not channel_id in PulseTrace: PulseTrace[channel_id] = Trace(trace_array, channel_id)
        #     else: PulseTrace[channel_id].arr = np.add(PulseTrace[channel_id].arr, trace_array)
        
        # for channel_id in PulseTrace:
        #     tr = PulseTrace[channel_id]
        #     x_axis = np.array(range(tr.length))/int(samp_rate)*1e6
        #     if tr.name != "MWswitch": ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color)
        #     ax.set_xlabel("Time ($\mu$s)")
        #     ax.legend(loc='best')
        
        # plt.show(block=False)

        # # testing quick plotting with matplotlib
        # ax1=ax
        # fig3, ax = plt.subplots()
        # x_axis = np.array(range(0, len(xLineData)))
        # for i in range(100):
        #     if np.mod(i,2) == 0:
        #         ax.plot(x_axis[::num_reads_per_iter], signal+i,'r', label='Sig')
        #     else:
        #         ax.plot(x_axis[1::num_reads_per_iter], ref+i,'k', label='Ref')
        #         # ax1.plot(x_axis[1::num_reads_per_iter], ref+i,'k', label='Ref')
        #     plt.pause(1e-3)

        # plt.show()
    # print('Duration =  ' + str(dt.datetime.now()-now))