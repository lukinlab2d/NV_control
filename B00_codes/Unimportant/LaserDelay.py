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

from typing import Dict, Union
from collections import namedtuple
import numpy as np
import datetime as dt
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc

import nidaqmx, datetime, ctypes, os, warnings, time, itertools
import nidaqmx.constants
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

####################################################################################################################
if __name__ == '__main__':
    print()
    print(dt.datetime.now())

    # clock speed is in MHz - is 'status' needed in the dictionary?
    LaserParam = {'delay_time': 2, 'channel':3}
    CtrCLKParam = {'delay_time': 2, 'channel':4}
    clock_speed = 500 # MHz
    settings = {'clock_speed': clock_speed, 'Laser': LaserParam, 'CtrCLK': CtrCLKParam, 'PB_type': 'USB',
                'min_pulse_dur': int(5*1e3/clock_speed)}

    pb = spc.B00PulseBlaster("SpinCorePB", settings=settings)

    num_loops = int(1e6); 
    laser_delay_in_ns = 1e3; laser_duration_in_ns = 1e3
    read_delay_in_ns = 1e3; read_duration_in_ns = 1e3
    
    pulse_sequence = []
    pulse_sequence += [spc.Pulse('Laser', laser_delay_in_ns, duration=int(laser_duration_in_ns))] # times are in ns
    pulse_sequence += [spc.Pulse('CtrCLK', read_delay_in_ns, duration=int(read_duration_in_ns))]
    total_duration_in_ns = (read_delay_in_ns + read_duration_in_ns + 100)*num_loops - read_delay_in_ns
    # half_cycle = 10
    # for ith_cycle in range(2000):
    #     pulse_sequence += [spc.Pulse('CtrCLK', 10e3+2*ith_cycle*half_cycle, duration=int(half_cycle))]
    pb.program_pb(pulse_sequence, num_loops=num_loops)

    samp_rate = 2e5
    total_duration_in_DAQ_cycles = int(total_duration_in_ns / 1e9 * samp_rate)
    laser_duration_in_DAQ_cycles = int(laser_duration_in_ns / 1e9 * samp_rate)
    read_duration_in_DAQ_cycles = int(read_duration_in_ns / 1e9 * samp_rate)
    print(total_duration_in_DAQ_cycles)
    print(laser_duration_in_DAQ_cycles)
    print(read_duration_in_DAQ_cycles)

    # the clock to sync AO and counter
    clktask = nidaqmx.Task() 
    clktask.co_channels.add_co_pulse_chan_freq(  # adding dig pulse train chan
        counter = "cDAQ1Mod1/ctr1",
        name_to_assign_to_channel = "",
        units = nidaqmx.constants.FrequencyUnits.HZ,
        idle_state = nidaqmx.constants.Level.LOW,
        initial_delay = 0.0,
        freq = samp_rate,
        duty_cycle = 0.5
        )
    clktask.timing.cfg_implicit_timing( # implicit timing by the hardware
        sample_mode = AcquisitionType.CONTINUOUS, # the clock should run continuously in principle
        samps_per_chan = int(total_duration_in_DAQ_cycles + 1) # does this matter?
        )

    ctrtask = nidaqmx.Task()
    ctrtask.ci_channels.add_ci_count_edges_chan( # define the counter
        counter = "cDAQ1Mod1/ctr0",
        name_to_assign_to_channel = "",
        edge = nidaqmx.constants.Edge.RISING,
        initial_count = 0,
        count_direction = nidaqmx.constants.CountDirection.COUNT_UP
        )
    ctrtask.timing.cfg_samp_clk_timing( # cfg sample clk timing
        rate = samp_rate,
        source = "/cDAQ1/Ctr1InternalOutput", # the clock defined above
        active_edge = nidaqmx.constants.Edge.RISING,
        sample_mode = AcquisitionType.FINITE,
        samps_per_chan = int(total_duration_in_DAQ_cycles + 1) 
        )
    
    clktask.triggers.start_trigger.cfg_dig_edge_start_trig(trigger_source="/cDAQ1Mod1/PFI1", trigger_edge=Edge.RISING)
    ctrtask.start()
    clktask.start()
    
    pb.start_pulse_seq()
    pb.wait()
    pb.stop_pulse_seq()

    xLineData = ctrtask.read(total_duration_in_DAQ_cycles + 1)  # +1 is to take the difference later
    ctrtask.wait_until_done()
    ctrtask.stop(); ctrtask.close()
    clktask.stop(); clktask.close()
    diffData = np.diff(xLineData); x_axis = np.array(range(0, len(diffData)))/samp_rate*1e6

    fig, ax = plt.subplots()
    ax.plot(x_axis, diffData)
    ax.set_xlabel('us')
    rate = sum(diffData)/x_axis[-1]*1e6/1e3
    print('sum(diffData) = ' + str(sum(diffData)))
    print('rate = ' + str(rate) + ' kc/s')

    laserPulse = np.concatenate((np.repeat([0],int(10e3/1e9*samp_rate)), 
                                np.repeat([1],laser_duration_in_DAQ_cycles),
                                np.repeat([0],read_duration_in_DAQ_cycles-laser_duration_in_DAQ_cycles),
                                np.repeat([0],int(100/1e9*samp_rate))))
    laserPulse = np.tile(laserPulse,num_loops-1)
    laserPulse = np.concatenate((np.repeat([1],laser_duration_in_DAQ_cycles),
                                np.repeat([0],read_duration_in_DAQ_cycles-laser_duration_in_DAQ_cycles),
                                np.repeat([0],int(100/1e9*samp_rate)),
                                laserPulse))
    x_axis2 = np.array(range(0, len(laserPulse)))/samp_rate*1e6
    # ax.plot(x_axis2, laserPulse, 'red')
    plt.show()






    # num_daq_reads = 0
    # for pulse in pulse_sequences[0]:
    #     if pulse.channel_id == 'CtrCLK':
    #         num_daq_reads += 1
    # signal = [0.0]
    # norms = np.repeat([0.0], (num_daq_reads - 1))
    # count_data = np.repeat([np.append(signal, norms)], len(pulse_sequences), axis=0)

    # print(signal)
    # print(norms)
    # print(count_data)
    # need a trigger pulse to do NIDAQ.read()

    

