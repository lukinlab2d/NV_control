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
import numpy as np
from nidaqmx.constants import *
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
from PlotPulse import *
from T1 import *

####################################################################################################################

start = 100; stop = 3000; num_sweep_points = 59
tausArray = np.linspace(start, stop, num_sweep_points)
uwPower = -50; uwFreq = 2.87e9

num_loops                                   = int(200000);            
laser_init_delay_in_ns                      = 1e3;      laser_init_duration_in_ns  = 1e3
AFG_delay_after_init_in_ns                  = 40;       pi_time                    = 100
laser_read_duration_in_ns                   = 2e3
read_signal_delay_after_read_laser_on_in_ns = 500;      read_signal_duration_in_ns = 300
read_ref_delay_after_read_laser_on_in_ns    = 1500;     read_ref_duration_in_ns    = 300


settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops,
            'uwPower': uwPower, 'uwFreq': uwFreq,
            'laser_init_delay_in_ns':                     laser_init_delay_in_ns,                      'laser_init_duration_in_ns': laser_init_duration_in_ns,
            'AFG_delay_after_init_in_ns':                 AFG_delay_after_init_in_ns,                  'pi_time':                   pi_time,
            'laser_read_duration_in_ns':                  laser_read_duration_in_ns,
            'read_signal_delay_after_read_laser_on_in_ns':read_signal_delay_after_read_laser_on_in_ns,'read_signal_duration_in_ns': read_signal_duration_in_ns,
            'read_ref_delay_after_read_laser_on_in_ns':   read_ref_delay_after_read_laser_on_in_ns,   'read_ref_duration_in_ns':    read_ref_duration_in_ns}

start = time.time()
T1Object = T1(settings=settings, ifPlotPulse=True) # this is implemented as an Instrument
T1Object.runScan()
print('Total time = ' + str(time.time() - start) + ' s')


