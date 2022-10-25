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
import matplotlib as mpl
import matplotlib.pyplot as plt
from ODMR import *  
from PlotPulse import *

####################################################################################################################

start = 2.8e9; stop = 3e9; num_sweep_points = 21
freqsArray = np.linspace(start, stop, num_sweep_points)
uwPower = -15

num_loops               = int(200000);            #samp_rate_DAQ              = 1e7 # not important
laser_init_delay_in_ns  = 1e3;                  laser_init_duration_in_ns  = 1e3; when_init_end = laser_init_delay_in_ns+laser_init_duration_in_ns
AFG_delay_in_ns         = when_init_end+40;     AFG_duration_in_ns         = 40; when_pulse_end = AFG_delay_in_ns+AFG_duration_in_ns
laser_read_delay_in_ns  = when_pulse_end;       laser_read_duration_in_ns  = 2e3
read_signal_delay_in_ns = when_pulse_end+500;   read_signal_duration_in_ns = 300
read_ref_delay_in_ns    = when_pulse_end+1500;  read_ref_duration_in_ns    = 300


settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
            'laser_init_delay_in_ns': laser_init_delay_in_ns,'laser_init_duration_in_ns': laser_init_duration_in_ns,
            'laser_read_delay_in_ns': laser_read_delay_in_ns,'laser_read_duration_in_ns': laser_read_duration_in_ns,
            'AFG_delay_in_ns':AFG_delay_in_ns, 'AFG_duration_in_ns':AFG_duration_in_ns,
            'read_signal_delay_in_ns':read_signal_delay_in_ns, 'read_signal_duration_in_ns':read_signal_duration_in_ns,
            'read_ref_delay_in_ns':read_ref_delay_in_ns, 'read_ref_duration_in_ns':read_ref_duration_in_ns}


ODMRObject = ODMR(settings=settings, ifPlotPulse=True) # this is implemented as an Instrument
ODMRObject.runScan()


