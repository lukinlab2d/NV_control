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
from ODMR_CW import *
from PlotPulse import *
import dataReader

####################################################################################################################
cw = 1
if cw != 1:
    # Dataset 1
    for i in np.linspace(800,2000,1):
        # Test for Pulsed ODMR
        start = 2.85e9; stop = 2.89e9; num_sweep_points = 81
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = -30

        num_loops               = int(1e6)
        laser_init_delay_in_ns  = 10;       laser_init_duration_in_ns = 1e3
        laser_to_AFG_delay      = 1000;     AFG_duration_in_ns        = 100
        laser_to_DAQ_delay      = 1400;     read_duration             = 300
        DAQ_to_laser_off_delay  = 10000

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                    'laser_init_delay_in_ns': laser_init_delay_in_ns,'laser_init_duration_in_ns': laser_init_duration_in_ns,
                    'laser_to_AFG_delay':     laser_to_AFG_delay ,   'AFG_duration_in_ns':        AFG_duration_in_ns,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay}

        start = time.time()
        ODMRObject = ODMR(settings=settings, ifPlotPulse=True) # this is implemented as an Instrument
        ODMRObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = ODMRObject.getDataFilename()
        dataReader.readData(dataFilename)
        ODMRObject.close()

else:
    # Test for CW ODMR
    # Dataset 
    for i in range(1):
        start = 2.855e9; stop = 2.885e9; num_sweep_points = 61
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = -50

        num_loops = int(1e6); wait_btwn_sig_ref = 20e3
        AFG_delay_in_ns = 1000; AFG_duration_in_ns = 10e3
        laser_delay_in_ns = 10; AFG_off_to_read_signal_off = 0

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                    'AFG_delay_in_ns':AFG_delay_in_ns, 'AFG_duration_in_ns':AFG_duration_in_ns, 
                    'AFG_off_to_read_signal_off':AFG_off_to_read_signal_off,
                    'wait_btwn_sig_ref':wait_btwn_sig_ref, 'laser_delay_in_ns':laser_delay_in_ns
                    }

        start = time.time()
        ODMRObject = ODMR_CW(settings=settings, ifPlotPulse=True) # implemented as Instrument. Turn of plot if pulse length > 1e5 ns
        ODMRObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = ODMRObject.getDataFilename()
        dataReader.readData(dataFilename)
        ODMRObject.close()






