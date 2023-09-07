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
from B00_codes.ODMR import *  
from B00_codes.ODMR_CW import *
from B00_codes.PlotPulse import *
import B00_codes.dataReader as dataReader

####################################################################################################################


# For NV tracking
if True: 
    if_tracking = 0
    xy_scan_read_time      = 5;      xy_scan_settle_time    = 0.2;  
    xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
    x_minus_range          = 0.05;   x_plus_range           = 0.05
    y_minus_range          = 0.05;   y_plus_range           = 0.05
    xy_displacement_limit  = 0.02;   num_of_scans           = 5;    tracking_period = 100

    xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
    x_minus_range          = 0.05;   x_plus_range           = 0.05
    z_minus_range          = 0.3;    z_plus_range           = 0.3
    xz_displacement_limit  = 0.1; 

trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
                    'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
                    'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                    'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
                    'xy_displacement_limit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
                    'tracking_period':        tracking_period,       'if_tracking':            if_tracking,
                    'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                    'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                    'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                    'xz_displacement_limit':  xz_displacement_limit,}

ifCW = 1

if ifCW != 1:
    for i in np.linspace(1,50,1):
        # Test for Pulsed ODMR
        start = 2.85e9; stop = 2.89e9; num_sweep_points = 51
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = -25; ifLooped = False
        laserInit_channel = 3; laserRead_channel = 3 # 532 is 3, 589 is 6

        num_loops               = int(0.5e6)
        laser_init_delay        = 0;          laser_init_duration    = 0*1000
        MWI_duration            = 36 #52 #420
        if laserRead_channel == 3:            laser_to_DAQ_delay     = 850  
        elif laserRead_channel == 6:          laser_to_DAQ_delay     = 1150 
        laser_to_MWI_delay      = laser_to_DAQ_delay + 150
        read_duration           = 200;        DAQ_to_laser_off_delay = 1000

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':    trackingSettings,}

        start = time.time()
        ODMRObject = ODMR(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        ODMRObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = ODMRObject.getDataFilename()
        if not ifLooped: dataReader.readData(dataFilename)
        ODMRObject.close()

        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-25/#007_ODMR_15-31-25/ODMRObject_sig_set.dat'
        # dataReader.readData(dataFilename)

else:
    # Test for CW ODMR
    for i in range(1):
        start = 2.83e9; stop = 2.91e9; num_sweep_points = 101
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = -20
        laserInit_channel = 3; laserRead_channel = 3 # 532 is 3, 589 is 6

        num_loops = int(0.2e6); wait_btwn_sig_ref = 1e3
        MWI_delay = 3e3; MWI_duration = 1e3
        laser_delay = 10; MWI_off_to_read_signal_off = 0

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                    'MWI_delay':MWI_delay,                                   'MWI_duration':MWI_duration, 
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'MWI_off_to_read_signal_off':MWI_off_to_read_signal_off,
                    'wait_btwn_sig_ref':wait_btwn_sig_ref,                   'laser_delay':laser_delay,
                    'trackingSettings': trackingSettings
                    }

        start = time.time()
        ODMRObject = ODMR_CW(settings=settings, ifPlotPulse=1) # implemented as Instrument. Turn of plot if pulse length > 1e5 ns
        ODMRObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = ODMRObject.getDataFilename()
        guess=(-2e6, 2.87e9, 0.02e9, 1)
        dataReader.readData(dataFilename, type='ODMR', ifFit=1, guess=guess);
        ODMRObject.close()

        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-17/#036_ODMR_22-34-35/ODMRObject_sig_set.dat'
        # dataReader.readData(dataFilename)






