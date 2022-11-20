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
from T2R import *
import dataReader

####################################################################################################################


for i in np.linspace(1000,2000,1):
    # T2R
    start = 100010; stop = 10; num_sweep_points = 501 # randomize
    tausArray = np.linspace(start, stop, num_sweep_points)
    uwPower = -30; uwFreq = 2.8705e9
    if True:
        print(uwFreq)

        # Test for pulsed ODMR
        num_loops               = int(1e6)
        laser_init_delay_in_ns  = 10;       laser_init_duration_in_ns = 1e3
        laser_to_AFG_delay      = 1000;     piOverTwo_time            = 50
        laser_to_DAQ_delay      = 1400;     read_duration             = 300
        DAQ_to_laser_off_delay  = 20000

        # For NV tracking
        if_tracking = 0
        xy_scan_read_time      = 5;      xy_scan_settle_time    = 0.2;  
        xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
        x_minus_range          = 0.1;    x_plus_range           = 0.1
        y_minus_range          = 0.05;   y_plus_range           = 0.05
        xy_displacement_limit  = 0.01;   num_of_scans           = 4;    tracking_period = 10

        xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
        x_minus_range          = 0.1;    x_plus_range           = 0.1
        z_minus_range          = 0.4;    z_plus_range           = 0.4
        xz_displacement_limit  = 0.05; 

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

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay_in_ns': laser_init_delay_in_ns,'laser_init_duration_in_ns': laser_init_duration_in_ns,
                    'laser_to_AFG_delay':     laser_to_AFG_delay ,   'piOverTwo_time':            piOverTwo_time,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings}
        

        start = time.time()
        T2RObject = T2R(settings=settings, ifPlotPulse=False) # this is implemented as an Instrument
        T2RObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = T2RObject.getDataFilename()
        dataReader.readData(dataFilename)
        T2RObject.close()