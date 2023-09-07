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
from B00_codes.PlotPulse import *
from B00_codes.CalibrateReadoutLength import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################


for i in np.linspace(1000,2000,1):
    # CalibrateReadoutLength
    start = 1010; stop = 20; num_sweep_points = 100
    ifRandomized = 0; ifLooped = 0; ifPlotPulse = True
    tausArray = np.linspace(start, stop, num_sweep_points)
    uwPower = -20; uwFreq = 2.870e9
    laserRead_channel = 3 # 532 is 3, 589 is 6
    if True:
        print(uwFreq)

        # Test for pulsed ODMR
        num_loops               = int(5e5)
        laser_init_delay        = 0;        laser_init_duration       = 0
        laser_to_MWI_delay      = 1000;     pi_time                   = 112
        if laserRead_channel == 3:          
            laser_to_DAQ_delay  = 850  
        elif laserRead_channel == 6 or laserRead_channel == 9:   
            laser_to_DAQ_delay  = 1150          
        DAQ_to_laser_off_delay  = 1000;     MWI_to_switch_delay       = 10 # cannot be between 0 and 10

        # For NV tracking
        if_tracking = 0
        xy_scan_read_time      = 5;      xy_scan_settle_time    = 0.2;  
        xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
        x_minus_range          = 0.1;    x_plus_range           = 0.1
        y_minus_range          = 0.05;   y_plus_range           = 0.05
        xy_displacement_limit  = 0.025;  num_of_scans           = 5;    tracking_period = 25

        xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
        x_minus_range          = 0.1;    x_plus_range           = 0.1
        z_minus_range          = 0.3;    z_plus_range           = 0.3
        xz_displacement_limit  = 0.15; 

        trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
                            'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
                            'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                            'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
                            'xy_displacement_liSmit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
                            'tracking_period':        tracking_period,       'if_tracking':            if_tracking,
                            'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                            'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                            'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                            'xz_displacement_limit':  xz_displacement_limit,}

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points,
                    'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,     'pi_time':             pi_time,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,     
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'MWI_to_switch_delay': MWI_to_switch_delay,
                    'ifRandomized':           ifRandomized,            'laserRead_channel': laserRead_channel,
                    'trackingSettings':    trackingSettings}

        # start = time.time()
        # CalibrateReadoutLengthObject = CalibrateReadoutLength(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
        # CalibrateReadoutLengthObject.runScan()
        # print('Total time = ' + str(time.time() - start) + ' s')
        
        # dataFilename = CalibrateReadoutLengthObject.getDataFilename()
        # dataReader.readData(dataFilename, typeNorm == REF_MINUS_SIG)
        # CalibrateReadoutLengthObject.close()

        dataFilename = 'C:/Users/lukin2dmaterials/data/2023-06-14/#011_CalibrateReadoutLength_17-11-46/CalibrateReadoutLengthObject_sig_set.dat'
        dataReader.readData(dataFilename, typeNorm=REF_MINUS_SIG)
        




