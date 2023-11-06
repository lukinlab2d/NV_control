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
from B00_codes.T2R import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################


for i in np.linspace(30,28,2):
    # T2R
    taus1 = np.linspace(10, 20010, 1001)
    tausArray = taus1

    # Params
    laserInit_channel            = 7;          laserRead_channel = 7 # 532 is 3, 589 is 6
    num_loops                    = int(1e6)
    laser_init_delay             = 0;          laser_init_duration       = 0
    pi_half                      = 20;         read_duration             = 250
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
    laser_to_MWI_delay           = laser_to_DAQ_delay + 150   
    DAQ_to_laser_off_delay       = 1000;       MWI_to_switch_delay       = 10 # cannot be between 0 and 10

    ifRandomized = 0; ifLooped = True; normalized_style = Q_FINAL
    uwPower = -15; uwFreq = 2824.6e6

    
    if True: 
        if_tracking = 2 # 2 is for the monty setup
        laserTrack_channel     = 7;       
        xy_scan_read_time      = 50;      xy_scan_settle_time    = 30;  
        xy_scan_resolution_hor = 20;      xy_scan_resolution_ver = 20
        x_minus_range          = 0.01;    x_plus_range           = 0.01
        y_minus_range          = 0.01;    y_plus_range           = 0.01
        xy_displacement_limit  = 0.01;    num_of_scans           = 2
        time_btwn_trackings    = 1000*60

        xz_scan_resolution_hor = 20;      xz_scan_resolution_ver = 25
        x_minus_range          = 0.02;    x_plus_range           = 0.02
        z_minus_range          = 0.1;     z_plus_range           = 0.1
        xz_displacement_limit  = 0.1; 

    trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
                        'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
                        'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                        'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
                        'xy_displacement_limit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
                        'tracking_period':        1e9,                   'if_tracking':            if_tracking,
                        'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                        'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                        'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                        'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}

    settings = {'num_loops':num_loops, 'tausArray': tausArray, 'uwPower':uwPower, 'uwFreq': uwFreq,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'pi_half':                   pi_half,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
                'MWI_to_switch_delay':    MWI_to_switch_delay,   'ifRandomized':              ifRandomized,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'normalized_style':     normalized_style}
    
    start = time.time()
    T2RObject = T2R(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    T2RObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = T2RObject.getDataFilename()
    if not ifLooped: dataReader.readData(dataFilename, typeNorm = normalized_style)
    T2RObject.close()