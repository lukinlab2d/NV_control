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
from B00_codes.PlotPulse import *
from B00_codes.CalibYellowIonizeRate import *
import B00_codes.dataReader as dataReader

####################################################################################################################


for i in np.linspace(1,5,5):
    # CalibYellowIonizeRate
    ifRandomized = 0; ifLooped = 0; ifPlotPulse = False
    ifMeaningfulRef=True; ifRefBright=True; ifRefInitAgain = True
    uwPower = -105; uwFreq = 2.870e9
    start = 0; stop = 30000e3; num_sweep_points = 31
    tausArray = np.linspace(start, stop, num_sweep_points)
    
    laserInit_channel = 3; laserRead_channel = 6; laserTrack_channel = 3 # 532 is 3, 589 is 6

    num_loops               = int(2e3)
    laser_init_delay        = 10;     laser_init_duration       = 5e2
    laser_to_MWI_delay      = 5e3;    pi_time                   = 60
    if laserRead_channel == 3:          
        laser_to_DAQ_delay  = 850  
    else:   
        laser_to_DAQ_delay  = 1150      
    DAQ_to_laser_off_delay  = 0;        MWI_to_switch_delay       = 10 # cannot be between 0 and 10
    read_duration           = 200e3;   ref_laser_to_read_delay   = 30000e3

    if True:
        # For NV tracking
        if_tracking = 1
        xy_scan_read_time      = 7;     xy_scan_settle_time    = 4;  
        xy_scan_resolution_hor = 30;     xy_scan_resolution_ver = 30
        x_minus_range          = 0.075;  x_plus_range           = 0.075
        y_minus_range          = 0.075;  y_plus_range           = 0.075
        xy_displacement_limit  = 0.05;   num_of_scans           = 3;    tracking_period = 1e9 # just dummy
        if np.mod(i,2)==0: time_btwn_trackings    = 2*60
        else: time_btwn_trackings    = 20*60

        xz_scan_resolution_hor = 30;     xz_scan_resolution_ver = 30
        x_minus_range          = 0.075;  x_plus_range           = 0.075
        z_minus_range          = 0.2;    z_plus_range           = 0.2
        xz_displacement_limit  = 0.2; 

        trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
                            'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
                            'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                            'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
                            'xy_displacement_limit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
                            'tracking_period':        tracking_period,       'if_tracking':            if_tracking,
                            'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                            'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                            'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                            'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points,
                    'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,     'pi_time':             pi_time,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,     'read_duration':       read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'MWI_to_switch_delay': MWI_to_switch_delay,
                    'ifRandomized':           ifRandomized,            'ifRefBright':         ifRefBright,
                    'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                    'trackingSettings':       trackingSettings,        'ifMeaningfulRef':     ifMeaningfulRef,
                    'ref_laser_to_read_delay':ref_laser_to_read_delay, 'ifRefInitAgain':      ifRefInitAgain,
                    'laserTrack_channel':     laserTrack_channel,}

        start = time.time()
        CalibYellowIonizeRateObject = CalibYellowIonizeRate(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
        CalibYellowIonizeRateObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        
        dataFilename = CalibYellowIonizeRateObject.getDataFilename()
        if ifPlotPulse: dataReader.readData(dataFilename)
        CalibYellowIonizeRateObject.close()

        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-23/#059_CalibYellowIonizeRate_17-02-21/CalibYellowIonizeRateObject_sig_set.dat'
        # dataReader.readDataFullData(dataFilename)
        




