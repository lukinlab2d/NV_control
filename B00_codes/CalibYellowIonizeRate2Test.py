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
from B00_codes.CalibYellowIonizeRate2 import *
import B00_codes.dataReader as dataReader

####################################################################################################################


for i in np.linspace(1,100,100):
# for duration in np.array((50, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000)):
    # CalibYellowIonizeRate2
    ifRandomized = 0;  ifPlotPulse = False
    ifMeaningfulRef=0; ifRefBright=0; ifRefInitAgain = 0
    uwPower = -105; uwFreq = 2.870e9
    start = 1000e6; stop = 1000e6; num_sweep_points = 1
    tausArray = np.linspace(start, stop, num_sweep_points)
    
    laserInit_channel = 3; laserRead_channel = 6; laserTrack_channel = 3 # 532 is 3, 589 is 6
   
    num_loops               = int(1e3)
    laser_init_delay        = 5e3;      laser_init_duration       = 50e3
    laser_to_MWI_delay      = 5e6;      pi_time                   = 60
    if laserRead_channel == 3:          
        laser_to_DAQ_delay  = 850  
    else:   
        laser_to_DAQ_delay  = 1150      
    DAQ_to_laser_off_delay  = 5e2;      MWI_to_switch_delay       = 10 # cannot be between 0 and 10
    read_duration           = 499e3;    ref_laser_to_read_delay   = 1000e6
    delay_between_reads     = 1e3
    
    if True:
        # For NV tracking
        if_tracking = 1
        xy_scan_read_time      = 5;      xy_scan_settle_time    = 3;  
        xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
        x_minus_range          = 0.02;   x_plus_range           = 0.02
        y_minus_range          = 0.02;   y_plus_range           = 0.02
        xy_displacement_limit  = 0.02;   num_of_scans           = 2;    tracking_period = 1e9 # just dummy
        if np.mod(i,2)==0: time_btwn_trackings    = 5*60
        else: time_btwn_trackings    = 5*60

        xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
        x_minus_range          = 0.02;   x_plus_range           = 0.02
        z_minus_range          = 0.1;    z_plus_range           = 0.1
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
                            'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}
        
        num_reads = int(start/(read_duration + delay_between_reads))
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
                    'laserTrack_channel':     laserTrack_channel,      'delay_between_reads': delay_between_reads,
                    'num_reads':              num_reads}

        start = time.time()
        CalibYellowIonizeRate2Object = CalibYellowIonizeRate2(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
        CalibYellowIonizeRate2Object.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        
        dataFilename = CalibYellowIonizeRate2Object.getDataFilename()
        # if ifPlotPulse: dataReader.readData(dataFilename)
        CalibYellowIonizeRate2Object.close()

        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-23/#059_CalibYellowIonizeRate2_17-02-21/CalibYellowIonizeRate2Object_sig_set.dat'
        # dataReader.readDataFullData(dataFilename)