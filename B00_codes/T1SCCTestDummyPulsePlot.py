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
from B00_codes.T1SCC import *
import B00_codes.dataReader as dataReader

####################################################################################################################

# T1SCC
ifRandomized = 0; ifPlotPulse = True; ifMeaningfulRef=True
laserInit_channel = 3; laserRead_channel   = 6; laserTrack_channel = 3 # 532 is 3, S589 is 6, W589 is 9
laserIon_channel  = 8; laserShelve_channel = 3 
uwPower = -105; uwFreq = 2870e6

# for pulse plot
start = 200; stop = 200; num_sweep_points = 1 # sweep t_ion
tausArray = np.linspace(start, stop, num_sweep_points)

num_loops               = int(10)
laser_init_delay        = 500;      laser_init_duration       = 500
laser_to_MWI_delay      = 200;      pi_time                   = 92
MWI_to_shelve_delay     = 50;       shelve_duration           = 50 # both tweakable
shelve_to_ion_delay     = 50;       ion_duration              = 100                                      # tweakable
ion_to_laserRead_delay  = 500                                      # check the glass PL effect
if laserRead_channel == 3:          
    laserRead_to_DAQ_delay  = 850;  DAQ_duration = 1000            
else:   
    laserRead_to_DAQ_delay  = 1150; DAQ_duration = 1000           # 750 delay is for picoharp
DAQ_to_laser_off_delay  = 50

if True:
    # For NV tracking
    if_tracking = 1
    xy_scan_read_time      = 5;      xy_scan_settle_time    = 3;  
    xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
    x_minus_range          = 0.02;   x_plus_range           = 0.02
    y_minus_range          = 0.02;   y_plus_range           = 0.02
    xy_displacement_limit  = 0.02;   num_of_scans           = 2;    tracking_period = 1e9 # just dummy
    time_btwn_trackings    = 12*60

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

    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,     'pi_time':             pi_time,
                'MWI_to_shelve_delay':    MWI_to_shelve_delay,     'shelve_duration':     shelve_duration,
                'shelve_to_ion_delay':    shelve_to_ion_delay,     'ion_duration':        ion_duration,
                'ion_to_laserRead_delay': ion_to_laserRead_delay,
                'laserRead_to_DAQ_delay': laserRead_to_DAQ_delay,  'DAQ_duration':        DAQ_duration,    
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  
                'ifRandomized':           ifRandomized,            'ifMeaningfulRef':     ifMeaningfulRef,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'laserTrack_channel':     laserTrack_channel,      'laserShelve_channel': laserShelve_channel,
                'laserIon_channel':       laserIon_channel,        
                'trackingSettings':       trackingSettings}

    start = time.time()
    T1SCCObject = T1SCC(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    T1SCCObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = T1SCCObject.getDataFilename()
    T1SCCObject.close()
        
    


    # # for pulse plot
    # start = 200; stop = 200; num_sweep_points = 1 # sweep t_ion
    # tausArray = np.linspace(start, stop, num_sweep_points)
    # # tausArray = np.array((300e6,400e6))
    
    # num_loops               = int(10)
    # laser_init_delay        = 500;      laser_init_duration       = 500
    # laser_to_MWI_delay      = 200;      pi_time                   = 92
    # MWI_to_shelve_delay     = 50;       shelve_duration           = 50 # both tweakable
    # shelve_to_ion_delay     = 50                                       # tweakable
    # ion_to_laserRead_delay  = 500                                      # check the glass PL effect
    # if laserRead_channel == 3:          
    #     laserRead_to_DAQ_delay  = 850;  DAQ_duration = 1000            
    # else:   
    #     laserRead_to_DAQ_delay  = 1150; DAQ_duration = 1000           # 750 delay is for picoharp
    # DAQ_to_laser_off_delay  = 50


