"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T2ESCC import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################
i=0
for rep in np.linspace(1,20,20):
    for reps in np.linspace(1,20,20):
        # T2ESCC
        ifRandomized = 0; ifPlotPulse = 1; ifMeaningfulRef=1; normalized_style = Q_FINAL
        laserInit_channel = 3; laserRead_channel   = 6; laserTrack_channel = 3 # 532 is 3, S589 is 6, W589 is 9
        laserIon_channel  = 8; laserShelve_channel = 3 
        uwPower = -105; uwFreq = 2870e6


         # for pulse plot
        start = 200; stop = 200; num_sweep_points = 1
        tausArray = np.linspace(start, stop, num_sweep_points)
        
        num_loops               = int(10)
        laser_init_delay        = 500;      laser_init_duration       = 500
        laser_to_MWS_delay      = 200;      pi_half                   = 92
        MWS_to_shelve_delay     = 50;       shelve_duration           = 50  # both tweakable
        shelve_to_ion_delay     = 50;       ion_duration              = 100 # tweakable
        ion_to_laserRead_delay  = 500                                       # check the glass PL effect
        if laserRead_channel == 3:          
            laserRead_to_DAQ_delay  = 850;  DAQ_duration = 1000            
        else:   
            laserRead_to_DAQ_delay  = 1150; DAQ_duration = 1000           # 750 delay is for picoharp
        DAQ_to_laser_off_delay  = 50;      MWI_to_switch_delay       = 10

        if True:
            # For NV tracking
            if_tracking = 1
            xy_scan_read_time      = 5;      xy_scan_settle_time    = 3;  
            xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
            x_minus_range          = 0.02;   x_plus_range           = 0.02
            y_minus_range          = 0.02;   y_plus_range           = 0.02
            xy_displacement_limit  = 0.02;   num_of_scans           = 2;    tracking_period = 1e9 # just dummy
            if np.mod(i,4)==0: time_btwn_trackings    = 20*60
            else: time_btwn_trackings    = 20*60

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
                        'laser_to_MWS_delay':     laser_to_MWS_delay ,     'pi_half':             pi_half,
                        'MWI_to_switch_delay':    MWI_to_switch_delay,
                        'MWS_to_shelve_delay':    MWS_to_shelve_delay,     'shelve_duration':     shelve_duration,
                        'shelve_to_ion_delay':    shelve_to_ion_delay,     'ion_duration':        ion_duration,
                        'ion_to_laserRead_delay': ion_to_laserRead_delay,
                        'laserRead_to_DAQ_delay': laserRead_to_DAQ_delay,  'DAQ_duration':        DAQ_duration,    
                        'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  
                        'ifRandomized':           ifRandomized,            'ifMeaningfulRef':     ifMeaningfulRef,
                        'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                        'laserTrack_channel':     laserTrack_channel,      'laserShelve_channel': laserShelve_channel,
                        'laserIon_channel':       laserIon_channel,        
                        'trackingSettings':       trackingSettings,        'normalized_style':     normalized_style}

            start = time.time()
            T2ESCCObject = T2ESCC(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
            T2ESCCObject.runScan()
            print('Total time = ' + str(time.time() - start) + ' s')
            
            dataFilename = T2ESCCObject.getDataFilename()
            T2ESCCObject.close()
    
