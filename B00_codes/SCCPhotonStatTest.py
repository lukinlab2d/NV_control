"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.SCCPhotonStat import *
import B00_codes.dataReader as dataReader

####################################################################################################################
i=0
for tr in np.array((1e6,2e6)):
    # SCCPhotonStat
    ifRandomized = 0; ifPlotPulse = False; ifMeaningfulRef=True
    laserInit_channel = 3; laserRead_channel   = 6; laserTrack_channel = 3 # 532 is 3, S589 is 6, W589 is 9
    laserIon_channel  = 8; laserShelve_channel = 3 
    uwPower = -25; uwFreq = 2870e6; sweepWhat = 'ti'

    tausArray = np.array((20,50,80,100,130,160,180,200,230,260,280)) # sweep t_i
    # tausArray =  np.array((10,20,60,100,150,200,250,300,350,400,450,500,550,650,750,900))# sweep t_sh
    
    num_loops               = int(10e3)
    laser_init_delay        = 5e3;      laser_init_duration       = 30e3
    laser_to_MWI_delay      = 5e3;      pi_time                   = 42
    MWI_to_shelve_delay     = 20;       shelve_duration           = 120   # both tweakable
    shelve_to_ion_delay     = 600;      ion_duration              = -1  # tweakable
    ion_to_laserRead_delay  = 0.5e6;    #tr                        = 5e6                                 
    if laserRead_channel == 3:          
        laserRead_to_DAQ_delay  = 850;  DAQ_duration = tr            
    else:   
        laserRead_to_DAQ_delay  = 1150; DAQ_duration = tr           # 750 delay is for picoharp
    DAQ_to_laser_off_delay  = 5e2;      i=i+1 
    
    if True:
        # For NV tracking
        if_tracking = 1
        xy_scan_read_time      = 5;       xy_scan_settle_time    = 3;  
        xy_scan_resolution_hor = 20;      xy_scan_resolution_ver = 20
        x_minus_range          = 0.02;    x_plus_range           = 0.02
        y_minus_range          = 0.02;    y_plus_range           = 0.02
        xy_displacement_limit  = 0.02;    num_of_scans           = 2;    tracking_period = 1e9 # just dummy
        if np.mod(i,4)==0: time_btwn_trackings    = 7*60
        else: time_btwn_trackings    = 7*60

        xz_scan_resolution_hor = 20;      xz_scan_resolution_ver = 20
        x_minus_range          = 0.02;    x_plus_range           = 0.02
        z_minus_range          = 0.1;     z_plus_range           = 0.1
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
                    'trackingSettings':       trackingSettings,        'sweepWhat':           sweepWhat}

        start = time.time()
        SCCPhotonStatObject = SCCPhotonStat(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
        SCCPhotonStatObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        
        dataFilename = SCCPhotonStatObject.getDataFilename()
        SCCPhotonStatObject.close()
    
    


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


