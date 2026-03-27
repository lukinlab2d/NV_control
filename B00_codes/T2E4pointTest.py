"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
import time
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T2E4point import *

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################

reps = int(30); ifLooped = (reps!=-1); srate=2.5e8
for i in np.linspace(1,reps,reps):
    # T2E4point
    # tausArray = np.linspace(20,5020,51)
    # tausArray = np.concatenate((np.linspace(20,200,10),np.round(np.logspace(np.log10(250),np.log10(4e3),12),-1)))
    # tausArray = np.round(np.logspace(np.log10(20),np.log10(60e3),26),-1)
    # tausArray = np.concatenate((np.linspace(20,156,35),np.round(np.logspace(np.log10(160),np.log10(4e3),12),-1)))
    tausArray = np.concatenate((np.linspace(20,156,35),np.round(np.logspace(np.log10(160),np.log10(10e3),14),-1)))

    # Params
    laserInit_channel            = 3;          laserRead_channel     = 3 # 532 is 3, 589 is 6
    num_loops                    = int(5e5);   phi_3rdRead           = 0
    laser_init_delay             = 0;          laser_init_duration   = 0
    pitime                       = 28;         pi2time               = pitime/2   
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
    laser_to_AWG_delay           = 5e3;        read_duration         = 300
    DAQ_to_laser_off_delay       = 5e3;        MW_to_read_delay      = 40

    AWG_output_delay = 1450; AWG_channel = 18; SRSnum = 1; SDGnum = 1; AWG_buffer = 10
    uwPower = -8; uwFreq = 2699.87e6; ifRandomized = 1

    if True: 
        if_tracking = 0 # 2 is for the monty setup
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

    settings = {'num_loops':num_loops, 'tausArray': tausArray, 'uwPower':uwPower, 'uwFreq': uwFreq,'srate':srate,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_AWG_delay':     laser_to_AWG_delay ,   'pi2time':                   pi2time,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
                'AWG_output_delay':    AWG_output_delay,   'ifRandomized':              ifRandomized,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'AWG_channel': AWG_channel, 'SRSnum': SRSnum, 'SDGnum': SDGnum, 'AWG_buffer': AWG_buffer,
                'MW_to_read_delay':MW_to_read_delay, 'pitime':pitime, 'phi_3rdRead':phi_3rdRead
                }
    
    start = time.time()
    T2E4pointObject = T2E4point(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    T2E4pointObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = T2E4pointObject.getDataFilename()
    T2E4pointObject.close()