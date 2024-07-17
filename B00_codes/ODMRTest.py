"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *

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

ifCW = 0

if ifCW == 0:
    reps = 1
    for i in range(reps):
        # Test for Pulsed ODMR
        start = 2850e6; stop = 2890e6; num_sweep_points = 51
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = -10; ifLooped = (reps != 1)
        laserInit_channel = 3; laserRead_channel = 3 # 532 is 3, 589 is 6
        SRSnum = 2; MWswitch_channel=11

        num_loops                    = int(1e5)
        laser_init_delay             = 0;          laser_init_duration    = 0*1000
        MWI_duration                 = 250
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
        laser_to_MWI_delay           = laser_to_DAQ_delay + 150
        read_duration                = 250;        DAQ_to_laser_off_delay = 25000

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':    trackingSettings,
                    'SRSnum':SRSnum, 'MWswitch_channel':MWswitch_channel,}

        start = time.time()
        ODMRObject = ODMR(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        ODMRObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = ODMRObject.getDataFilename()
        guess=(-2e6, 2.87e9, 0.02e9, 1)
        if not ifLooped: dataReader.readData(dataFilename, type='ODMR', ifFit=1, guess=guess)
        ODMRObject.close()

        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-25/#007_ODMR_15-31-25/ODMRObject_sig_set.dat'
        # dataReader.readData(dataFilename)

else:
    # Test for CW ODMR
    reps = 1
    for i in range(reps):
        start = 2850e6; stop = 2890e6; num_sweep_points = 41
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = 10; ifLooped = (reps != 1)
        laserInit_channel = 3; laserRead_channel = 3
        SRSnum = 2; MWswitch_channel=11

        num_loops = int(100e3); wait_btwn_sig_ref = 30e3
        MWI_delay = 1e3; MWI_duration = 1e3
        laser_delay = 10; MWI_off_to_read_signal_off = 0

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                    'MWI_delay':MWI_delay,                                   'MWI_duration':MWI_duration, 
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'MWI_off_to_read_signal_off':MWI_off_to_read_signal_off,
                    'wait_btwn_sig_ref':wait_btwn_sig_ref,                   'laser_delay':laser_delay,
                    'trackingSettings': trackingSettings, 'SRSnum':SRSnum, 'MWswitch_channel':MWswitch_channel,
                    }

        start = time.time()
        ODMRObject = ODMR_CW(settings=settings, ifPlotPulse=not(ifLooped)) # implemented as Instrument. Turn of plot if pulse length > 1e5 ns
        ODMRObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = ODMRObject.getDataFilename()
        guess=(-2e6, 2.87e9, 0.02e9, 1)
        if not ifLooped: dataReader.readData(dataFilename, type='ODMR', ifFit=1, guess=guess)
        ODMRObject.close()

        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-17/#036_ODMR_22-34-35/ODMRObject_sig_set.dat'
        # dataReader.readData(dataFilename)






