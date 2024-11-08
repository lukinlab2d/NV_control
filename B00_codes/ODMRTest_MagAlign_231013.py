"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import os
import numpy as np
from nidaqmx.constants import *
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
import matplotlib as mpl
import matplotlib.pyplot as plt
from B00_codes.ODMR import *  
from B00_codes.ODMR_CW import *
from B00_codes.PlotPulse import *
from B00_codes.Confocal import *
import B00_codes.dataReader as dataReader

from Newport.AGPiezo import *

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

ifCW = 1

if ifCW == 0:
    for i in np.linspace(1,50,50):
        # Test for Pulsed ODMR
        start = 2.855e9; stop = 2.885e9; num_sweep_points = 76
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = -35; ifLooped = 1
        laserInit_channel = 7; laserRead_channel = 7 # 532 is 3, 589 is 6

        num_loops                    = int(1e6)
        laser_init_delay             = 0;          laser_init_duration    = 0*1000
        MWI_duration                 = 900
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
        laser_to_MWI_delay           = laser_to_DAQ_delay + 150
        read_duration                = 250;        DAQ_to_laser_off_delay = 1000

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':    trackingSettings,}

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
    end = 7
    for rep in range(5):
        for x in range(end):
            for y in range(end):
                print('x = ' + str(x)); print('y = ' + str(y))
                start = 2400e6; stop = 3340e6; num_sweep_points = 236
                freqsArray = np.linspace(start, stop, num_sweep_points)
                uwPower = -15; ifLooped = 1
                laserInit_channel = 7; laserRead_channel = 7 # 532 is 3, 589 is 6, 2nd 532 is 7

                num_loops = int(200e3); wait_btwn_sig_ref = 1e3
                MWI_delay = 3e3; MWI_duration = 1e3
                laser_delay = 10; MWI_off_to_read_signal_off = 0

                settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower,
                            'MWI_delay':MWI_delay,                                   'MWI_duration':MWI_duration, 
                            'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                            'MWI_off_to_read_signal_off':MWI_off_to_read_signal_off,
                            'wait_btwn_sig_ref':wait_btwn_sig_ref,                   'laser_delay':laser_delay,
                            'trackingSettings': trackingSettings
                            }

                for i in range(2):
                    start = time.time()
                    ODMRObject = ODMR_CW(settings=settings, ifPlotPulse=0) # implemented as Instrument. Turn of plot if pulse length > 1e5 ns
                    ODMRObject.runScan()
                    print('Total time = ' + str(time.time() - start) + ' s')

                    dataFilename = ODMRObject.getDataFilename()
                    guess=(-2e6, 2.87e9, 0.02e9, 1)
                    if not ifLooped: dataReader.readData(dataFilename, type='ODMR', ifFit=1, guess=guess)
                    ODMRObject.close()
                
                AGPiezoObject = AGPiezo(COMchannel='COM7')
                ax = 1
                AGPiezoObject.SetStepAmplitudeNegative(ax, 50); AGPiezoObject.SetStepAmplitudePositive(ax, 50)
                if y != end-1: 
                    AGPiezoObject.RelativeMove(ax, -10000)
                    AGPiezoObject.Wait(20)
                else:
                    AGPiezoObject.RelativeMove(ax, 65000)
                    AGPiezoObject.Wait(120)
                AGPiezoObject.Close()

            AGPiezoObject = AGPiezo(COMchannel='COM7')
            ax = 2
            AGPiezoObject.SetStepAmplitudeNegative(ax, 50); AGPiezoObject.SetStepAmplitudePositive(ax, 50)
            if x != end-1: 
                AGPiezoObject.RelativeMove(ax, -10000)
                AGPiezoObject.Wait(20)
            else:
                AGPiezoObject.RelativeMove(ax, 65000)
                AGPiezoObject.Wait(120)
            AGPiezoObject.Close()

        laserTrack_channel     = 7;       if_tracking = 1
        xy_scan_read_time      = 50;      xy_scan_settle_time    = 30;  
        xy_scan_resolution_hor = 20;      xy_scan_resolution_ver = 20
        x_minus_range          = 0.01;    x_plus_range           = 0.01
        y_minus_range          = 0.01;    y_plus_range           = 0.01
        xy_displacement_limit  = 0.01;    num_of_scans           = 2
        time_btwn_trackings    = 1*60

        xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 25
        x_minus_range          = 0.02;   x_plus_range           = 0.02
        z_minus_range          = 0.1;    z_plus_range           = 0.1
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


        cfcObject = Confocal(settings=trackingSettings, laserChannel=laserTrack_channel)
        x1, y1, z = cfcObject.optimize_xy(direction=1)
        x2, y2, z = cfcObject.optimize_xy(direction=-1)
        cfcObject.set_coordinate_fnc((x1+x2)/2, (y1+y2)/2, z)
        time.sleep(1)
        cfcObject.close()  

# dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-17/#036_ODMR_22-34-35/ODMRObject_sig_set.dat'
# dataReader.readData(dataFilename)






