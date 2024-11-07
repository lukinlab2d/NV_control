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
from B00_codes.ODMRAWG import *  
from B00_codes.ODMR_CWAWG import *
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

ifCW = 0

if ifCW == 0:
    # Test for Pulsed ODMR - sweeping magnet loc
    endx = 5; endy=5; count = 0
    start = 2550e6; stop = 3210e6; num_sweep_points = 166
    freqsArray = np.linspace(start, stop, num_sweep_points)

    SDGnum=1; SRSnum=1; uwPower = 0; ifLooped = 1#(reps != 1)
    laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

    num_loops                    = int(2e5)
    laser_init_delay             = 0;       laser_init_duration    = 0
    pitime                       = 34;      laser_to_AWG_delay     = 0
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    AWG_output_delay             = 1450;    AWGbuffer = 1
    read_duration                = 300;     DAQ_to_laser_off_delay = 400
    ifSingleGreenRead            = 0

    settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 
                'SRSnum':SRSnum, 'uwPower':uwPower, 'SDGnum':SDGnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                'AWG_output_delay':       AWG_output_delay,      'AWGbuffer': AWGbuffer,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                'AWG_channel':    AWG_channel,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':    trackingSettings,
                'ifSingleGreenRead': ifSingleGreenRead}

    for x in range(endx):
        for y in range(endy):
            print('x = ' + str(x)); print('y = ' + str(y)); count += 1

            if True:
                start = time.time()
                ODMRAWGObject = ODMRAWG(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
                ODMRAWGObject.runScan()
                print('Total time = ' + str(time.time() - start) + ' s')
                ODMRAWGObject.close()
            
            AGPiezoObject = AGPiezo(COMchannel='COM7')
            ax = 1
            AGPiezoObject.SetStepAmplitudeNegative(ax, 50); AGPiezoObject.SetStepAmplitudePositive(ax, 50)
            if y != endy-1: 
                AGPiezoObject.RelativeMove(ax, -15000)
                AGPiezoObject.Wait(20)
            else:
                if count != endx*endy:
                    AGPiezoObject.RelativeMove(ax, 65000)
                    AGPiezoObject.Wait(120)
            AGPiezoObject.Close()
        
        AGPiezoObject = AGPiezo(COMchannel='COM7')
        ax = 2
        AGPiezoObject.SetStepAmplitudeNegative(ax, 50); AGPiezoObject.SetStepAmplitudePositive(ax, 50)
        if x != endx-1: 
            AGPiezoObject.RelativeMove(ax, -15000)
            AGPiezoObject.Wait(20)
        else:
            if count != endx*endy:
                AGPiezoObject.RelativeMove(ax, 65000)
                AGPiezoObject.Wait(120)
        AGPiezoObject.Close()


else:
    # Test for CW ODMR
    end = 10; count = 0
    for rep in range(1):
        for x in range(end):
            for y in range(end):
                print('x = ' + str(x)); print('y = ' + str(y)); count += 1
                start1 = 2590e6; stop1 = 2670e6; num_sweep_points1 = 21
                f1 = np.linspace(start1, stop1, num_sweep_points1)
                start3 = 3090e6; stop3 = 3170e6; num_sweep_points3 = 21
                f3 = np.linspace(start3, stop3, num_sweep_points3)
                freqsArray = np.concatenate((f1,f3))
                
                # freqsArray = np.linspace(start, stop, num_sweep_points)
                
                laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18
                SRSnum = 1; SDGnum = 1; uwPower = -8; ifLooped = 1

                num_loops = int(2e5); wait_btwn_sig_ref = 2e3; AWGbuffer = 1
                AWG_output_delay = 1450; MW_duration = 5e2
                laser_delay = 10; MW_off_to_read_signal_off = 0

                settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 
                        'SRSnum':SRSnum, 'uwPower':uwPower, 'SDGnum':SDGnum, 'AWGbuffer': AWGbuffer,
                        'AWG_output_delay':AWG_output_delay, 'MW_duration':MW_duration, 
                        'laserInit_channel':laserInit_channel, 'laserRead_channel':laserRead_channel,'AWG_channel':AWG_channel,
                        'MW_off_to_read_signal_off':MW_off_to_read_signal_off,
                        'wait_btwn_sig_ref':wait_btwn_sig_ref, 'laser_delay':laser_delay,
                        'trackingSettings': trackingSettings,
                        }

                for i in range(1):
                    start = time.time()
                    ODMRAWGObject = ODMR_CWAWG(settings=settings, ifPlotPulse=not(ifLooped))  # implemented as Instrument. Turn of plot if pulse length > 1e5 ns
                    ODMRAWGObject.runScan()
                    print('Total time = ' + str(time.time() - start) + ' s')
                    ODMRAWGObject.close()
                
                AGPiezoObject = AGPiezo(COMchannel='COM7')
                ax = 1
                AGPiezoObject.SetStepAmplitudeNegative(ax, 50); AGPiezoObject.SetStepAmplitudePositive(ax, 50)
                if y != end-1: 
                    AGPiezoObject.RelativeMove(ax, -8000)
                    AGPiezoObject.Wait(20)
                else:
                    if count != end**2:
                        AGPiezoObject.RelativeMove(ax, 65000)
                        AGPiezoObject.Wait(120)
                AGPiezoObject.Close()

            AGPiezoObject = AGPiezo(COMchannel='COM7')
            ax = 2
            AGPiezoObject.SetStepAmplitudeNegative(ax, 50); AGPiezoObject.SetStepAmplitudePositive(ax, 50)
            if x != end-1: 
                AGPiezoObject.RelativeMove(ax, -8000)
                AGPiezoObject.Wait(20)
            else:
                if count != end**2:
                    AGPiezoObject.RelativeMove(ax, 65000)
                    AGPiezoObject.Wait(120)
            AGPiezoObject.Close()

        # if True:
        #     laserTrack_channel     = 7;       if_tracking = 0
        #     xy_scan_read_time      = 50;      xy_scan_settle_time    = 30;  
        #     xy_scan_resolution_hor = 20;      xy_scan_resolution_ver = 20
        #     x_minus_range          = 0.01;    x_plus_range           = 0.01
        #     y_minus_range          = 0.01;    y_plus_range           = 0.01
        #     xy_displacement_limit  = 0.01;    num_of_scans           = 2
        #     time_btwn_trackings    = 1*60

        #     xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 25
        #     x_minus_range          = 0.02;   x_plus_range           = 0.02
        #     z_minus_range          = 0.1;    z_plus_range           = 0.1
        #     xz_displacement_limit  = 0.1; 

        # trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
        #                     'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
        #                     'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
        #                     'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
        #                     'xy_displacement_limit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
        #                     'tracking_period':        1e9,                   'if_tracking':            if_tracking,
        #                     'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
        #                     'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
        #                     'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
        #                     'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}


        # cfcObject = Confocal(settings=trackingSettings, laserChannel=laserTrack_channel)
        # x1, y1, z = cfcObject.optimize_xy(direction=1)
        # x2, y2, z = cfcObject.optimize_xy(direction=-1)
        # cfcObject.set_coordinate_fnc((x1+x2)/2, (y1+y2)/2, z)
        # time.sleep(1)
        # cfcObject.close()  






