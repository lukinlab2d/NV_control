"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *

from B00_codes.ODMRAWGFast import *  
from B00_codes.ODMR_CWAWG import *
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
        # Test for Pulsed ODMR Fast
        a1 = 2834e6; b1 = 2854e6; num_sweep_points = 21
        f1 = np.linspace(a1, b1, num_sweep_points)
        a2 = 2912e6; b2 = 2932e6; num_sweep_points = 21
        f2 = np.linspace(a2, b2, num_sweep_points)
        f3 = np.linspace(2795e6,a1-3e6,4); f4 = np.linspace(b2+3e6,2975e6,4)
        f5 = np.linspace(b1+3e6,a2-3e6, 4)
        freqsArray = np.concatenate((f3,f1,f5,f2,f4))

        # freqsArray = np.linspace(2800e6,2960e6,89)

        SDGnum=1; SRSnum=1; uwPower = -2; ifLooped = 1#(reps != 1)
        laserInit_channel=3; laserRead_channel=3; AWG_channel=18

        pi_increment = 0
        num_loops                    = int(2e5);          pitime                 = 78
        laser_init_delay             = 0;                 laser_init_duration    = 0
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}   
        AWG_output_delay             = 1450;              AWGbuffer              = 1
        read_duration                = 300;               DAQ_to_laser_off_delay = 400
        padding                      = 900+pi_increment;  MW_to_DAQ_delay        = 0
        padding_green1               = 100;               AWG_delay              = 1800

        if True:
            settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 
                        'SRSnum':SRSnum, 'uwPower':uwPower, 'SDGnum':SDGnum,
                        'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                        'pitime':         pitime, 'MW_to_DAQ_delay':MW_to_DAQ_delay,
                        'laser_to_DAQ_delay':     laser_to_DAQ_delay_directory.get(laserRead_channel, 0),    'read_duration':       read_duration,
                        'AWG_output_delay':       AWG_output_delay,      'AWGbuffer': AWGbuffer, 'AWG_delay': AWG_delay,
                        'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                        'AWG_channel':    AWG_channel, 'padding':padding,'padding_green1':padding_green1,
                        'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':    trackingSettings}

            start = time.time()
            ODMRAWGFastObject = ODMRAWGFast(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
            ODMRAWGFastObject.runScan()
            print('Total time = ' + str(time.time() - start) + ' s')

            dataFilename = ODMRAWGFastObject.getDataFilename()
            guess=(-2e6, 2.87e9, 0.02e9, 1)
            # if not ifLooped: dataReader.readData(dataFilename, type='ODMR', ifFit=0, guess=guess)
            ODMRAWGFastObject.close()

            # dataFilename = 'C:/Users/lukin2dmaterials/data/2024-04-22/#168_ODMRAWGFast_23-39-22/ODMRAWGFastObject_sig_set.dat'
            # dataReader.readData(dataFilename)

elif ifCW == 1:
    # Test for CW ODMR
    reps = 1
    for i in range(reps):
        start = 2767e6; stop = 2797e6; num_sweep_points = 41
        freqsArray = np.linspace(start, stop, num_sweep_points)
        uwPower = -8; ifLooped = 1#(reps != 1)
        laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18
        SRSnum = 1; SDGnum = 1

        num_loops = int(10e5); wait_btwn_sig_ref = 1e3; AWGbuffer = 1
        AWG_output_delay = 1445; MW_duration = 1e3
        laser_delay = 10; MW_off_to_read_signal_off = 0
        
        if True:
            settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 
                        'SRSnum':SRSnum, 'uwPower':uwPower, 'SDGnum':SDGnum, 'AWGbuffer': AWGbuffer,
                        'AWG_output_delay':AWG_output_delay, 'MW_duration':MW_duration, 
                        'laserInit_channel':laserInit_channel, 'laserRead_channel':laserRead_channel,'AWG_channel':AWG_channel,
                        'MW_off_to_read_signal_off':MW_off_to_read_signal_off,
                        'wait_btwn_sig_ref':wait_btwn_sig_ref, 'laser_delay':laser_delay,
                        'trackingSettings': trackingSettings,
                        }

            start = time.time()
            ODMRAWGObject = ODMR_CWAWG(settings=settings, ifPlotPulse=not(ifLooped)) # implemented as Instrument. Turn of plot if pulse length > 1e5 ns
            ODMRAWGObject.runScan()
            print('Total time = ' + str(time.time() - start) + ' s')

            dataFilename = ODMRAWGObject.getDataFilename()
            guess=(-2e6, 2.87e9, 0.02e9, 1)
            if not ifLooped: dataReader.readData(dataFilename, type='ODMR', ifFit=0, guess=guess)
            ODMRAWGObject.close()

            # dataFilename = 'C:/Users/lukin2dmaterials/data/2024-04-22/#083_ODMR_CWAWG_13-01-12/ODMRAWGObject_sig_set.dat'
            # dataReader.readData(dataFilename)
        
        






