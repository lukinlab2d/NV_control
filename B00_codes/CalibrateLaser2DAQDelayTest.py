"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.CalibrateLaser2DAQDelay import *

####################################################################################################################
 # For NV tracking
if_tracking = 0
xy_scan_read_time      = 5;      xy_scan_settle_time    = 0.2;  
xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
x_minus_range          = 0.1;    x_plus_range           = 0.1
y_minus_range          = 0.05;   y_plus_range           = 0.05
xy_displacement_limit  = 0.01;   num_of_scans           = 4;    tracking_period = 10

xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
x_minus_range          = 0.1;    x_plus_range           = 0.1
z_minus_range          = 0.4;    z_plus_range           = 0.4
xz_displacement_limit  = 0.05; 

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


for i in np.linspace(1000,2000,1):
    # CalibrateLaser2DAQDelay
    start = -100; stop = 500; num_sweep_points = 61
    tausArray = np.linspace(start, stop, num_sweep_points)

    laserRead_channel = 10 # 532 is 3, 589 is 6
    
    if True:
        num_loops        = int(0.5e6)
        laser_init_delay = 5e3;     laser_init_duration = 300
        read_duration    = 50;      #laser_to_DAQ_delay = 260

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops,
                    'tausArray': tausArray,
                    'laser_init_delay': laser_init_delay,'laser_init_duration': laser_init_duration,
                    'laserRead_channel': laserRead_channel,
                    'read_duration':    read_duration,   'trackingSettings':    trackingSettings}
        
        start = time.time()
        CalibrateLaser2DAQDelayObject = CalibrateLaser2DAQDelay(settings=settings, ifPlotPulse=True, 
                                                              ifRandom=False, ifSweepLaserInitDelay=False) # this is implemented as an Instrument
        CalibrateLaser2DAQDelayObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        CalibrateLaser2DAQDelayObject.close()