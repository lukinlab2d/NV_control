"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
from PlotPulse import *
from CalibrateLaser2MWDelay import *
import dataReader

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
    # Blink laser
    tausArray = np.ones(300)*10
    
    if True:
        num_loops        = int(1e6)
        laser_init_delay = 100;     laser_init_duration = 20
        read_duration    = 800;      

        settings = {'num_loops':num_loops,
                    'tausArray': tausArray,
                    'laser_init_delay': laser_init_delay,'laser_init_duration': laser_init_duration,
                    'read_duration':    read_duration,   'trackingSettings':    trackingSettings}
        
        start = time.time()
        CalibrateLaser2MWDelayObject = CalibrateLaser2MWDelay(settings=settings, ifPlotPulse=True, 
                                                              ifRandom=False, ifSweepLaserInitDelay=False) # this is implemented as an Instrument
        CalibrateLaser2MWDelayObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = CalibrateLaser2MWDelayObject.getDataFilename()
        # dataReader.readData(dataFilename)
        CalibrateLaser2MWDelayObject.close()