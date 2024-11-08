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
from B00_codes.PlotPulse import *
from B00_codes.XY8 import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################

for i in np.linspace(1,5,5):
    # XY8
    ifRandomized = 0; ifLooped = True; normalized_style = Q_FINAL; ifStartInY = 0
    uwPower = -25; uwFreq = 2.871e9

    start = 20; stop = 20000; num_sweep_points = 31 # tau must be divisible by 4
    tausArray = np.linspace(start, stop, num_sweep_points)    

    taus1 = np.linspace(20,4220,71)
    taus2 = np.linspace(4400,20000,79)
    tausArray = np.concatenate((taus1, taus2))
    if True:
        print(uwFreq)

        # Test for pulsed ODMR
        num_loops               = int(1e6)
        laser_init_delay        = 0;        laser_init_duration       = 0
        laser_to_MWI_delay     = 1000;     piOverTwo_time            = 24
        laser_to_DAQ_delay      = 850;      read_duration             = 200
        DAQ_to_laser_off_delay  = 2500;     MWI_to_switch_delay       = 10 # cannot be between 0 and 10

        # For NV tracking
        if_tracking = 1
        xy_scan_read_time      = 5;      xy_scan_settle_time    = 1;  
        xy_scan_resolution_hor = 40;     xy_scan_resolution_ver = 20
        x_minus_range          = 0.07;   x_plus_range           = 0.05
        y_minus_range          = 0.05;   y_plus_range           = 0.05
        xy_displacement_limit  = 0.05;   num_of_scans           = 5;    tracking_period = 10

        xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
        x_minus_range          = 0.05;   x_plus_range           = 0.05
        z_minus_range          = 0.3;    z_plus_range           = 0.3
        xz_displacement_limit  = 0.25; 

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

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'piOverTwo_time':            piOverTwo_time,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
                    'MWI_to_switch_delay':    MWI_to_switch_delay,   'ifRandomized':              ifRandomized,
                    'normalized_style':       normalized_style,      'ifStartInY':                ifStartInY}
        

        start = time.time()
        XY8Object = XY8(settings=settings, ifPlotPulse=not(ifLooped), tausArray=tausArray) # this is implemented as an Instrument
        XY8Object.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = XY8Object.getDataFilename()
        if not ifLooped: dataReader.readData(dataFilename, typeNorm = normalized_style, type='XY8')
        XY8Object.close()