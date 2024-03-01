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
from T2E import *
from XY8 import *
import dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################

# For NV tracking
if_tracking = 1
xy_scan_read_time      = 4;      xy_scan_settle_time    = 1;  
xy_scan_resolution_hor = 40;     xy_scan_resolution_ver = 20
x_minus_range          = 0.07;   x_plus_range           = 0.05
y_minus_range          = 0.05;   y_plus_range           = 0.05
xy_displacement_limit  = 0.05;   num_of_scans           = 3;    tracking_period = 15

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

# Params
num_loops               = int(0.7e6)
laser_init_delay        = 0;        laser_init_duration       = 0
laser_to_MWI_delay      = 1000;     piOverTwo_time            = 24
laser_to_DAQ_delay      = 300;      read_duration             = 200
DAQ_to_laser_off_delay  = 1000;     MWI_to_switch_delay       = 10 # cannot be between 0 and 10

# for pigtail green test
for i in np.linspace(200,400,6):
    for j in np.linspace(600,2000,8):

        # T2E
        ifRandomized = 0; ifLooped = True; normalized_style = Q_FINAL
        uwPower = -25; uwFreq = 2.870e9

        # for pigtail green test
        laser_to_DAQ_delay      = i;      read_duration             = 200
        DAQ_to_laser_off_delay  = j;        MWI_to_switch_delay       = 10 # cannot be between 0 and 10

        taus1 = np.linspace(20,2020,26)
        taus2 = np.linspace(2200,9800,20)
        taus3 = np.linspace(10000,30000,11)
        taus4 = np.linspace(32000,48000,5)
        # tausArray = taus1
        tausArray = np.concatenate((taus1, taus2, taus3, taus4))
        start = tausArray[0]; stop = tausArray[-1]; num_sweep_points = len(tausArray)

        if True:
            print(uwFreq)
            settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
            'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
            'laser_to_MWI_delay':     laser_to_MWI_delay ,   'piOverTwo_time':            piOverTwo_time,
            'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
            'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
            'MWI_to_switch_delay':    MWI_to_switch_delay,   'ifRandomized':              ifRandomized,
            'normalized_style':     normalized_style}

            start = time.time()
            T2EObject = T2E(settings=settings, ifPlotPulse=not(ifLooped), tausArray=tausArray) # this is implemented as an Instrument
            T2EObject.runScan()
            print('Total time = ' + str(time.time() - start) + ' s')

            dataFilename = T2EObject.getDataFilename()
            if not ifLooped: dataReader.readData(dataFilename, typeNorm = normalized_style)
            T2EObject.close()




    # # XY8
    # ifRandomized = 0; ifLooped = True; normalized_style = Q_FINAL; ifStartInY = 0
    # uwPower = -25; uwFreq = 2.870e9

    # taus1 = np.linspace(20,2960,50)
    # taus2 = np.linspace(3000,14000,12)
    # tausArray = np.concatenate((taus1, taus2))
    # start = tausArray[0]; stop = tausArray[-1]; num_sweep_points = len(tausArray)

    # if True:
        # print(uwFreq)
        # settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
        #             'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
        #             'laser_to_MWI_delay':     laser_to_MWI_delay ,   'piOverTwo_time':            piOverTwo_time,
        #             'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
        #             'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
        #             'MWI_to_switch_delay':    MWI_to_switch_delay,   'ifRandomized':              ifRandomized,
        #             'normalized_style':       normalized_style,      'ifStartInY':                ifStartInY}
        

        # start = time.time()
        # XY8Object = XY8(settings=settings, ifPlotPulse=not(ifLooped), tausArray=tausArray) # this is implemented as an Instrument
        # XY8Object.runScan()
        # print('Total time = ' + str(time.time() - start) + ' s')

        # dataFilename = XY8Object.getDataFilename()
        # if not ifLooped: dataReader.readData(dataFilename, typeNorm = normalized_style, type='XY8')
        # XY8Object.close()