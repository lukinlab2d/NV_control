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
from B00_codes.XY8AWG import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################

mode = 0
reps = 1
ifLooped = False
for i in range(reps):
    tausArray = np.linspace(100, 200, 2)

    # Params
    laserInit_channel            = 7;          laserRead_channel = 7 # 532 is 3, 589 is 6
    num_loops                    = int(1e6)
    laser_init_delay             = 0;          laser_init_duration       = 0
    pi2time                      = 20;         pitime = 40; read_duration  = 250; numxy8 = 1
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
    laser_to_AWG_delay           = 1000
    AWG_output_delay = 1450; AWG_channel = 18; SRSnum = 1; SDGnum = 1; AWG_buffer = 1
    DAQ_to_laser_off_delay       = 1000

    ifRandomized = 0; normalized_style = Q_FINAL
    uwPower = -15; uwFreq = 2824.6e6

    
    if True: 
        if_tracking = 2 # 2 is for the monty setup
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

    settings = {'num_loops':num_loops, 'tausArray': tausArray, 'uwPower':uwPower, 'uwFreq': uwFreq, 'mode':  mode,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_AWG_delay':     laser_to_AWG_delay ,   'pi2time':    pi2time, 'pitime':  pitime, 'numxy8':  numxy8,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
                'AWG_output_delay':    AWG_output_delay,   'ifRandomized':              ifRandomized,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'AWG_channel': AWG_channel, 'SRSnum': SRSnum, 'SDGnum': SDGnum, 'AWG_buffer': AWG_buffer
                }
    
    start = time.time()
    XY8AWGObject = XY8AWG(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    XY8AWGObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = XY8AWGObject.getDataFilename()
    if not ifLooped: dataReader.readData(dataFilename, typeNorm = normalized_style)
    XY8AWGObject.close()