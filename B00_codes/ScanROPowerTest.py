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
import matplotlib as mpl
import matplotlib.pyplot as plt
from B00_codes.ScanROPower import *  
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

reps = 1; ifLooped = (reps != 1); ifFit = 0
laserInit_channel = 7; laserRead_channel = 5 # 532 is 3, 589 is 6, Vel is 5
for i in np.linspace(1,reps,reps):
    # Test for ScanROPower
    start = 20; stop = 65; num_sweep_points = 91
    velCurrentArray = np.linspace(start, stop, num_sweep_points)
    vel_vpz = 7.7;     vel_wvl = 637.23

    if False: 
        SRSnum = 1; MWPower = -20; MWI_duration = 72; MWFreq  = 2747.88e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    else:
        SRSnum = 2; MWPower = -17; MWI_duration = 52; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11
    
    num_loops                    = int(0.5e4)
    laser_init_delay             = 1e3;        laser_init_duration    = 8e3
    MW_to_read_delay             = 1e2
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1750}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 300;        read_laser_duration = 200

    settings = {'velCurrentArray': velCurrentArray,    'num_loops':num_loops,  'MWPower':MWPower,    'MWFreq': MWFreq,
                'SRSnum':   SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                'read_laser_duration':    read_laser_duration,   'trackingSettings':    trackingSettings,
                'MW_to_read_delay':       MW_to_read_delay,   
                'vel_vpz':                vel_vpz,               'vel_wvl':             vel_wvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}

    start = time.time()
    ScanROPowerObject = ScanROPower(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    ScanROPowerObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = ScanROPowerObject.getDataFilename()
    guess=(-2e6, 2.87e9, 0.02e9, 1)
    # if not ifLooped: dataReader.readData(dataFilename, type='ScanROPower', ifFit=ifFit, guess=guess)
    ScanROPowerObject.close()

    # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-25/#007_ScanROPower_15-31-25/ScanROPowerObject_sig_set.dat'
    # dataReader.readData(dataFilename)

