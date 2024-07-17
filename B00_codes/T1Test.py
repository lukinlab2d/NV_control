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
from B00_codes.T1 import *
import B00_codes.dataReader as dataReader

####################################################################################################################

# T1
laserInit_channel = 3; laserRead_channel   = 3; laserTrack_channel = 3 # 532 is 3, S589 is 6, W589 is 9
ifRandomized = 0; ifLooped = 0; ifPlotPulse=1
# tausArray = np.round(np.logspace(3,np.log10(25e6),21),-1)
tausArray = np.linspace(10,20,2)
uwPower = -20; uwFreq = 2.870e9

num_loops               = int(20e3); ifShimon                  = 0  # do ref under the same green as read
laser_init_delay        = 5e3;      laser_init_duration       = 30e3
laser_to_MWI_delay      = 5e3;      pi_time                   = 50
laser_to_DAQ_delay      = 850;      read_duration             = 200
DAQ_to_laser_off_delay  = 500;      MWI_to_switch_delay       = 10 # cannot be between 0 and 10
sig_to_ref_delay_Shimon = 7500

# For NV tracking
if_tracking = 1
xy_scan_read_time      = 6;      xy_scan_settle_time    = 4;  
xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
x_minus_range          = 0.03;   x_plus_range           = 0.03
y_minus_range          = 0.03;   y_plus_range           = 0.03
xy_displacement_limit  = 0.03;   num_of_scans           = 2;    tracking_period = 1e9 # just dummy
time_btwn_trackings    = 12*60

xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
x_minus_range          = 0.03;   x_plus_range           = 0.03
z_minus_range          = 0.1;    z_plus_range           = 0.1
xz_displacement_limit  = 0.1; 

if True:

    trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
                        'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
                        'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                        'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
                        'xy_displacement_limit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
                        'tracking_period':        tracking_period,       'if_tracking':            if_tracking,
                        'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                        'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                        'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                        'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings': time_btwn_trackings}

    settings = {'tausArray': tausArray,
                'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,     'pi_time':             pi_time,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,     'read_duration':       read_duration,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'MWI_to_switch_delay': MWI_to_switch_delay,
                'ifRandomized':           ifRandomized,            'ifShimon':            ifShimon,
                'sig_to_ref_delay_Shimon':sig_to_ref_delay_Shimon, 'trackingSettings':    trackingSettings,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'laserTrack_channel':     laserTrack_channel}

start = time.time()
T1Object = T1(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
T1Object.runScan()
print('Total time = ' + str(time.time() - start) + ' s')

dataFilename = T1Object.getDataFilename()
dataReader.readData(dataFilename)
T1Object.close()
        




