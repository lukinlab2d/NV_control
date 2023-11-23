"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
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
from B00_codes.ScanROFreq import *  
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

reps = 1; ifLooped = (reps != 1); ifFit = 0; laserInit_channel = 7; 
num_of_cavity_conditioning = 2; ifInitWvl = 0
for i in np.linspace(1,reps,reps):
    # Test for ScanROFreq
    if True: 
        SRSnum = 1; MWPower = -18; MWI_duration = 52; MWFreq  = 2747.88e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
        velNum = 1; vel_current = 50; vel_wvl = 637.22; laserRead_channel = 5
        start = 72; stop = 72.6; num_sweep_points = 61
    else:
        SRSnum = 2; MWPower = -16.5; MWI_duration = 52; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11
        velNum = 2; vel_current = 67; vel_wvl = 636.83; laserRead_channel = 14
        start = 62.4; stop = 63.4; num_sweep_points = 41
    vpzArray = np.linspace(start, stop, num_sweep_points)

    
    num_loops                    = int(3e4);   ifInitVpz = 0
    laser_init_delay             = 1e2;        laser_init_duration    = 8e3
    MW_to_read_delay             = 1e2
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 300;        read_laser_duration = 200

    settings = {'vpzArray': vpzArray,    'num_loops':num_loops,  'MWPower':MWPower,    'MWFreq': MWFreq,
                'SRSnum':   SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                'read_laser_duration':    read_laser_duration,   'trackingSettings':    trackingSettings,
                'MW_to_read_delay':       MW_to_read_delay,   
                'vel_current':            vel_current,           'vel_wvl':             vel_wvl, 'velNum': velNum,
                'ifInitVpz':ifInitVpz,   'num_of_cavity_conditioning': num_of_cavity_conditioning,
                'ifInitWvl':ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}

    start = time.time()
    ScanROFreqObject = ScanROFreq(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    ScanROFreqObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = ScanROFreqObject.getDataFilename()
    guess=(-2e6, 2.87e9, 0.02e9, 1)
    ScanROFreqObject.close()

    
    # if not ifLooped: dataReader.readData(dataFilename, type='ScanROFreq', ifFit=ifFit, guess=guess)
    # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-25/#007_ScanROFreq_15-31-25/ScanROFreqObject_sig_set.dat'
    # dataReader.readData(dataFilename)