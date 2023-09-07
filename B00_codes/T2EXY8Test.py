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
from B00_codes.PlotPulse import *
from B00_codes.T2E import *
from B00_codes.XY8 import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################

# For NV tracking

if True:
    # For NV tracking
    if_tracking = 1
    xy_scan_read_time      = 5;      xy_scan_settle_time    = 3;  
    xy_scan_resolution_hor = 20;     xy_scan_resolution_ver = 20
    x_minus_range          = 0.02;   x_plus_range           = 0.02
    y_minus_range          = 0.02;   y_plus_range           = 0.02
    xy_displacement_limit  = 0.02;   num_of_scans           = 2;    tracking_period = 1e9 # just dummy
    time_btwn_trackings    = 20*60

    xz_scan_resolution_hor = 20;     xz_scan_resolution_ver = 20
    x_minus_range          = 0.02;   x_plus_range           = 0.02
    z_minus_range          = 0.1;    z_plus_range           = 0.1
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
                    'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}

# Params
laserInit_channel       = 3;        laserRead_channel = 3 # 532 is 3, 589 is 6
num_loops               = int(0.3e6)
laser_init_delay        = 0;        laser_init_duration       = 0
pi_half                 = 22;       read_duration             = 200
if laserRead_channel == 3:          laser_to_DAQ_delay        = 850  
elif laserRead_channel == 6:        laser_to_DAQ_delay        = 1150 
laser_to_MWI_delay      = laser_to_DAQ_delay + 150   
DAQ_to_laser_off_delay  = 1000;     MWI_to_switch_delay       = 10 # cannot be between 0 and 10

ifRandomized = 0; ifLooped = True; normalized_style = Q_FINAL; 
uwPower = -25; uwFreq = 2870e6

for i in np.linspace(1,5,5):
    for j in np.linspace(1,300,300):
        # T2E   
        taus1 = np.linspace(20,200020,51)
        taus2 = np.linspace(210000, 490000, 15)
        tausArray = taus1
        # tausArray = np.concatenate((taus1, taus2))
        
        # taus1 = np.linspace(20,2020,26)
        # taus2 = np.linspace(2200,9800,20)
        # taus3 = np.linspace(10000,30000,11)
        # taus4 = np.linspace(32000,48000,5)

        if True:
            print(uwFreq)
            settings = {'num_loops':num_loops, 'tausArray': tausArray, 'uwPower':uwPower, 'uwFreq': uwFreq,
            'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
            'laser_to_MWI_delay':     laser_to_MWI_delay ,   'pi_half':                   pi_half,
            'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
            'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
            'MWI_to_switch_delay':    MWI_to_switch_delay,   'ifRandomized':              ifRandomized,
            'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
            'normalized_style':     normalized_style}

            start = time.time()
            T2EObject = T2E(settings=settings, ifPlotPulse=False) # this is implemented as an Instrument
            T2EObject.runScan()
            print('Total time = ' + str(time.time() - start) + ' s')

            dataFilename = T2EObject.getDataFilename()
            if not ifLooped: dataReader.readData(dataFilename, typeNorm = normalized_style)
            T2EObject.close()


    # XY8
    ifStartInY = 0

    taus1 = np.linspace(20,2960,50)
    taus2 = np.linspace(3000,14000,12)
    tausArray = np.concatenate((taus1, taus2))
    start = tausArray[0]; stop = tausArray[-1]; num_sweep_points = len(tausArray)

    if True:
        print(uwFreq)
        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'pi_half':                   pi_half,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
                    'MWI_to_switch_delay':    MWI_to_switch_delay,   'ifRandomized':              ifRandomized,
                    'normalized_style':       normalized_style,      'ifStartInY':                ifStartInY,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel}
        

        start = time.time()
        XY8Object = XY8(settings=settings, ifPlotPulse=not(ifLooped), tausArray=tausArray) # this is implemented as an Instrument
        XY8Object.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = XY8Object.getDataFilename()
        if not ifLooped: dataReader.readData(dataFilename, typeNorm = normalized_style, type='XY8')
        XY8Object.close()