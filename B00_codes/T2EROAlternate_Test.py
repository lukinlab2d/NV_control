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
from B00_codes.T2ERO import *
from B00_codes.XY8 import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################
reps = 10000;  ifRandomized = 1; ifLooped = (reps!=1); normalized_style = Q_FINAL
laserInit_channel = 7;  ifScanVpz   = 0;  ifInitVpz    = 1;   ifInitWvl         = 0
vel_vpz_start     = 57; vel_vpz_end = 62; vel_vpz_step = 0.1; vel_vpz_step_time = 0.3

vel_vpz_target_1        = 78.6;        vel_vpz_target_2         = 76.5

# Params
num_loops                    = int(3e3)
laser_init_delay             = 1e2;        laser_init_duration    = 8e3
MW_to_read_delay             = 1e2;        MWI_to_switch_delay    = 30
laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 300;        read_laser_duration    = 200

if_tracking = 1; threshold_repumpVpz = 9; threshold_scanVpz = 11
num_loops_track = 5e3; num_of_cavity_conditioning = 1

# T2ERO   
taus1 = np.linspace(100, 5100, 51)
tausArray = taus1

for i in np.linspace(1,reps,reps):    
    # NV1
    velNum = 1; vel_current = 56.5; vel_wvl = 637.22; laserRead_channel = 5 
    SRSnum = 1; MWPower = -20; pi_half = 34; MWFreq  = 2747.88e6   #NV D1
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2

    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
    start = vel_vpz_target_1 - 1.6; stop = vel_vpz_target_1 + 1.6; num_sweep_points = 65; 
    vpzArray = np.linspace(start, stop, num_sweep_points)

    ROtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
                'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                    'SRSnum':   SRSnum,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        2*pi_half,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'read_laser_duration':    read_laser_duration,   
                    'MW_to_read_delay':       MW_to_read_delay,   
                    'vel_current':            vel_current,           'vel_wvl':             vel_wvl,
                    'velNum': velNum,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
                    'ifInitWvl':ifInitWvl,
                    'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}
    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_to_switch_delay':       MWI_to_switch_delay,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'normalized_style':       normalized_style,       'pi_half':        pi_half,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_start': vel_vpz_start, 'vel_vpz_target': vel_vpz_target_1, 'vel_vpz_end': vel_vpz_end,
                'vel_vpz_step': vel_vpz_step, 'vel_vpz_step_time': vel_vpz_step_time, 'ifScanVpz': ifScanVpz, 
                'ifInitVpz':ifInitVpz, 'ifInitWvl':ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'ROtrackingSettings': ROtrackingSettings, 'ifRandomized': ifRandomized }

    start = time.time()
    T2EROObject = T2ERO(settings=settings, ifPlotPulse=(reps == 1)) # this is implemented as an Instrument
    T2EROObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if T2EROObject.hasTracked:
        vel_vpz_target_1 = T2EROObject.vpz + 0.1
    T2EROObject.close()

    ##############################################################################################
    # NV2
    velNum = 2; vel_current = 67; vel_wvl = 636.83; laserRead_channel = 14
    SRSnum = 2; MWPower = -17; pi_half = 26; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
    MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11

    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
    start = vel_vpz_target_2 - 1.6; stop = vel_vpz_target_2 + 1.6; num_sweep_points = 65; 
    vpzArray = np.linspace(start, stop, num_sweep_points)

    ROtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
                'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                    'SRSnum':   SRSnum,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        2*pi_half,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'read_laser_duration':    read_laser_duration,   
                    'MW_to_read_delay':       MW_to_read_delay,   
                    'vel_current':            vel_current,           'vel_wvl':             vel_wvl,
                    'velNum': velNum,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
                    'ifInitWvl':ifInitWvl,
                    'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}
    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_to_switch_delay':       MWI_to_switch_delay,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'normalized_style':       normalized_style,       'pi_half':        pi_half,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_start': vel_vpz_start, 'vel_vpz_target': vel_vpz_target_2, 'vel_vpz_end': vel_vpz_end,
                'vel_vpz_step': vel_vpz_step, 'vel_vpz_step_time': vel_vpz_step_time, 'ifScanVpz': ifScanVpz, 
                'ifInitVpz':ifInitVpz, 'ifInitWvl':ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'ROtrackingSettings': ROtrackingSettings, 'ifRandomized': ifRandomized }

    start = time.time()
    T2EROObject = T2ERO(settings=settings, ifPlotPulse=(reps == 1)) # this is implemented as an Instrument
    T2EROObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if T2EROObject.hasTracked:
        vel_vpz_target_2 = T2EROObject.vpz + 0.1
    T2EROObject.close()