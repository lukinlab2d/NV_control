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
from B00_codes.T2ERODualNV import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 100000;  ifLooped = (reps != 1); laserInit_channel = 7; normalized_style = Q_FINAL
ifInitWvl = 0; ifInitVpz = 0; ifRandomized = 1

# T2ERODualNV
taus1 = np.linspace(80, 5e5+80, 51); tausArray = taus1       

vel_vpz_target = 76.18; vel_vpz_target2 = 64.1
pi_half        = 26;    pi_half2        = pi_half # PB isn't happy if 0 < pi_half2-pi_half < 10
ifNeedVel1     = 0;     ifNeedVel2      = 0

num_loops                    = int(1e5)
laser_init_delay             = 1e2;        laser_init_duration    = 8e3
MW_to_read_delay             = 1e2;        MWI_to_switch_delay    = 30
laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 300;        read_laser_duration    = 200
shift_btwn_2NV_MW            = 0;          shift_btwn_2NV_read    = 0

if_tracking = 0; threshold_repumpVpz = 7; threshold_scanVpz = 9
if_tracking2= 0; threshold_repumpVpz2= 8; threshold_scanVpz2= 12
num_loops_track = 5e3; num_of_cavity_conditioning = 1

for i in np.linspace(1, reps, reps):
    # NV 1
    velNum = 1; vel_current = 56.5; vel_wvl = 637.22;  laserRead_channel = 5
    SRSnum = 1; MWPower = -18; MWFreq  = 2747.88e6   #NV D1
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2

    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
    start = vel_vpz_target - 1.6; stop = vel_vpz_target + 1.6; num_sweep_points = 65; 
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
    
    # NV2
    velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.83; laserRead2_channel = 14
    SRSnum2 = 2; MWPower2 = -16.5; MWFreq2  = 2838.26e6   #NV D2, 2nd MW path
    MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11

    laser_to_DAQ_delay2           = laser_to_DAQ_delay_directory.get(laserRead2_channel, 0)   
    start2 = vel_vpz_target2 - 1.6; stop2 = vel_vpz_target2 + 1.6; num_sweep_points = 65; 
    vpzArray2 = np.linspace(start2, stop2, num_sweep_points)

    ROtrackingSettings2 = {'if_tracking': if_tracking2, 'threshold_repumpVpz': threshold_repumpVpz2, 'threshold_scanVpz': threshold_scanVpz2,
                         'vpzArray': vpzArray2,    'num_loops':num_loops_track,  'MWPower':MWPower2,    'MWFreq': MWFreq2,
                            'SRSnum':   SRSnum2,
                            'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                            'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        2*pi_half2,
                            'laser_to_DAQ_delay':     laser_to_DAQ_delay2 ,   'read_duration':       read_duration,
                            'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead2_channel,
                            'read_laser_duration':    read_laser_duration,   
                            'MW_to_read_delay':       MW_to_read_delay,   
                            'vel_current':            vel_current2,           'vel_wvl':             vel_wvl2,
                            'velNum': velNum2,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
                            'ifInitWvl':ifInitWvl,
                            'MWI_channel': MWI2_channel,  'MWQ_channel': MWQ2_channel,  'MWswitch_channel': MWswitch2_channel,}

    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay,    'MWI_to_switch_delay':       MWI_to_switch_delay,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel, 
                'laser_to_DAQ_delay2':    laser_to_DAQ_delay2,   'laserRead2_channel':        laserRead2_channel,
                'normalized_style':       normalized_style,      'pi_half': pi_half, 'pi_half2': pi_half2,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum, 'ifNeedVel1': ifNeedVel1,
                'vel_current2':  vel_current2, 'vel_wvl2': vel_wvl2, 'velNum2': velNum2, 'ifNeedVel2': ifNeedVel2,
                'vel_vpz_target': vel_vpz_target, 'vel_vpz_target2': vel_vpz_target2, 
                'ifInitVpz':ifInitVpz,    'ifInitWvl': ifInitWvl, 'ifRandomized': ifRandomized,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'MWI2_channel': MWI2_channel,  'MWQ2_channel': MWQ2_channel,  'MWswitch2_channel': MWswitch2_channel,
                'ROtrackingSettings': ROtrackingSettings, 'ROtrackingSettings2': ROtrackingSettings2,
                'shift_btwn_2NV_MW':shift_btwn_2NV_MW, 'shift_btwn_2NV_read': shift_btwn_2NV_read}

    start = time.time()
    T2ERODualNVObject = T2ERODualNV(settings=settings, ifPlotPulse=not(ifLooped)) 
    T2ERODualNVObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if T2ERODualNVObject.hasTracked1:
        vel_vpz_target = T2ERODualNVObject.vpz + 0.1
    if T2ERODualNVObject.hasTracked2:
        vel_vpz_target2 = T2ERODualNVObject.vpz2 + 0.1
    T2ERODualNVObject.close()
