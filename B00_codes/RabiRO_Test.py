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
from B00_codes.RabiRO import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 2;  ifLooped = (reps != 1); laserInit_channel = 7; ifInitWvl = 0
vel_vpz_start = 57;  vel_vpz_end       = 62;                   ifNeedVel = 0
vel_vpz_step  = 0.1; vel_vpz_step_time = 0.3; ifScanVpz   = 0; ifInitVpz = 0

for i in np.linspace(1, reps, reps):
    # RabiRO
    start = 12; stop = 2412; num_sweep_points = 61 
    tausArray = np.linspace(start, stop, num_sweep_points)        
    
    if True:
        velNum = 1; vel_current = 56.5; vel_wvl = 637.22
        vel_vpz_target = 69.835; laserRead_channel = 5
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.83
        vel_vpz_target          = 75; laserRead_channel = 14
    if False: 
        SRSnum = 1; MWPower = -20; MWI_duration = 68; MWFreq  = 2747.88e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    else:
        SRSnum = 2; MWPower = -42; MWI_duration = 640; MWFreq  = 3007.6e6   #NV D2, 2nd MW path
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11
    
    num_loops                    = int(1e5)
    laser_init_delay             = 1e2;        laser_init_duration    = 8e3
    MW_to_read_delay             = 1e2
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 300;        read_laser_duration    = 200

    if_tracking = 0; threshold_repumpVpz = 18; threshold_scanVpz = 23
    num_loops_track = 5e3; num_of_cavity_conditioning = 1 # if already cond' for many times, use more avgs. Not that longer
    start = vel_vpz_target - 1.6; stop = vel_vpz_target + 1.6; num_sweep_points = 65; 
    vpzArray = np.linspace(start, stop, num_sweep_points)

    ROtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
                          'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                            'SRSnum':   SRSnum,
                            'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                            'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
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
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_start': vel_vpz_start, 'vel_vpz_target': vel_vpz_target, 'vel_vpz_end': vel_vpz_end,
                'vel_vpz_step': vel_vpz_step, 'vel_vpz_step_time': vel_vpz_step_time, 'ifScanVpz': ifScanVpz, 'ifInitVpz':ifInitVpz,
                'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'ROtrackingSettings': ROtrackingSettings, }

    start = time.time()
    RabiROObject = RabiRO(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    RabiROObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if RabiROObject.hasTracked:
        vel_vpz_target = RabiROObject.vpz






    # if not ifLooped: 
    #     dataFilename = RabiROObject.getDataFilename()
    #     guess=(0.2, 100, 0, 0.9, 600)
    #     dataReader.readData(dataFilename, type='RabiDecay', guess=guess, ifFit=1)
    # RabiROObject.close()
    

