"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.RabiRR import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 5;  ifLooped = (reps != 1); laserInit_channel = 7; ifInitWvl = 0
vel_vpz_start = 57;  vel_vpz_end       = 62
vel_vpz_step  = 0.1; vel_vpz_step_time = 0.3; ifScanVpz   = 0; ifInitVpz = 0

# RabiRRAlternate
start = 12; stop = 252; num_sweep_points = 61 
tausArray = np.linspace(start, stop, num_sweep_points)        

vel_vpz_target_1        = 75.5;        vel_vpz_target_2         = 75.7
ifNeedVel1              = 0;           ifNeedVel2               = 0

num_loops                    = int(1e4)
laser_init_delay             = 1e2;        laser_init_duration    = 8e3
MW_to_read_delay             = 1e2
laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1750, 14:900}
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 12300;        read_laser_duration  = 12200

if_tracking = 0; threshold_repumpVpz = 13; threshold_scanVpz = 20
num_loops_track = 5e3; num_of_cavity_conditioning = 1

for i in np.linspace(1, reps, reps):
    # NV 1
    velNum = 1; vel_current = 62.2; vel_wvl = 637.22;  laserRead_channel = 5
    SRSnum = 1; MWPower = -17.6; MWI_duration = 52; MWFreq  = 2747.88e6   #NV D1
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2

    # SRSnum = 3; MWPower = -23; MWI_duration = 56; MWFreq  = 3007.65e6  #NV D1, ms+1, 3rd MW path
    # MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 15 # IQ are dummy

    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)

    start = vel_vpz_target_1 - 1.6; stop = vel_vpz_target_1 + 1.6; num_sweep_points = 65; 
    vpzArray = np.linspace(start, stop, num_sweep_points)

    RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
                         'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                            'SRSnum':   SRSnum,
                            'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                            'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                            'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                            'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                            'read_laser_duration':    read_laser_duration,   
                            'MW_to_read_delay':       MW_to_read_delay,   
                            'vel_current':            vel_current,           'vel_wvl':             vel_wvl,
                            'velNum': velNum, 'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
                            'ifInitWvl':ifInitWvl,
                            'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}
    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_start': vel_vpz_start, 'vel_vpz_target': vel_vpz_target_1, 'vel_vpz_end': vel_vpz_end,
                'vel_vpz_step': vel_vpz_step, 'vel_vpz_step_time': vel_vpz_step_time, 'ifScanVpz': ifScanVpz, 'ifInitVpz':ifInitVpz,
                'ifInitWvl': ifInitWvl, 'ifNeedVel':ifNeedVel1,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'RRtrackingSettings': RRtrackingSettings, }

    start = time.time()
    RabiRRObject = RabiRR(settings=settings, ifPlotPulse=not(ifLooped)) 
    RabiRRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if RabiRRObject.hasTracked:
        vel_vpz_target_1 = RabiRRObject.vpz + 0.1
    RabiRRObject.close()

    # ###############################################################################################
    # # NV2
    # velNum = 2; vel_current = 67; vel_wvl = 636.83; laserRead_channel = 14
    # SRSnum = 2; MWPower = -16.5; MWI_duration = 52; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
    # MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11

    # # SRSnum = 4; MWPower = -9; MWI_duration = 52; MWFreq  = 2932.8e6  #NV D2, ms+1, 4th MW path
    # # MWI_channel = 0; MWQ_channel = 0; MWswitch_channel = 16 # IQ are dummy

    # laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   

    # start = vel_vpz_target_2 - 1.6; stop = vel_vpz_target_2 + 1.6; num_sweep_points = 65; 
    # vpzArray = np.linspace(start, stop, num_sweep_points)

    # RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
    #                     'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
    #                         'SRSnum':   SRSnum,
    #                         'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
    #                         'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
    #                         'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
    #                         'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
    #                         'read_laser_duration':    read_laser_duration,   
    #                         'MW_to_read_delay':       MW_to_read_delay,   
    #                         'vel_current':            vel_current,           'vel_wvl':             vel_wvl,
    #                         'velNum': velNum,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
    #                         'ifInitWvl':ifInitWvl,
    #                         'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}
    # settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
    #             'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
    #             'laser_to_MWI_delay':     laser_to_MWI_delay ,   
    #             'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
    #             'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
    #             'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
    #             'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
    #             'vel_vpz_start': vel_vpz_start, 'vel_vpz_target': vel_vpz_target_2, 'vel_vpz_end': vel_vpz_end,
    #             'vel_vpz_step': vel_vpz_step, 'vel_vpz_step_time': vel_vpz_step_time, 'ifScanVpz': ifScanVpz, 'ifInitVpz':ifInitVpz,
    #             'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel2,
    #             'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
    #             'RRtrackingSettings': RRtrackingSettings, }

    # start = time.time()
    # RabiRRObject = RabiRR(settings=settings, ifPlotPulse=not(ifLooped)) 
    # RabiRRObject.runScan()
    # print('Total time = ' + str(time.time() - start) + ' s')
    # if RabiRRObject.hasTracked:
    #     vel_vpz_target_2 = RabiRRObject.vpz + 0.1
    # RabiRRObject.close()