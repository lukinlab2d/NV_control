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
reps = 1;  ifLooped = (reps != 1); laserInit_channel = 3; ifAWG=1; ifInitWvl = 0
vel_vpz_start = 57;  vel_vpz_end       = 62;                       ifNeedVel = 0
vel_vpz_step  = 0.1; vel_vpz_step_time = 0.3; ifScanVpz   = 0;     ifInitVpz = 0; ifFakeRabi=0

for i in np.linspace(1, reps, reps):
    # RabiRR
    start = 2; stop = 122; num_sweep_points = 61
    tausArray = np.linspace(start, stop, num_sweep_points)        
    
    if True:
        velNum = 1; vel_current = 62.7; vel_wvl = 637.20; laserRead_channel = 5; vel_vpz_target = 28.58
    else:
        velNum = 2; vel_current = 67.0; vel_wvl = 636.88; laserRead_channel = 14;vel_vpz_target = 27.04; 
    if True: 
        SRSnum = 1; MWPower = 0.5; MWI_duration = 20; MWFreq  = 2597.9e6 #NV D1
        SDGnum = 1; AWG_channel = 18 
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    else:
        SRSnum = 2; MWPower = -2; MWI_duration = 20; MWFreq  = 2785.2e6   #NV D2, 2nd MW path
        SDGnum = 2; AWG_channel = 19
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11
    
    num_loops                    = int(1e5)
    laser_init_delay             = 1e2;        laser_init_duration    = 8e3
    MW_to_read_delay             = 1e2
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MW_delay            = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 2e3;        read_laser_duration    = read_duration
    AWG_buffer                   = 10;         AWG_output_delay       = 1450     

    if_tracking = 0; threshold_repumpVpz = 18; threshold_scanVpz = 23
    num_loops_track = 5e3; num_of_cavity_conditioning = 1 # if already cond' for many times, use more avgs. Not that longer
    start = vel_vpz_target - 1.6; stop = vel_vpz_target + 1.6; num_sweep_points = 65; 
    vpzArray = np.linspace(start, stop, num_sweep_points)

    RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
                          'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                            'SRSnum':   SRSnum,
                            'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                            'laser_to_MW_delay':     laser_to_MW_delay ,   'MWI_duration':        MWI_duration,
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
                'laser_to_MW_delay':     laser_to_MW_delay ,   
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_start': vel_vpz_start, 'vel_vpz_target': vel_vpz_target, 'vel_vpz_end': vel_vpz_end,
                'vel_vpz_step': vel_vpz_step, 'vel_vpz_step_time': vel_vpz_step_time, 'ifScanVpz': ifScanVpz, 'ifInitVpz':ifInitVpz,
                'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'RRtrackingSettings': RRtrackingSettings, 'ifFakeRabi':ifFakeRabi, 'MWI_duration':MWI_duration,
                'ifAWG': ifAWG, 'SDGnum': SDGnum, 'AWG_channel':AWG_channel, 'AWG_buffer':AWG_buffer, 'AWG_output_delay':AWG_output_delay}

    start = time.time()
    RabiRRObject = RabiRR(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    RabiRRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if RabiRRObject.hasTracked:
        vel_vpz_target = RabiRRObject.vpz


