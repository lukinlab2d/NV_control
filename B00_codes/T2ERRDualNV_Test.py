"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T2ERRDualNV import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 100000;  ifLooped = (reps != 1); laserInit_channel = 7; normalized_style = Q_FINAL
ifInitWvl = 0; ifInitVpz = 0; ifRandomized = 1; ifNeedVel1 = 0; ifNeedVel2 = 0
if_tracking = 0; threshold_repumpVpz = -1; threshold_scanVpz = -1
if_tracking2= 0; threshold_repumpVpz2= -1; threshold_scanVpz2 = -1
num_loops_track = -1; num_of_cavity_conditioning = -1

# T2ERRDualNV
taus1 = np.linspace(100,20100,201); tausArray = taus1       

num_loops                    = int(1e5)
laser_init_delay             = 1e2;        laser_init_duration    = 8e3
MW_to_read_delay             = 1e2;        MWI_to_switch_delay    = 30
laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 300;        read_laser_duration    = 200
shift_btwn_2NV_MW            = 0;          shift_btwn_2NV_read    = 0

for i in np.linspace(1, reps, reps):
    # NV1
    velNum = 1; vel_current = 62.2; vel_wvl = 637.22; vel_vpz_target = -1; laserRead_channel = 5
    SRSnum  = 1; MWPower  = -7.6; pi_half  = 26; MWFreq   = 2747.88e6   #NV D1 ms-1
    SRSnum3 = 3; MWPower3 = -13;  pi_half3 = 26; MWFreq3  = 3007.65e6   #NV D1 ms+1
    MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch3_channel = 15

    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
    start = vel_vpz_target - 1.6; stop = vel_vpz_target + 1.6; num_sweep_points = 65; 
    vpzArray = np.linspace(start, stop, num_sweep_points)

    RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
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
    velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.83; vel_vpz_target2 = -1; laserRead2_channel = 14
    SRSnum2 = 2; MWPower2 = -6.1; pi_half2  = 52; MWFreq2  = 2838.26e6   #NV D2, ms-1
    SRSnum4 = 4; MWPower4 = 3;    pi_half4 = 52; MWFreq4  = 2932.8e6    #NV D2 ms+1
    MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11; MWswitch4_channel = 16

    laser_to_DAQ_delay2           = laser_to_DAQ_delay_directory.get(laserRead2_channel, 0)   
    start2 = vel_vpz_target2 - 1.6; stop2 = vel_vpz_target2 + 1.6; num_sweep_points = 65; 
    vpzArray2 = np.linspace(start2, stop2, num_sweep_points)

    RRtrackingSettings2 = {'if_tracking': if_tracking2, 'threshold_repumpVpz': threshold_repumpVpz2, 'threshold_scanVpz': threshold_scanVpz2,
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
                'RRtrackingSettings': RRtrackingSettings, 'RRtrackingSettings2': RRtrackingSettings2,
                'shift_btwn_2NV_MW':shift_btwn_2NV_MW, 'shift_btwn_2NV_read': shift_btwn_2NV_read}

    start = time.time()
    T2ERRDualNVObject = T2ERRDualNV(settings=settings, ifPlotPulse=not(ifLooped)) 
    T2ERRDualNVObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if T2ERRDualNVObject.hasTracked1:
        vel_vpz_target = T2ERRDualNVObject.vpz + 0.1
    if T2ERRDualNVObject.hasTracked2:
        vel_vpz_target2 = T2ERRDualNVObject.vpz2 + 0.1
    T2ERRDualNVObject.close()
