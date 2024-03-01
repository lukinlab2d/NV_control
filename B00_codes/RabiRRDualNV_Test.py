"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.RabiRRDualNV import *
import B00_codes.dataReader as dataReader
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 200000;  ifLooped = (reps != 1); laserInit_channel = 3; vel_vpz_target = 76.18; vel_vpz_target2 = 62.7
ifInitWvl = 0; ifInitVpz = 0; ifFakeRabi=1

# RabiRRDualNV
# start = 12; stop = 252; num_sweep_points = 61 
start = 1; stop = 1000000; num_sweep_points = 1000000
tausArray = np.linspace(start, stop, num_sweep_points)        

ifNeedVel1           = 0;            ifNeedVel2              = 0

num_loops                    = int(7e3)
laser_init_delay             = 1e2;        laser_init_duration    = 30e3
MW_to_read_delay             = 1e2
laser_to_DAQ_delay_directory = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 12300;      read_laser_duration    = 12200
shift_btwn_2NV_MW            = 0;          shift_btwn_2NV_read    = 15e3

if_tracking = 0; threshold_repumpVpz = 9; threshold_scanVpz = 13
if_tracking2 = 0; threshold_repumpVpz2 = 9; threshold_scanVpz2 = 13
num_loops_track = 5e3; num_of_cavity_conditioning = 1

for i in np.linspace(1, reps, reps):
    # NV 1
    velNum = 1; vel_current = 62.2; vel_wvl = 637.22; laserRead_channel = 5
    SRSnum = 1; MWPower = -2.7; MWI_duration = 44; MWFreq  = 2747.88e6   
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2

    # SRSnum = 3; MWPower = -6;  MWI_duration = 44; MWFreq  = 3007.65e6   #NV D1 ms+1
    # MWI_channel = 0; MWQ_channel = 0; MWswitch_channel = 15

    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)

    start = vel_vpz_target - 1.6; stop = vel_vpz_target + 1.6; num_sweep_points = 65; 
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
                            'velNum': velNum,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
                            'ifInitWvl':ifInitWvl,
                            'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}
    
    # NV2
    velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.83; laserRead2_channel = 14
    SRSnum2 = 2; MWPower2 = -1; MWI_duration2 = 44; MWFreq2  = 2838.26e6 
    MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11

    # SRSnum2 = 4; MWPower2 = 10;   MWI_duration2 = 44; MWFreq2  = 2932.8e6   #NV D2 ms+1
    # MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16

    laser_to_DAQ_delay2           = laser_to_DAQ_delay_directory.get(laserRead2_channel, 0)   

    start2 = vel_vpz_target2 - 1.6; stop2 = vel_vpz_target2 + 1.6; num_sweep_points = 65; 
    vpzArray2 = np.linspace(start2, stop2, num_sweep_points)
    RRtrackingSettings2 = {'if_tracking': if_tracking2, 'threshold_repumpVpz': threshold_repumpVpz2, 'threshold_scanVpz': threshold_scanVpz2,
                         'vpzArray': vpzArray2,    'num_loops':num_loops_track,  'MWPower':MWPower2,    'MWFreq': MWFreq2,
                            'SRSnum':   SRSnum2,
                            'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                            'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration2,
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
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':MWI_duration,'MWI_duration2':MWI_duration2,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel, 
                'laser_to_DAQ_delay2':     laser_to_DAQ_delay2,  'laserRead2_channel': laserRead2_channel,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum, 'ifNeedVel1': ifNeedVel1,
                'vel_current2':  vel_current2, 'vel_wvl2': vel_wvl2, 'velNum2': velNum2, 'ifNeedVel2': ifNeedVel2,
                'vel_vpz_target': vel_vpz_target, 'vel_vpz_target2': vel_vpz_target2, 
                'ifInitVpz':ifInitVpz,    'ifInitWvl': ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'MWI2_channel': MWI2_channel,  'MWQ2_channel': MWQ2_channel,  'MWswitch2_channel': MWswitch2_channel,
                'RRtrackingSettings': RRtrackingSettings, 'RRtrackingSettings2': RRtrackingSettings2,
                'shift_btwn_2NV_MW':shift_btwn_2NV_MW, 'shift_btwn_2NV_read': shift_btwn_2NV_read,
                'ifFakeRabi': ifFakeRabi,}
    # AG = AG33522A()
    # AG.disable_PM()
    # AG.disable_RFOutput()


    start = time.time()
    RabiRRDualNVObject = RabiRRDualNV(settings=settings, ifPlotPulse=not(ifLooped)) 
    RabiRRDualNVObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if RabiRRDualNVObject.hasTracked1:
        vel_vpz_target = RabiRRDualNVObject.vpz + 0.1
    if RabiRRDualNVObject.hasTracked2:
        vel_vpz_target2 = RabiRRDualNVObject.vpz2 + 0.1
    RabiRRDualNVObject.close()
