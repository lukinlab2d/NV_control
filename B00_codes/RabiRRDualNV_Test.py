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
# RabiRRDualNV
reps=int(1); ifLooped=(reps!=1); laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17
vel_vpz_target = 28.76; vel_vpz_target2 = 62.7; ifNeedVel1 = 0; ifNeedVel2 = 0; ifHiloExtra=1
ifInitWvl=0; ifInitVpz=0; ifFakeRabi=0; ifAWG=ifAWG2=ifIQ=1; ifIonInit=0
start = 4; stop = 80; num_sweep_points = 20
tausArray = np.linspace(start, stop, num_sweep_points)        
num_loops                    = int(5e4)
laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:170, 14:900}
iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
laser_init_delay             = 5e2;        laser_ionInit_duration = 5e2
laser_ionInit_to_green_delay = 5e2;        laser_init_duration    = 30e3
MW_to_read_delay             = 40
laser_to_MWI_delay           = laserSwitch_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 2e3;        read_laser_duration    = read_duration    
read_duration2               = 2e3;        read_laser_duration2   = read_duration2    
shift_btwn_2NV_MW            = 100;        shift_btwn_2NV_read    = read_duration+1.7e3
AWG_buffer                   = 10;         AWG_output_delay       = 1450
AWG_buffer2                  = 10;         AWG_output_delay2      = 1450

MWfreqs1 = np.linspace(2598.44e6, 2598.44e6,reps); MWfreqs2 = np.linspace(2788.70e6,2788.70e6,reps)
for pow in np.linspace(-10.5,-10.5,reps):
    #################### NV 1 ############################################################################
    velNum = 1; vel_current = 62.7; vel_wvl = 637.22; laserRead_channel = 5
    SRSnum = 1; MWPower = -6.2; MWI_duration = 40; MWFreq  = 2598.44e6 #2953.76e6  #-5.7 #-4.0 2598.44e6
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    SDGnum = 1; AWG_channel = 18
    laser_to_DAQ_delay = laserSwitch_delay_directory.get(laserRead_channel, 0)
    # SRSnum = 3; MWPower = 6;  MWI_duration = 40; MWFreq  = 2886.745e6 #3161.27e6   #NV D1 ms+1
    # MWI_channel = 0; MWQ_channel = 0; MWswitch_channel = 15
    ##################### NV2 ############################################################################
    velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.88; laserRead2_channel = 14
    SRSnum2 = 2; MWPower2 = -8.8; MWI_duration2 = 44; MWFreq2  = 2788.70e6 #-9.8 #-8.8
    MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11
    SDGnum2 = 2; AWG2_channel = 19
    laser_to_DAQ_delay2 = laserSwitch_delay_directory.get(laserRead2_channel, 0)   
    # SRSnum2 = 4; MWPower2 = -105; MWI_duration2 = 40; MWFreq2  = 3037.20e6   #NV D2 ms+1
    # MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16
    ######################################################################################################
    # Obsolete tracking settings
    if True:
        if_tracking = 0; threshold_repumpVpz = 9; threshold_scanVpz = 13
        if_tracking2 = 0; threshold_repumpVpz2 = 9; threshold_scanVpz2 = 13
        num_loops_track = 5e3; num_of_cavity_conditioning = 1
        # NV1
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
    ###############################################################################
    if ifFakeRabi>0:
        start = 1; stop = int(ifFakeRabi*1e5); num_sweep_points = int(ifFakeRabi*1e5)
        tausArray = np.linspace(start, stop, num_sweep_points)
        num_loops = int(7e4)

    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'laserIon_channel':       laserIon_channel,
                'laser_ionInit_duration': laser_ionInit_duration, 'laser_ionInit_to_green_delay':laser_ionInit_to_green_delay,
                'iznLaserSwitch_delay':iznLaserSwitch_delay,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':MWI_duration,'MWI_duration2':MWI_duration2,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel, 
                'laser_to_DAQ_delay2':     laser_to_DAQ_delay2,  'laserRead2_channel': laserRead2_channel,
                'hiLoMWPwr_channel':      hiLoMWPwr_channel,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum, 'ifNeedVel1': ifNeedVel1,
                'vel_current2':  vel_current2, 'vel_wvl2': vel_wvl2, 'velNum2': velNum2, 'ifNeedVel2': ifNeedVel2,
                'vel_vpz_target': vel_vpz_target, 'vel_vpz_target2': vel_vpz_target2, 
                'ifInitVpz':ifInitVpz,    'ifInitWvl': ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'MWI2_channel': MWI2_channel,  'MWQ2_channel': MWQ2_channel,  'MWswitch2_channel': MWswitch2_channel,
                'RRtrackingSettings': RRtrackingSettings, 'RRtrackingSettings2': RRtrackingSettings2,
                'shift_btwn_2NV_MW':shift_btwn_2NV_MW, 'shift_btwn_2NV_read': shift_btwn_2NV_read,
                'ifFakeRabi': ifFakeRabi, 'ifAWG': ifAWG, 'ifAWG2': ifAWG2, 'ifIonInit':ifIonInit,
                'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
                'SDGnum2': SDGnum2, 'AWG2_channel':AWG2_channel, 'AWG_buffer2':AWG_buffer2, 'AWG_output_delay2':AWG_output_delay2,
                'read_duration2':read_duration2, 'read_laser_duration2':read_laser_duration2,'ifIQ':ifIQ,
                'ifHiloExtra':ifHiloExtra,}
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
