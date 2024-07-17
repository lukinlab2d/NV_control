"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalRR import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 1;  ifLooped=(reps!=1); laserInit_channel = 3; ifNeedSRS=0; ifInitVpz=0; ifInitWvl=0
ifNeedVel1=1; ifNeedVel2=1; vel_vpz_target=50; vel_vpz_target2=70

# ConfocalRR
xArray = np.linspace(0.720,0.754,18); yArray = np.linspace(0.124,0.100,13)  

num_loops                    = int(100e3);  settleTime             = 5e9
laser_init_delay             = 1e2;         laser_init_duration    = 8e3
MW_to_read_delay             = 1e2
laser_to_DAQ_delay_directory = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 12300;      read_laser_duration    = 12200
shift_btwn_2NV_MW            = 0;          shift_btwn_2NV_read    = 15e3

for i in np.linspace(1, reps, reps):
    # NV 1
    velNum = 1; vel_current = 62.2; vel_wvl = 637.22; laserRead_channel = 5
    laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)

    SRSnum = 1; MWPower = -2.7; MWI_duration = 44; MWFreq  = 2747.88e6   
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    # SRSnum = 3; MWPower = -6;  MWI_duration = 44; MWFreq  = 3007.65e6   #NV D1 ms+1
    # MWI_channel = 0; MWQ_channel = 0; MWswitch_channel = 15
    
    # NV2
    velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.83; laserRead2_channel = 14
    laser_to_DAQ_delay2 = laser_to_DAQ_delay_directory.get(laserRead2_channel, 0)

    SRSnum2 = 2; MWPower2 = -1; MWI_duration2 = 44; MWFreq2  = 2838.26e6 
    MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11
    # SRSnum2 = 4; MWPower2 = 10;   MWI_duration2 = 44; MWFreq2  = 2932.8e6   #NV D2 ms+1
    # MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16
    
    settings = {'xArray': xArray, 'yArray':yArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2,
                'settleTime': settleTime,
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
                'shift_btwn_2NV_MW':shift_btwn_2NV_MW, 'shift_btwn_2NV_read': shift_btwn_2NV_read,
                'ifNeedSRS': ifNeedSRS}

    start = time.time()
    ConfocalRRObject = ConfocalRR(settings=settings, ifPlotPulse=not(ifLooped)) 
    ConfocalRRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    ConfocalRRObject.close()
