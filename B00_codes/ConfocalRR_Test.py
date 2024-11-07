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
reps=1;  ifLooped=(reps!=1); laserInit_channel = 3; ifInitVpz=0; ifInitWvl=0
ifNeedVel1=0; ifNeedVel2=0; vel_vpz_target=50; vel_vpz_target2=70
ifNeedSRS=0; ifAWG=0

# ConfocalRR
xArray = np.linspace(0.05, 0.2, 31); yArray = np.linspace(-0.32, -0.17, 31)  

num_loops                    = int(5e4);   settleTime             = 4e6
laser_init_delay             = 1e2;        laser_init_duration    = 10e3
MW_to_read_delay             = 1e2
laser_to_DAQ_delay_directory = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 2e3;        read_laser_duration    = read_duration
shift_btwn_2NV_MW            = 0;          shift_btwn_2NV_read    = read_duration+1.7e3
AWG_buffer                   = 10;         AWG_output_delay       = 1450

for i in np.linspace(1, reps, reps):
    # NV 1
    velNum = 1; vel_current = 62.7; vel_wvl = 637.20; laserRead_channel = 5
    laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)

    SRSnum = 1; MWPower = -7.1; MWI_duration = 48; MWFreq  = 2598.1e6      
    SDGnum = 1; AWG_channel = 18
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    # SRSnum = 3; MWPower = -6;  MWI_duration = 44; MWFreq  = 3007.65e6   #NV D1 ms+1
    # MWI_channel = 0; MWQ_channel = 0; MWswitch_channel = 15
    
    # NV2
    velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.88; laserRead2_channel = 14
    laser_to_DAQ_delay2 = laser_to_DAQ_delay_directory.get(laserRead2_channel, 0)

    SRSnum2 = 2; MWPower2 = -11.4; MWI_duration2 = 48; MWFreq2  = 2789.2e6 
    SDGnum2 = 2; AWG2_channel = 19
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
                'ifNeedSRS': ifNeedSRS,
                'ifAWG': ifAWG, 'SDGnum': SDGnum, 'AWG_channel':AWG_channel, 'AWG_buffer':AWG_buffer, 'AWG_output_delay':AWG_output_delay}

    start = time.time()
    ConfocalRRObject = ConfocalRR(settings=settings, ifPlotPulse=not(ifLooped)) 
    ConfocalRRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    ConfocalRRObject.close()
