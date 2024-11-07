"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ODMR_RR import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 1;  ifLooped = (reps != 1); laserInit_channel = 3
ifNeedVel = 0; ifInitVpz = 0; ifInitWvl = 0; ifAWG = 1

for i in np.linspace(1, reps, reps):
    # Pulsed ODMR_RR
    freqsArray = np.linspace(2590e6, 2605e6, 51)        
    
    if True:
        velNum = 1; vel_current = 62.7; vel_wvl = 637.2; laserRead_channel = 5; vel_vpz_target = 69.8
        SRSnum = 1; MWPower = -30; MWI_duration = 640; MWFreq  = 2598.1e6   #NV D1
        SDGnum = 1; AWG_channel = 18
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.88; laserRead_channel = 14; vel_vpz_target = 75
        SRSnum = 2; MWPower = -32; MWI_duration = 640; MWFreq  = 2789.2e6   #NV D2
        SDGnum = 2; AWG_channel = 19
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11        
    
    num_loops                    = int(1.5e5)
    laser_init_delay             = 1e2;        laser_init_duration    = 8e3
    MW_to_read_delay             = 1e2      
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 2e3;        read_laser_duration    = read_duration
    AWG_buffer                   = 10;         AWG_output_delay       = 1450

    settings = {'freqsArray': freqsArray, 'num_loops':num_loops, 'MWPower':MWPower, 'SRSnum': SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration': MWI_duration,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_target': vel_vpz_target, 'ifInitVpz':ifInitVpz,
                'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'ifAWG': ifAWG, 'AWG_buffer':AWG_buffer, 'AWG_output_delay':AWG_output_delay,
                'SDGnum': SDGnum, 'AWG_channel':AWG_channel}

    start = time.time()
    ODMR_RRObject = ODMR_RR(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    ODMR_RRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
