"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.SCCRRPhotonStatIrber import *
import B00_codes.dataReader as dataReader

####################################################################################################################
reps = -1
for zzzz in np.linspace(1,2,1):
    for iii in np.linspace(1,2,1):
        # SCCRRPhotonStatIrber
        ifRandomized = 0;  ifPlotPulse = (reps==1);  ifMWReadLowDutyCycle = 0; ifNeedVel = 0
        ifFancySpinInit = 0; nSpinInit = -1; sweepWhich = 'ti'; ifMWDuringRead = 1; ifMW2DuringRead = 1
        laserInit_channel = 3; laserIon_channel = 10; hiLoMWPwr_channel = 17; ifInitVpz   = 0; ifInitWvl = 0
        if iii==1: 
            velNum = 1; vel_current = 62.2; vel_wvl = 637.22; vel_vpz_target = 72.3
            
            SRSnum = 1;  MWPower = -2.7; pi_time  = 44; MWFreq   = 2747.88e6   #NV D1 ms-1
            MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2; laserRead_channel = 5
            
            SRSnum2 = 3; MWPower2 = -6;  pi_time2 = 40; MWFreq2  = 3007.65e6   #NV D1 ms+1
            MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 15

            # SRSnum = 3;  MWPower = -6; pi_time  = 44; MWFreq   = 3007.65e6   #NV D1 ms-1
            # MWI_channel = 0; MWQ_channel = 0; MWswitch_channel = 15; laserRead_channel = 5
            
            # SRSnum2 = 1; MWPower2 = -2.7;  pi_time2 = 44; MWFreq2  = 2747.88e6   #NV D1 ms+1
            # MWI2_channel = 1; MWQ2_channel = 0; MWswitch2_channel = 2
        else:
            velNum = 2; vel_current = 67; vel_wvl = 636.83; vel_vpz_target = 76.56
            ifMWDuringRead = 1
            SRSnum  = 2; MWPower = -1; pi_time  = 44; MWFreq   = 2838.26e6   #NV D2, ms-1
            MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11; laserRead_channel = 14
            ifMW2DuringRead = 1
            SRSnum2 = 4; MWPower2 = 10;   pi_time2 = 40; MWFreq2  = 2932.8e6   #NV D2 ms+1 #power=10
            MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16

        tausArray = np.array((15e3,12e3, 10.5e3, 9000, 7500,6000,4500,3000)) # sweep ti
        # tausArray = np.array((100,300,400,500,600,700,800,1000,1200,1500)) # sweep MWmix_duration_short 
        # tausArray = np.array((0.8e3, 1e3, 1.5e3,2e3,3e3,4e3))#,5e3,6e3,7e3,8e3,10e3)) # sweep delay_between_MWmix 
        # tausArray = np.array((1,5,10,30,50,75,100)) # sweep nSpinInit

        num_loops                    = int(1.5e3)
        laser_init_delay             = 1e3;       laser_init_duration       = int(1.3e6)
        laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
        RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
        laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 500
        pi_to_ion_delay              = 100;       ion_duration              = -1
        ion_to_read_delay            = 3e6;       DAQ_duration              = 2e6
        DAQ_to_laser_off_delay       = 1e2;       laserRead_to_MWmix        = RRLaserSwitch_delay
        iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
        MWmix_duration_short         = -1;        delay_between_MWmix       = -1
        spinInit_RR_duration         = 20e3;      spinInit_RR_to_pi_delay   = RRLaserSwitch_delay+100
        spinInit_pi_to_RR_delay      = 100

        settings = {'tausArray': tausArray,
                    'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
                    'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                    'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
                    'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
                    'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                    'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                    'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                    'laser_to_pi_delay':      laser_to_pi_delay,       'DAQ_duration': DAQ_duration,
                    'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'pi_time': pi_time, 'pi_time2':pi_time2,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'ion_to_read_delay': ion_to_read_delay,
                    'laserRead_to_MWmix':     laserRead_to_MWmix,      'iznLaserSwitch_delay':iznLaserSwitch_delay,
                    'ifRandomized':           ifRandomized,            'pi_to_ion_delay': pi_to_ion_delay,
                    'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                    'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,
                    'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    
                    'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                    'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                    'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                    'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
                    'hiLoMWPwr_channel':      hiLoMWPwr_channel}

        start = time.time()
        SCCRRPhotonStatIrberObject = SCCRRPhotonStatIrber(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
        SCCRRPhotonStatIrberObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = SCCRRPhotonStatIrberObject.getDataFilename()
        SCCRRPhotonStatIrberObject.close()
    




