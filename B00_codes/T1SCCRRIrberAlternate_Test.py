"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T1SCCRRIrber import *
import B00_codes.dataReader as dataReader

###########################################################################################e#########################
reps = -1
for iii in range(int(1e5)):
    # T1SCCRRIrber
    ifRandomized = 0;  ifPlotPulse = (reps==1);  ifMWReadLowDutyCycle = 0; ifNeedVel = 0
    ifFancySpinInit = 0; nSpinInit = -1; sweepWhich = 'pi_to_ion_delay'
    laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0
    
    tausArray = np.round(np.logspace(3,np.log10(100e6),21),-2) # sweep pi_to_ion_delay

    # #############################################################################################
    # if True:
    #     # NV1
    #     elNum = 1; vel_current = 62.2; vel_wvl = 637.22; vel_vpz_target = 72.3
    #     ifMWDuringRead = 1
    #     SRSnum = 1;  MWPower = -2.7; pi_time  = 44; MWFreq   = 2747.88e6   #NV D1 ms-1
    #     MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2; laserRead_channel = 5
    #     ifMW2DuringRead = 1
    #     SRSnum2 = 3; MWPower2 = -6;  pi_time2 = 44; MWFreq2  = 3007.65e6   #NV D1 ms+1
    #     MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 15

    #     num_loops                    = int(2e3)
    #     laser_init_delay             = 1e3;       laser_init_duration       = int(1e5)
    #     laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
    #     RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
    #     laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 200
    #     pi_to_ion_delay              = -1;        ion_duration              = 4800
    #     ion_to_read_delay            = 1.5e6;     DAQ_duration              = 2e6
    #     DAQ_to_laser_off_delay       = 1e2;       laserRead_to_MWmix        = RRLaserSwitch_delay
    #     iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    #     MWmix_duration_short         = 300;       delay_between_MWmix       = 2e3
    #     spinInit_RR_duration         = 15e3;      spinInit_RR_to_pi_delay   = RRLaserSwitch_delay+200
    #     spinInit_pi_to_RR_delay      = 50

    #     settings = {'tausArray': tausArray,
    #                 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
    #                 'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
    #                 'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
    #                 'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
    #                 'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
    #                 'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
    #                 'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
    #                 'laser_to_pi_delay':      laser_to_pi_delay,       'DAQ_duration': DAQ_duration,
    #                 'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'pi_time': pi_time, 'pi_time2':pi_time2,
    #                 'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'ion_to_read_delay': ion_to_read_delay,
    #                 'laserRead_to_MWmix':     laserRead_to_MWmix,      'iznLaserSwitch_delay':iznLaserSwitch_delay,
    #                 'ifRandomized':           ifRandomized,            'pi_to_ion_delay': pi_to_ion_delay,
    #                 'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
    #                 'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,
    #                 'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,  
    #                 'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
    #                 'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
    #                 'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
    #                 'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
    #                 'hiLoMWPwr_channel':      hiLoMWPwr_channel}

    #     start = time.time()
    #     T1SCCRRIrberObject = T1SCCRRIrber(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    #     T1SCCRRIrberObject.runScan()
    #     print('Total time = ' + str(time.time() - start) + ' s')
        
    #     dataFilename = T1SCCRRIrberObject.getDataFilename()
    #     T1SCCRRIrberObject.close()

    #############################################################################################
    # NV2
    velNum = 2; vel_current = 67; vel_wvl = 636.83; vel_vpz_target = 76.56
    ifMWDuringRead = 1
    SRSnum  = 2; MWPower = -1; pi_time  = 44; MWFreq   = 2838.26e6   #NV D2, ms-1
    MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11; laserRead_channel = 14
    ifMW2DuringRead = 1
    SRSnum2 = 4; MWPower2 = 10;   pi_time2 = 44; MWFreq2  = 2932.8e6   #NV D2 ms+1
    MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16
    
    num_loops                    = int(2e3)
    laser_init_delay             = 1e3;       laser_init_duration       = int(1e5)
    laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
    laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 200
    pi_to_ion_delay              = -1;        ion_duration              = 4800
    ion_to_read_delay            = 1.5e6;     DAQ_duration              = 2e6
    DAQ_to_laser_off_delay       = 1e2;       laserRead_to_MWmix        = RRLaserSwitch_delay
    iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    MWmix_duration_short         = 300;       delay_between_MWmix       = 2e3
    spinInit_RR_duration         = 15e3;      spinInit_RR_to_pi_delay   = RRLaserSwitch_delay+200
    spinInit_pi_to_RR_delay      = 50

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
    T1SCCRRIrberObject = T1SCCRRIrber(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    T1SCCRRIrberObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = T1SCCRRIrberObject.getDataFilename()
    T1SCCRRIrberObject.close()




