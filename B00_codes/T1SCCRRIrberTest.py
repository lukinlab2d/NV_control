"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T1SCCRRIrber import *

###########################################################################################e#########################
reps = 1
for iii in range(int(reps)):
#     # T1SCCRRIrber
#     ifRandomized=0; ifPlotPulse=(reps==-1); ifMWReadLowDutyCycle=0; ifNeedVel=0
#     ifFancySpinInit=0; nSpinInit=0; sweepWhich = 'pi_to_ion_delay'
#     laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0
#     ifMWDuringRead=1; ifMW2DuringRead=1; ifAWG=1; ifIQ=ifAWG; ifHiloExtra=1
#     if True: 
#         velNum = 1; vel_current = 62.7; vel_wvl = 637.20; vel_vpz_target = -1; laserRead_channel = 5
#         SRSnum  = 1; MWPower  = -7.1; pi_time  = 40;  MWFreq   = 2598.1e6 #NV D1 ms-1
#         SRSnum2 = 3; MWPower2 = 3;    pi_time2 = -1;  MWFreq2  = 3162e6   #NV D1 ms+1
#         MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch2_channel = 15
#         SDGnum = 1; AWG_channel = 18; srate = 2.5e8
#     else:
#         velNum = 2; vel_current = 67; vel_wvl = 636.88; vel_vpz_target = -1; laserRead_channel = 14
#         SRSnum  = 2; MWPower  = -9.5; pi_time  = 126; MWFreq   = 2789.2e6  #NV D2, ms-1
#         SRSnum2 = 4; MWPower2 = -2;   pi_time2 = -1;  MWFreq2  = 3037.2e6  #NV D2 ms+1
#         MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11; MWswitch2_channel = 16
#         SDGnum = 2; AWG_channel = 19
    
#     tausArray = np.round(np.logspace(np.log10(3e3),np.log10(27e6),21),-2) # sweep pi_to_ion_delay

#     num_loops                    = int(2e3)
#     laser_init_delay             = 1e2;       laser_init_duration       = int(1050e3)
#     laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:170, 14:900}   
#     laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 5e2
#     pi_to_ion_delay              = -1;        ion_duration              = 4e3 
#     ion_to_read_delay            = 1e5;       DAQ_duration              = 2.8e6 
#     DAQ_to_laser_off_delay       = 1e2;       pi_to_hilo_extra_delay    = 3.5e3
#     AWG_buffer                   = 80;        AWG_output_delay          = 1450
#     RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)    
#     iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
#     laserRead_to_MWmix           = 0
#     MWmix_duration_short         = 0;         delay_between_MWmix       = 0
#     spinInit_RR_duration         = 0;         spinInit_RR_to_pi_delay   = 0*(RRLaserSwitch_delay+200)
#     spinInit_pi_to_RR_delay      = 0

#     settings = {'tausArray': tausArray,
#                 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
#                 'MWswitch_channel': MWswitch_channel, 'hiLoMWPwr_channel':hiLoMWPwr_channel,
#                 'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
#                 'MWswitch2_channel': MWswitch2_channel,
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
#                 'ifAWG':ifAWG, 'ifIQ':ifIQ,'pi_to_hilo_extra_delay':pi_to_hilo_extra_delay,'ifHiloExtra':ifHiloExtra,
#                 'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
#                 # 'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 
#                 'srate':srate}

#     start = time.time()
#     T1SCCRRIrberObject = T1SCCRRIrber(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
#     T1SCCRRIrberObject.runScan()
#     print('Total time = ' + str(time.time() - start) + ' s')
#     T1SCCRRIrberObject.close()




    # T1SCCRRIrber
    ifRandomized=0; ifPlotPulse=(reps==-1); ifMWReadLowDutyCycle=0; ifNeedVel=0
    ifFancySpinInit=0; nSpinInit=0; sweepWhich = 'pi_to_ion_delay'
    laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0
    ifMWDuringRead=1; ifMW2DuringRead=1; ifAWG=1; ifIQ=ifAWG; ifHiloExtra=1
    if False: 
        velNum = 1; vel_current = 62.7; vel_wvl = 637.20; vel_vpz_target = -1; laserRead_channel = 5
        SRSnum  = 1; MWPower  = -6; pi_time  = 138; MWFreq   = 2598.1e6 #NV D1 ms-1
        SRSnum2 = 3; MWPower2 = 2;  pi_time2 = -1;  MWFreq2  = 3162e6   #NV D1 ms+1
        MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch2_channel = 15
        SDGnum = 1; AWG_channel = 18
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.88; vel_vpz_target = -1; laserRead_channel = 14
        SRSnum  = 2; MWPower  = -11.4; pi_time = 40;  MWFreq   = 2789.2e6  #NV D2, ms-1
        SRSnum2 = 4; MWPower2 = -4;    pi_time2 = -1; MWFreq2  = 3037.2e6  #NV D2 ms+1
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11; MWswitch2_channel = 16
        SDGnum = 2; AWG_channel = 19; srate = 2.5e8
    
    tausArray = np.round(np.logspace(np.log10(3e3),np.log10(27e6),21),-2) # sweep pi_to_ion_delay

    num_loops                    = int(2e3)
    laser_init_delay             = 1e6;       laser_init_duration       = int(1050e3)
    laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:170, 14:900}
    laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 5e2
    pi_to_ion_delay              = -1;        ion_duration              = 13e3 
    ion_to_read_delay            = 3e3;       DAQ_duration              = 26e5
    DAQ_to_laser_off_delay       = 1e2;       pi_to_hilo_extra_delay    = 3.5e3
    AWG_buffer                   = 80;        AWG_output_delay          = 1450
    iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0) 
    laserRead_to_MWmix           = 0
    MWmix_duration_short         = 0;         delay_between_MWmix       = 0
    spinInit_RR_duration         = 0;         spinInit_RR_to_pi_delay   = 0*(RRLaserSwitch_delay+200)
    spinInit_pi_to_RR_delay      = 0

    settings = {'tausArray': tausArray,
                'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
                'MWswitch_channel': MWswitch_channel, 'hiLoMWPwr_channel':hiLoMWPwr_channel,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
                'MWswitch2_channel': MWswitch2_channel,
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
                'ifAWG':ifAWG, 'ifIQ':ifIQ, 'pi_to_hilo_extra_delay':pi_to_hilo_extra_delay,'ifHiloExtra':ifHiloExtra,
                'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
                # 'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 
                'srate':srate}

    start = time.time()
    T1SCCRRIrberObject = T1SCCRRIrber(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    T1SCCRRIrberObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    T1SCCRRIrberObject.close()




