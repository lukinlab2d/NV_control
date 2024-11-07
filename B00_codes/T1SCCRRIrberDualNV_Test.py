"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T1SCCRRIrberDualNV import *
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A
###########################################################################################e#########################
reps = 5e5
for iii in range(int(reps)):
    # T1SCCRRIrberDualNV
    ifRandomized = 0;  ifPlotPulse=(reps==-1);  ifMWReadLowDutyCycle = 0; ifNeedVel = 0
    ifFancySpinInit = 0; sweepWhich = 'pi_pulse'#'pi_to_ion_delay'
    ifMWDuringRead=1; ifMW2DuringRead=1; ifAWG=1; ifIQ=ifAWG; ifHiloExtra=1
    laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0
    ifRndPhaseNoise = 1; AGBW = 100e3; AGfreq = 1.5625e6; AGamp = 0.25 # beware of heating!!
    
    tausArray = np.round(np.logspace(np.log10(3e3),np.log10(15e6),21),-2) # sweep pi_to_ion_delay
    ##############################################################################################################
    if True:
        # NV1
        velNum = 1; vel_current = 62.7; vel_wvl = 637.20; vel_vpz_target = -1; laserRead_channel = 5
        SRSnum  = 1; MWPower  = -6.4; pi_time  = -1;  MWFreq   = 2598.1e6 #NV D1 ms-1
        SRSnum3 = 3; MWPower3 = 3;    pi_time3 = -1;  MWFreq3  = 3162e6   #NV D1 ms+1
        MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch3_channel = 15
        SDGnum = 1; AWG_channel = 18; srate = 2.5e8

        # NV2
        velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.88; vel_vpz_target2 = -1; laserRead2_channel = 14
        SRSnum2 = 2; MWPower2 = -10.2; pi_time2 = -1;  MWFreq2  = 2789.2e6  #NV D2, ms-1
        SRSnum4 = 4; MWPower4 = -4;    pi_time4 = -1;  MWFreq4  = 3037.2e6  #NV D2 ms+1
        MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11; MWswitch4_channel = 16
        SDGnum2 = 2; AWG2_channel = 19; srate2 = 2.5e8
    ##############################################################################################################
    num_loops                    = int(2e3);  nSpinInit                 = 0
    laser_init_delay             = 1e2;       laser_init_duration       = int(2e6)
    laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:170, 14:900}
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)   
    RRLaser2Switch_delay         = laserSwitch_delay_directory.get(laserRead2_channel, 0)        
    laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 5e2
    pi_to_ion_delay              = 5e3;        pi_to_hilo_extra_delay    = 3.5e3
    ion_duration                 = 4000;      ion_duration2             = 13e3
    ion_to_read_delay            = 3e3;       DAQ_duration              = 28e5
    DAQ_to_laser_off_delay       = 1e2;       shift_btwn_2NV_read       = DAQ_duration+3e3
    AWG_buffer                   = 40;        AWG_output_delay          = 1450
    iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    laserRead_to_MWmix           = 0
    spinInit_RR_duration         = 0;         spinInit_RR_to_pi_delay   = 0*(RRLaserSwitch_delay+200)
    spinInit_pi_to_RR_delay      = 0;        
    MWmix_duration_short         = 0;         delay_between_MWmix       = 0

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
                'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,'ion_duration2': ion_duration2,
                'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,  
                'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
                'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'shift_btwn_2NV_read':shift_btwn_2NV_read,
                'laserRead2_channel':     laserRead2_channel,
                'MWPower3':MWPower3, 'MWFreq3': MWFreq3, 'SRSnum3': SRSnum3,
                'MWPower4':MWPower4, 'MWFreq4': MWFreq4, 'SRSnum4': SRSnum4,
                'MWswitch3_channel': MWswitch3_channel,'MWswitch4_channel': MWswitch4_channel,
                'RRLaser2Switch_delay':RRLaser2Switch_delay,
                'ifRndPhaseNoise':ifRndPhaseNoise, 'AGBW':AGBW, 'AGfreq':AGfreq, 'AGamp':AGamp,
                'ifAWG':ifAWG, 'ifIQ':ifIQ, 'pi_to_hilo_extra_delay':pi_to_hilo_extra_delay,'ifHiloExtra':ifHiloExtra,
                'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
                'SDGnum2': SDGnum2, 'AWG2_channel':AWG2_channel, 'srate':srate, 'srate2':srate2
                }

    ####### Random-phase noise ######
    if ifRndPhaseNoise:
        AG = AG33522A()
        AG.disable_PM()
        AG.disable_RFOutput()
        
        AG.set_PMsource()
        AG.set_PMfunction(function='NOIS')
        AG.set_PMdeviation()
        AG.set_noiseBandwidth(bandwidth=AGBW)

        AG.apply(function='SIN', freq=AGfreq, amplitude=AGamp, DCoffset=0)
        AG.enable_PM()
    else:
        AG = AG33522A()
        AG.disable_PM()
        AG.disable_RFOutput()
    #################################
    start = time.time()
    T1SCCRRIrberDualNVObject = T1SCCRRIrberDualNV(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    T1SCCRRIrberDualNVObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = T1SCCRRIrberDualNVObject.getDataFilename()
    T1SCCRRIrberDualNVObject.close()

    # if ifRndPhaseNoise:
    #     AG.disable_PM()
    #     AG.disable_RFOutput()
