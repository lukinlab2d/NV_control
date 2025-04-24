"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T1SCCRRIrberDualNVSameFreq import *
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A
###########################################################################################e#########################
reps = 1e5
for iii in range(int(reps)):
    # T1SCCRRIrberDualNVSameFreq
    ifRandomized = 1;  ifPlotPulse=(reps==-1);  ifMWReadLowDutyCycle = 0; ifNeedVel = 0
    ifFancySpinInit = 0; sweepWhich = 'pi_to_ion_delay'; ifIznNV2later=0
    ifMWDuringRead=1; ifMW2DuringRead=1; ifAWG=1; ifIQ=ifAWG; ifHiloExtra=1
    laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0
    ifRndPhaseNoise = 0; AGBW = 100e3; AGfreq = 1.5625e6; AGamp = 0.25 # beware of heating!!
    
    tausArray = np.round(np.logspace(np.log10(3e3),np.log10(30e6),14),-2) # sweep pi_to_ion_delay
    tausArray[-1] = 25e6
    ##############################################################################################################
    if True:
        # NV1
        velNum = 1; vel_current = 62.7; vel_wvl = 637.22; vel_vpz_target = -1; laserRead_channel = 5
        SRSnum  = 1; MWPower  = -6.2; pi_time  = 40; MWFreq  = 2953.76e6; MWPnoOvlap  = -6.2  #NV D1 ms-1
        SRSnum3 = 3; MWPower3 = 6.00; pi_time3 = 40; MWFreq3 = 2886.745e6                     #NV D1 ms+1
        MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch3_channel = 15
        SDGnum = 1; AWG_channel = 18; srate = 250000000; amp_MW_mix = 0.7
        # NV2
        velNum2 = 2; vel_current2 = 67.1; vel_wvl2 = 636.88; vel_vpz_target2=-1; laserRead2_channel = 14
        
    ##############################################################################################################
    num_loops                    = int(1e4);  nSpinInit                 = 0
    laser_init_delay             = 1e2;       laser_init_duration       = int(1.5e6)
    laserSwitch_delay_directory  = {3:860, 6:1160, 9:1160, 7:900, 5:1660, 10:170, 14:800}
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)   
    RRLaser2Switch_delay         = laserSwitch_delay_directory.get(laserRead2_channel, 0)        
    laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 5e2
    pi_to_ion_delay              = 2e2;       pi_to_hilo_extra_delay    = 1.5e3
    ion_duration                 = 4e3;       ion_duration2             = 1.5e3
    ion_to_read_delay            = 3e3;       DAQ_duration              = 3e6
    DAQ_to_laser_off_delay       = 1e2;       shift_btwn_2NV_read       = DAQ_duration+3e3
    AWG_buffer                   = 40;        AWG_output_delay          = 1450
    iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    laserRead_to_MWmix           = 0;         shift_btwn_2NV_MW         = 80
    spinInit_RR_duration         = 0;         spinInit_RR_to_pi_delay   = 0*(RRLaserSwitch_delay+200)
    spinInit_pi_to_RR_delay      = 0;         pi_pulse_bf_read          = None
    MWmix_duration_short         = 0;         delay_between_MWmix       = 0

    settings = {'tausArray': tausArray,
                'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
                'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                 'ifMW2DuringRead':ifMW2DuringRead,
                'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_pi_delay':      laser_to_pi_delay,       'DAQ_duration': DAQ_duration,
                'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'pi_time': pi_time,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'ion_to_read_delay': ion_to_read_delay,
                'laserRead_to_MWmix':     laserRead_to_MWmix,      'iznLaserSwitch_delay':iznLaserSwitch_delay,
                'ifRandomized':           ifRandomized,            'pi_to_ion_delay': pi_to_ion_delay,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,'ion_duration2': ion_duration2,
                'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    'shift_btwn_2NV_MW':shift_btwn_2NV_MW,
                'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
                'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'shift_btwn_2NV_read':shift_btwn_2NV_read,
                'laserRead2_channel':     laserRead2_channel,      'pi_time3':pi_time3, 'ifIznNV2later':ifIznNV2later,
                'MWPower3':MWPower3, 'MWFreq3': MWFreq3, 'SRSnum3': SRSnum3,
                'MWswitch3_channel': MWswitch3_channel,
                'RRLaser2Switch_delay':RRLaser2Switch_delay, 'pi_pulse_bf_read':pi_pulse_bf_read,
                'ifRndPhaseNoise':ifRndPhaseNoise, 'AGBW':AGBW, 'AGfreq':AGfreq, 'AGamp':AGamp,
                'ifAWG':ifAWG, 'ifIQ':ifIQ, 'pi_to_hilo_extra_delay':pi_to_hilo_extra_delay,'ifHiloExtra':ifHiloExtra,
                'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
                'srate':srate, 'amp_MW_mix':amp_MW_mix,
                }

    ####### Random-phase noise ######
    if ifRndPhaseNoise==1:
        AG = AG33522A()
        AG.disable_PM()
        AG.disable_RFOutput()
        AG.disable_PM(channel=2)
        AG.disable_RFOutput(channel=2)
        
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
        AG.disable_PM(channel=2)
        AG.disable_RFOutput(channel=2)
    #################################
    start = time.time()
    T1SCCRRIrberDualNVSameFreqObject = T1SCCRRIrberDualNVSameFreq(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    T1SCCRRIrberDualNVSameFreqObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = T1SCCRRIrberDualNVSameFreqObject.getDataFilename()
    T1SCCRRIrberDualNVSameFreqObject.close()

    if np.mod(iii,3)==2:
        time.sleep(60)
    
    if ifRndPhaseNoise==1:
        AG = AG33522A()
        AG.disable_PM()
        AG.disable_RFOutput()
        AG.disable_PM(channel=2)
        AG.disable_RFOutput(channel=2)
