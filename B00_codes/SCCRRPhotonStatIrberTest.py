"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.SCCRRPhotonStatIrber import *
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A

####################################################################################################################
reps = 1; ifPlotPulse=(reps==-1)
for whichNV in np.linspace(2,2,1):
    # SCCRRPhotonStatIrber
    ifNeedVel=0; ifInitVpz=0; ifInitWvl=0; ifFancySpinInit=0; 
    ifAWG=1; ifMWDuringRead=1; ifMW2DuringRead=1; ifMWReadLowDutyCycle=0; ifIQ=ifAWG
    ifRandomized = 0;  sweepWhich = 'rabi'; 
    laserInit_channel = 3; laserIon_channel = 10; hiLoMWPwr_channel = 17; 
    ifRndPhaseNoise = 0; AGBW = 100e3; AGfreq = 1506024; AGamp = 0.25 # beware of heating!!
    if np.mod(whichNV,2)==1: 
        velNum = 1; vel_current = 62.7; vel_wvl = 637.20; vel_vpz_target=-1; laserRead_channel = 5
        SRSnum = 1; MWPower = -5.5; pi_time  = 36; MWFreq = 2598.1e6   #NV D1 ms-1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
        SDGnum = 1; AWG_channel = 18; srate = 2.5e8
        SRSnum2= 3; MWPower2 = 3; pi_time2 = -1; MWFreq2= 3162e6   #NV D1 ms+1
        MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 15
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.88; vel_vpz_target=-1; laserRead_channel = 14
        SRSnum = 2; MWPower = -9.4; pi_time  = 36; MWFreq = 2789.2e6   #NV D2, ms-1
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11
        SDGnum = 2; AWG_channel = 19; srate = 2.5e8
        SRSnum2= 4; MWPower2 = -4;  pi_time2 = -1; MWFreq2= 3037.2e6   #NV D2 ms+1 #power=10
        MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16
    
    # tausArray = np.array((1e3,1.5e3,2e3,3e3,4e3,5e3,6e3,8e3,10e3,12e3,15e3,22e3,30e3)) # sweep ti
    # tausArray = np.array((2e6,2.5e6,3e6,4e6,5e6,6e6)) # sweep tr
    tausArray = np.linspace(32,44,4)         # rabi
    # tausArray = np.array((1e6,1.5e6,2e6,2.5e6,3e6,5e6,7e6,10e6)) # tinit
    # tausArray = np.array((5e2,1e3,3e3, 5e3, 1e4, 3e4, 5e4,1e5,3e5,5e5,7e5,9e5,11e5,14e5,17e5,20e5)) #i2r
    num_loops                    = int(1.5e3); nSpinInit                 = 0
    laser_init_delay             = 1e2;        laser_init_duration       = int(2.5e6)
    laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1650, 10:170, 14:800}
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
    laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 5e2
    pi_to_ion_delay              = 5e2;        ion_duration              = 5e3
    ion_to_read_delay            = 3e2;        DAQ_duration              = 3e6
    DAQ_to_laser_off_delay       = 1e2;        laserRead_to_MWmix        = 0
    iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    AWG_buffer                   = 40;         AWG_output_delay          = 1450

    MWmix_duration_short         = 0;          delay_between_MWmix       = 0
    spinInit_RR_duration         = 0;          spinInit_RR_to_pi_delay   = 0*(RRLaserSwitch_delay+100)
    spinInit_pi_to_RR_delay      = 0
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
                'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    'srate': int(srate),
                'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
                'hiLoMWPwr_channel':      hiLoMWPwr_channel,
                'ifAWG': ifAWG, 'ifIQ': ifIQ,
                'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,}
    ####### Random-phase noise ######
    if ifRndPhaseNoise==1:
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
    SCCRRPhotonStatIrberObject = SCCRRPhotonStatIrber(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    SCCRRPhotonStatIrberObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    SCCRRPhotonStatIrberObject.close()

    # if ifRndPhaseNoise==1:
    #     AG.disable_PM()
    #     AG.disable_RFOutput()



