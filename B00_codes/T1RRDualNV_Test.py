"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T1RRDualNV import *
import B00_codes.dataReader as dataReader
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A

###########################################################################################e#########################
reps = -1
for iii in range(int(1e5)):
    # T1RRDualNV
    ifRandomized = 1;  ifPlotPulse = (reps==1);  ifNeedVel = 0
    ifFancySpinInit = 0; sweepWhich = 'pi_to_ion_delay'
    laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0
    ifRndPhaseNoise = 0; AGBW = 1500e3; AGfreq = 3e6; AGamp = 0.5 # beware of heating!!
    
    tausArray = np.round(np.logspace(3,np.log10(100e6),21),-2) # sweep pi_to_ion_delay

    ##############################################################################################################
    if True:
        # NV1
        velNum = 1; vel_current = 62.2; vel_wvl = 637.22; vel_vpz_target = 72.3; laserRead_channel = 5

        SRSnum  = 1; MWPower  = -2.7; pi_time  = 44; MWFreq   = 2747.88e6   #NV D1 ms-1
        SRSnum3 = 3; MWPower3 = -6;   pi_time3 = 44; MWFreq3  = 3007.65e6   #NV D1 ms+1
        MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch3_channel = 15
        
        # NV2
        velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.83; vel_vpz_target2 = 76.56; laserRead2_channel = 14
        
        SRSnum2 = 2; MWPower2 = -1;   pi_time2 = 44; MWFreq2  = 2838.26e6   #NV D2, ms-1
        SRSnum4 = 4; MWPower4 = 10;   pi_time4 = 44; MWFreq4  = 2932.8e6    #NV D2 ms+1
        MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11; MWswitch4_channel = 16
    ##############################################################################################################

    num_loops                    = int(2e3);  nSpinInit                 = -1
    laser_init_delay             = 1e2;       laser_init_duration       = int(30e3)
    laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)   
    RRLaser2Switch_delay         = laserSwitch_delay_directory.get(laserRead2_channel, 0)        
    laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 200
    pi_to_ion_delay              = -1;        DAQ_duration              = 10e3
    DAQ_to_laser_off_delay       = 1e2;       shift_btwn_2NV_read       = DAQ_duration+1e3
    spinInit_RR_duration         = 15e3;      spinInit_RR_to_pi_delay   = RRLaserSwitch_delay+200
    spinInit_pi_to_RR_delay      = 50;        

    settings = {'tausArray': tausArray,
                'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 
                'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
                'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_pi_delay':      laser_to_pi_delay,       'DAQ_duration': DAQ_duration,
                'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'pi_time': pi_time, 'pi_time2':pi_time2,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  
                'ifRandomized':           ifRandomized,            'pi_to_ion_delay': pi_to_ion_delay,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'laserIon_channel':       laserIon_channel,        
                'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
                'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'shift_btwn_2NV_read':shift_btwn_2NV_read,
                'laserRead2_channel':     laserRead2_channel,
                'MWPower3':MWPower3, 'MWFreq3': MWFreq3, 'SRSnum3': SRSnum3,
                'MWPower4':MWPower4, 'MWFreq4': MWFreq4, 'SRSnum4': SRSnum4,
                'MWswitch3_channel': MWswitch3_channel,'MWswitch4_channel': MWswitch4_channel,
                'RRLaser2Switch_delay':RRLaser2Switch_delay}

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
    T1RRDualNVObject = T1RRDualNV(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    T1RRDualNVObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = T1RRDualNVObject.getDataFilename()
    T1RRDualNVObject.close()

    if ifRndPhaseNoise:
        AG.disable_PM()
        AG.disable_RFOutput()
