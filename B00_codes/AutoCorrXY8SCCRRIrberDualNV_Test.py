"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.AutoCorrXY8SCCRRIrberDualNV import *
import dropbox
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2; pi = np.pi

###########################################################################################e#########################
reps = 1e5
for iii in range(int(reps)):
    # AutoCorrXY8SCCRRIrberDualNV
    ifRandomized=0; ifPlotPulse=(reps==-1); ifMWReadLowDutyCycle=0; ifNeedVel=0; ifStartInY=0
    ifFancySpinInit=0; normalized_style = Q_FINAL; ifAntiCorrel=0; ifSinDetect=1; ifAWG=1
    ifMWDuringRead=1; ifMW2DuringRead=1; ifJustRef_CorrACorr=1; if_tracking=0 # ifSinDetect, specify phase of last pulse
    laserInit_channel=3; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0
    ifTestSig=1; ifRndPhaseNoise=1; AGBW = 100e3; AGfreq = 1.5625e6; AGamp = 0.18 # beware of heating!!
    sweepWhich='shift_btwn_2NV_MW'; NXY8=3

    offset = 0; tausArray = np.linspace(offset,offset+1920,25)
    ##############################################################################################################
    if True:
        # NV1
        velNum = 1; vel_current = 62.7; vel_wvl = 637.20; vel_vpz_target = -1; laserRead_channel = 5
        SRSnum  = 1; MWPower  = -6.0; pi_time  = 40;  MWFreq   = 2598.1e6 #NV D1 ms-1
        SRSnum3 = 3; MWPower3 = 3;    pi_time3 = -1;  MWFreq3  = 3162e6   #NV D1 ms+1
        MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch3_channel = 15
        SDGnum = 1; AWG_channel = 18; srate = 2.5e8; amp_MW_mix = 1#0.83
        # NV2
        velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.88; vel_vpz_target2 = -1; laserRead2_channel = 14
        SRSnum2 = 2; MWPower2 = -9.8;  pi_time2 = 40;  MWFreq2  = 2789.2e6  #NV D2, ms-1
        SRSnum4 = 4; MWPower4 = -4;    pi_time4 = -1;  MWFreq4  = 3037.2e6  #NV D2 ms+1
        MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11; MWswitch4_channel = 16
        SDGnum2 = 2; AWG2_channel = 19; srate2 = 2.5e8; amp_MW_mix2 = 1#0.83
    ##############################################################################################################
    num_loops                    = int(10e3); nSpinInit                 = 0
    laser_init_delay             = 1e2;       laser_init_duration       = int(2.5e6)
    laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1650, 10:170, 14:800}
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)   
    RRLaser2Switch_delay         = laserSwitch_delay_directory.get(laserRead2_channel, 0)        
    laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 5e2
    pi_to_ion_delay              = 5e2;       MWI_to_switch_delay       = 30
    ion_duration                 = 6e3;       ion_duration2             = 6e3
    ion_to_read_delay            = 3e3;       DAQ_duration              = 2.5e6
    DAQ_to_laser_off_delay       = 1e2;       shift_btwn_2NV_read       = DAQ_duration+3e3
    AWG_buffer                   = 40;        AWG_output_delay          = 1450       
    iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    laserRead_to_MWmix           = 0;         phi_IQ = 0.05*pi; phi_IQ2 = 0.05*pi 
    tauExtra                     = 0;         phi_IQ_antiCorrPulse      = 0.5*pi
    spinInit_RR_duration         = 0;         spinInit_RR_to_pi_delay   = 40
    spinInit_pi_to_RR_delay      = 0;         tau                       = int(1/(2*AGfreq)*1e9-40)
    shift_btwn_2NV_MW            = 80;        shift_btwn_2NV_read       = DAQ_duration+3e3
    MWmix_duration_short         = -1;        delay_between_MWmix       = -1
    tcorr                        = 640

    if True:
        ######################################### Scan params ########################################
        threshold_scanVpz            = 2.4;       threshold_scanVpz2        = 3.3
        num_loops_track              = int(2e4);   
        laser_init_delay_track       = 1e2;       laser_init_duration_track = 15e3
        MW_to_read_delay_track       = 1e2
        laser_to_DAQ_delay_directory = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
        laser_to_MWI_delay_track     = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
        read_duration_track          = 12300;      read_laser_duration_track  = 12200
        ######################################################
        time_sleep_after_scan        = 35;         wvl_correction             = 8e-6
        scan_lower_margin = 0.15; scan_upper_margin = 0.1; num_point_scan = 26

    # Create a Dropbox client
    tkFile = 'C:/Users/lukin2dmaterials/data/tk.txt'; lines = []
    with open(tkFile, 'r') as file:
        for line in file:
            if "\n" in line: line = line[0:-1]
            lines.append(line)
    access_token = lines[0]
    app_key = lines[1]
    app_secret = lines[2]
    refresh_token = lines[3]
    #######################################
    dbx = dropbox.Dropbox(oauth2_access_token=access_token,
                          app_key=app_key, 
                          app_secret=app_secret, 
                          oauth2_refresh_token=refresh_token)
    RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_scanVpz': threshold_scanVpz,
                    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                    'SRSnum':   SRSnum,
                    'laser_init_delay':       laser_init_delay_track,     'laser_init_duration': laser_init_duration_track,
                    'laser_to_MWI_delay':     laser_to_MWI_delay_track,   'MWI_duration':        pi_time,
                    'laser_to_DAQ_delay':     RRLaserSwitch_delay,   'read_duration':       read_duration_track,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'read_laser_duration':    read_laser_duration_track,   
                    'MW_to_read_delay':       MW_to_read_delay_track,   
                    'vel_current':            vel_current,           'vel_wvl':             vel_wvl,
                    'velNum': velNum,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':1,
                    'ifInitWvl':ifInitWvl,    'dbx':dbx,
                    'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                    'time_sleep_after_scan': time_sleep_after_scan, 'wvl_correction': wvl_correction,
                    'scan_lower_margin': scan_lower_margin, 'scan_upper_margin': scan_upper_margin, 'num_point_scan':num_point_scan,}
    RRtrackingSettings2 = {'if_tracking': if_tracking, 'threshold_scanVpz': threshold_scanVpz2,
                    'num_loops':num_loops_track,  'MWPower':MWPower2,    'MWFreq': MWFreq2,
                    'SRSnum':   SRSnum2,
                    'laser_init_delay':       laser_init_delay_track,     'laser_init_duration': laser_init_duration_track,
                    'laser_to_MWI_delay':     laser_to_MWI_delay_track,   'MWI_duration':        pi_time2,
                    'laser_to_DAQ_delay':     RRLaser2Switch_delay,   'read_duration':       read_duration_track,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead2_channel,
                    'read_laser_duration':    read_laser_duration_track,   
                    'MW_to_read_delay':       MW_to_read_delay_track,   
                    'vel_current':            vel_current2,           'vel_wvl':             vel_wvl2,
                    'velNum': velNum2,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':1,
                    'ifInitWvl':ifInitWvl,    'dbx':dbx,
                    'MWI_channel': MWI2_channel,  'MWQ_channel': MWQ2_channel,  'MWswitch_channel': MWswitch2_channel,
                    'time_sleep_after_scan': time_sleep_after_scan, 'wvl_correction': wvl_correction,
                    'scan_lower_margin': scan_lower_margin, 'scan_upper_margin': scan_upper_margin, 'num_point_scan':num_point_scan,}
    settings = {'tausArray': tausArray,
                'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
                'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
                'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
                'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                'velNum1': velNum2, 'vel_current2': vel_current2, 'vel_wvl2': vel_wvl2, 'vel_vpz_target2': vel_vpz_target2,
                'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_pi_delay':      laser_to_pi_delay,       'DAQ_duration': DAQ_duration,
                'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'pi_time': pi_time, 'pi_time3':pi_time3,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'ion_to_read_delay': ion_to_read_delay,
                'laserRead_to_MWmix':     laserRead_to_MWmix,      'iznLaserSwitch_delay':iznLaserSwitch_delay,
                'ifRandomized':           ifRandomized,            'pi_to_ion_delay': pi_to_ion_delay,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,'ion_duration2': ion_duration2,
                'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    'MWI_to_switch_delay': MWI_to_switch_delay,
                'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'normalized_style':normalized_style,
                'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'shift_btwn_2NV_read':shift_btwn_2NV_read,
                'laserRead2_channel':     laserRead2_channel,      'NXY8':NXY8,
                'MWPower3':MWPower3, 'MWFreq3': MWFreq3, 'SRSnum3': SRSnum3,
                'MWPower4':MWPower4, 'MWFreq4': MWFreq4, 'SRSnum4': SRSnum4,
                'MWswitch3_channel': MWswitch3_channel,'MWswitch4_channel': MWswitch4_channel,
                'RRLaser2Switch_delay':RRLaser2Switch_delay,
                'RRtrackingSettings':RRtrackingSettings, 'RRtrackingSettings2':RRtrackingSettings2,
                'ifAntiCorrel':ifAntiCorrel, 'ifSinDetect':ifSinDetect, 'ifJustRef_CorrACorr':ifJustRef_CorrACorr,
                'shift_btwn_2NV_MW':shift_btwn_2NV_MW,'ifStartInY':ifStartInY,
                'ifAWG':ifAWG,'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
                'SDGnum2': SDGnum2, 'AWG2_channel':AWG2_channel, 'srate':srate, 'srate2':srate2,
                'phi_IQ':phi_IQ, 'phi_IQ2':phi_IQ2, 'phi_IQ_antiCorrPulse':phi_IQ_antiCorrPulse,
                'tauExtra':tauExtra, 'tcorr':tcorr,
                'ifRndPhaseNoise':ifRndPhaseNoise, 'AGBW':AGBW, 'AGfreq':AGfreq, 'AGamp':AGamp,'ifTestSig':ifTestSig,
                'tau':tau, 'sweepWhich':sweepWhich,'amp_MW_mix':amp_MW_mix,'amp_MW_mix2':amp_MW_mix2}
    ####### Random-phase noise ######
    if ifTestSig==1:
        AG = AG33522A()
        AG.disable_PM()
        AG.disable_RFOutput()
        
        if ifRndPhaseNoise==1:
            AG.set_PMsource()
            AG.set_PMfunction(function='NOIS')
            AG.set_PMdeviation()
            AG.set_noiseBandwidth(bandwidth=AGBW)

        AG.apply(function='SIN', freq=AGfreq, amplitude=AGamp, DCoffset=0)
        if ifRndPhaseNoise==1: AG.enable_PM()
    else:
        AG = AG33522A()
        AG.disable_PM()
        AG.disable_RFOutput()
    #################################
    start = time.time()
    AutoCorrXY8SCCRRIrberDualNVObject = AutoCorrXY8SCCRRIrberDualNV(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    AutoCorrXY8SCCRRIrberDualNVObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    AutoCorrXY8SCCRRIrberDualNVObject.close()

    # if ifTestSig==1:
    #     AG.disable_PM()
    #     AG.disable_RFOutput()