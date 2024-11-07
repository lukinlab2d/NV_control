"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T2ERRDualNV_NewTrack import *
import time, dropbox
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3
pi = np.pi
####################################################################################################################
reps=int(1);  ifLooped=(reps!=1); laserInit_channel=3; normalized_style=Q_FINAL; 
ifInitWvl=0; ifInitVpz=0; ifNeedVel1=0; ifNeedVel2=0
ifRandomized=0; if_tracking=0; ifAWG=1
ifRndPhaseNoise = 0; AGBW = 50e3; AGfreq = 3.3e6; AGamp = 0*0.2 # beware of heating!!
for i in np.linspace(1, reps, reps):
    # T2ERRDualNV_NewTrack
    tausArray = np.round(np.logspace(2,np.log10(1e6),41),-1)

    num_loops                    = int(1e5);   phi_IQ = phi_IQ2       = pi/2 #rad, angle of last pi/2
    laser_init_delay             = 1e2;        laser_init_duration    = 8e3
    MW_to_read_delay             = 1e2;        MWI_to_switch_delay    = 30
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 2e3;        read_laser_duration    = read_duration
    read_duration2               = 2e3;        read_laser_duration2   = read_duration2
    shift_btwn_2NV_MW            = 80;         shift_btwn_2NV_read    = read_duration+1.7e3
    AWG_buffer                   = 1;          AWG_output_delay       = 1450  
    AWG_buffer2                  = 1;          AWG_output_delay2      = 1450  

    if True:
        ########### NV1 ##############
        velNum = 1; vel_current = 62.7; vel_wvl = 637.20; vel_vpz_target = -1; laserRead_channel = 5
        SRSnum  = 1; MWPower  = -6.0;  pi_half  = 20; MWFreq   = 2598.1e6   #NV D1 ms-1
        SRSnum3 = 3; MWPower3 = -100;  pi_half3 = 10; MWFreq3  = 3000e6     #NV D1 ms+1
        SDGnum  = 1; AWG_channel = 18 
        MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch3_channel = 15
        laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 

        ############ NV2 #############
        velNum2 = 2; vel_current2 = 67; vel_wvl2 = 636.88; vel_vpz_target2 = -1; laserRead2_channel = 14
        SRSnum2 = 2; MWPower2 = -9.8; pi_half2  = 20; MWFreq2  = 2789.2e6   #NV D2, ms-1
        SRSnum4 = 4; MWPower4 = -100; pi_half4  = 10; MWFreq4  = 3100e6     #NV D2, ms+1
        SDGnum2 = 2; AWG2_channel = 19
        MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11; MWswitch4_channel = 16
        laser_to_DAQ_delay2 = laser_to_DAQ_delay_directory.get(laserRead2_channel, 0) 

    ######################################### Line track params ########################################
    threshold_scanVpz            = 1.2;        threshold_scanVpz2        = 0
    num_loops_track              = int(2e4);   
    laser_init_delay_track       = 1e2;        laser_init_duration_track  = 8e3
    MW_to_read_delay_track       = 1e2
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1750, 14: 900}
    laser_to_MWI_delay_track     = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration_track          = 2000;       read_laser_duration_track  = 2000
    #######################################
    time_sleep_after_scan        = 40;         wvl_correction             = 8e-6
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
    dbx = dropbox.Dropbox(oauth2_access_token=access_token,
                          app_key=app_key, 
                          app_secret=app_secret, 
                          oauth2_refresh_token=refresh_token)
    
    ##############################################################################
    RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_scanVpz': threshold_scanVpz,
                    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                    'SRSnum':   SRSnum,
                    'laser_init_delay':       laser_init_delay_track,     'laser_init_duration': laser_init_duration_track,
                    'laser_to_MWI_delay':     laser_to_MWI_delay_track,   'MWI_duration':        2*pi_half,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,   'read_duration':       read_duration_track,
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
                    'laser_to_MWI_delay':     laser_to_MWI_delay_track,   'MWI_duration':        2*pi_half2,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay2,   'read_duration':       read_duration_track,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead2_channel,
                    'read_laser_duration':    read_laser_duration_track,   
                    'MW_to_read_delay':       MW_to_read_delay_track,   
                    'vel_current':            vel_current2,           'vel_wvl':             vel_wvl2,
                    'velNum': velNum2,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':1,
                    'ifInitWvl':ifInitWvl,    'dbx':dbx,
                    'MWI_channel': MWI2_channel,  'MWQ_channel': MWQ2_channel,  'MWswitch_channel': MWswitch2_channel,
                    'time_sleep_after_scan': time_sleep_after_scan, 'wvl_correction': wvl_correction,
                    'scan_lower_margin': scan_lower_margin, 'scan_upper_margin': scan_upper_margin, 'num_point_scan':num_point_scan,}
    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay,    'MWI_to_switch_delay':       MWI_to_switch_delay,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel, 
                'laser_to_DAQ_delay2':    laser_to_DAQ_delay2,   'laserRead2_channel':        laserRead2_channel,
                'normalized_style':       normalized_style,      'pi_half': pi_half, 'pi_half2': pi_half2,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum, 'ifNeedVel1': ifNeedVel1,
                'vel_current2':  vel_current2, 'vel_wvl2': vel_wvl2, 'velNum2': velNum2, 'ifNeedVel2': ifNeedVel2,
                'vel_vpz_target': vel_vpz_target, 'vel_vpz_target2': vel_vpz_target2, 
                'ifInitVpz':ifInitVpz,    'ifInitWvl': ifInitWvl, 'ifRandomized': ifRandomized,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'MWI2_channel': MWI2_channel,  'MWQ2_channel': MWQ2_channel,  'MWswitch2_channel': MWswitch2_channel,
                'RRtrackingSettings': RRtrackingSettings, 'RRtrackingSettings2': RRtrackingSettings2,
                'shift_btwn_2NV_MW':shift_btwn_2NV_MW, 'shift_btwn_2NV_read': shift_btwn_2NV_read, 'ifAWG':ifAWG,
                'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
                'SDGnum2': SDGnum2, 'AWG2_channel':AWG2_channel, 'AWG_buffer2':AWG_buffer2, 'AWG_output_delay2':AWG_output_delay2,
                'read_duration2':read_duration2, 'read_laser_duration2':read_laser_duration2,'phi_IQ':phi_IQ,'phi_IQ2':phi_IQ2}

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
    #################################
    start = time.time()
    T2ERRDualNV_NewTrackObject = T2ERRDualNV_NewTrack(settings=settings, ifPlotPulse=not(ifLooped)) 
    T2ERRDualNV_NewTrackObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    T2ERRDualNV_NewTrackObject.close()

    if ifRndPhaseNoise:
        AG.disable_PM()
        AG.disable_RFOutput()