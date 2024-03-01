"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.PLSCCRRDualNV import *
import dropbox
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

###########################################################################################e#########################
reps = int(1e5)
# PLSCCRRDualNV
ifRandomized = 0;  ifPlotPulse = (reps==1);  ifMWReadLowDutyCycle = 0; ifNeedVel = 0
ifFancySpinInit = 0; normalized_style = Q_FINAL; ifAntiCorrel = 0; ifSinDetect = 0
ifMWDuringRead = 1; ifMW2DuringRead = 1; if_tracking = 1
laserInit_channel=7; laserIon_channel=10; hiLoMWPwr_channel=17; ifInitVpz=0; ifInitWvl=0

taus1 = np.linspace(301,reps,reps-300); tausArray = taus1

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
num_loops                    = int(2e3);  nSpinInit                 = 0
laser_init_delay             = 1e3;       laser_init_duration       = int(1e6)
laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)   
RRLaser2Switch_delay         = laserSwitch_delay_directory.get(laserRead2_channel, 0)        
laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 200
pi_to_ion_delay              = 50;        ion_duration              = 2700
ion_to_read_delay            = 1.5e6;     DAQ_duration              = 1e6
DAQ_to_laser_off_delay       = 1e2;       laserRead_to_MWmix        = RRLaserSwitch_delay
iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
spinInit_RR_duration         = 15e3;      spinInit_RR_to_pi_delay   = RRLaserSwitch_delay+200
spinInit_pi_to_RR_delay      = 50;        shift_btwn_2NV_read       = DAQ_duration+3e3
MWmix_duration_short         = -1;        delay_between_MWmix       = -1
MWI_to_switch_delay          = 30
######################################### Scan params ########################################
threshold_scanVpz            = 2.5;         threshold_scanVpz2        = 3
num_loops_track              = int(2e4);   
laser_init_delay_track       = 1e2;        laser_init_duration_track  = 20e3
MW_to_read_delay_track       = 1e2
laser_to_DAQ_delay_directory = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:900}
laser_to_MWI_delay_track     = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration_track          = 12300;      read_laser_duration_track  = 12200
#######################################
time_sleep_after_scan        = 45;         wvl_correction             = 8e-6
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
            'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,
            'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    'MWI_to_switch_delay': MWI_to_switch_delay,
            'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
            'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
            'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
            'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'normalized_style':normalized_style,
            'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'shift_btwn_2NV_read':shift_btwn_2NV_read,
            'laserRead2_channel':     laserRead2_channel,
            'MWPower3':MWPower3, 'MWFreq3': MWFreq3, 'SRSnum3': SRSnum3,
            'MWPower4':MWPower4, 'MWFreq4': MWFreq4, 'SRSnum4': SRSnum4,
            'MWswitch3_channel': MWswitch3_channel,'MWswitch4_channel': MWswitch4_channel,
            'RRLaser2Switch_delay':RRLaser2Switch_delay,
            'RRtrackingSettings':RRtrackingSettings, 'RRtrackingSettings2':RRtrackingSettings2,
            'ifAntiCorrel':ifAntiCorrel, 'ifSinDetect':ifSinDetect}

start = time.time()
PLSCCRRDualNVObject = PLSCCRRDualNV(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
PLSCCRRDualNVObject.runScan()
print('Total time = ' + str(time.time() - start) + ' s')

dataFilename = PLSCCRRDualNVObject.getDataFilename()
PLSCCRRDualNVObject.close()