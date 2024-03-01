"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
import matplotlib as mpl
import matplotlib.pyplot as plt
from B00_codes.ScanRRFreqPumpFree import *  
from B00_codes.PlotPulse import *
import B00_codes.dataReader as dataReader
import dropbox

####################################################################################################################

reps = 1; ifLooped=(reps!=1); timeSleepReadWvl=0; hiLoMWPwr_channel=17
num_of_cavity_conditioning=1; ifInitWvl=0; ifInitVpz=0; ifNeedVel=1

if True:
    # NV1
    SRSnum  = 1; MWPower  = -2.7; pi_time  = 44; MWFreq   = 2747.88e6   #NV D1 ms-1
    SRSnum3 = 3; MWPower3 = -6;   pi_time3 = 44; MWFreq3  = 3007.65e6   #NV D1 ms+1
    MWI_channel  = 1; MWQ_channel  = 0; MWswitch_channel  = 2; MWswitch3_channel = 15
    
    # NV2
    SRSnum2 = 2; MWPower2 = -1;   pi_time2 = 44; MWFreq2  = 2838.26e6   #NV D2, ms-1
    SRSnum4 = 4; MWPower4 = 10;   pi_time4 = 44; MWFreq4  = 2932.8e6    #NV D2 ms+1
    MWI2_channel = 12; MWQ2_channel = 13; MWswitch2_channel = 11; MWswitch4_channel = 16

for i in np.linspace(1,reps,reps):
    # Test for ScanRRFreqPumpFree
    # if np.mod(i,2) == 1: 
    if False:       
        velNum = 1; vel_current = 62.2; vel_wvl = 637.22; vel_vpz_target = 72.3; laserRead_channel = 5

        start = 53; stop = 89; num_sweep_points = 901
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.83; vel_vpz_target = 76.56; laserRead_channel = 14

        start  = 69.4; stop  = 86; num_sweep_points  = 351
        # start  = 72.8; stop  = 86; num_sweep_points  = 326
        start2 = stop+10; stop2 = start2+5; num_sweep_points2 = 126 # fine scan

    vpzArray = np.linspace(start, stop, num_sweep_points)
    vpzArray2 = np.linspace(start2, stop2, num_sweep_points2); vpzArray=np.concatenate((vpzArray,vpzArray2))

    num_loops                    = int(2e4);  
    read_laser_delay             = 1e3;        read_laser_duration    = 100e3
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1750, 14: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    DAQ_duration                 = read_laser_duration-laser_to_DAQ_delay     

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

    settings = {'vpzArray': vpzArray,    'num_loops':num_loops,  'MWPower':MWPower,    'MWFreq': MWFreq,
                'SRSnum':   SRSnum,      'dbx':dbx,
                'read_laser_delay':       read_laser_delay,      'read_laser_duration': read_laser_duration,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'DAQ_duration':       DAQ_duration,
                'laserRead_channel':      laserRead_channel,   
                'vel_current':            vel_current,           'vel_wvl':             vel_wvl, 'velNum': velNum,
                'ifInitVpz':ifInitVpz,   'num_of_cavity_conditioning': num_of_cavity_conditioning,
                'ifInitWvl':ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'timeSleepReadWvl': timeSleepReadWvl,
                
                'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2,
                'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
                'vel_vpz_target': vel_vpz_target,
                'ifNeedVel': ifNeedVel,
                'hiLoMWPwr_channel':      hiLoMWPwr_channel,      
                'MWPower3':MWPower3, 'MWFreq3': MWFreq3, 'SRSnum3': SRSnum3,
                'MWPower4':MWPower4, 'MWFreq4': MWFreq4, 'SRSnum4': SRSnum4,
                'MWswitch3_channel': MWswitch3_channel,'MWswitch4_channel': MWswitch4_channel}

    start = time.time()
    ScanRRFreqPumpFreeObject = ScanRRFreqPumpFree(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    ScanRRFreqPumpFreeObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = ScanRRFreqPumpFreeObject.getDataFilename()
    ScanRRFreqPumpFreeObject.close()