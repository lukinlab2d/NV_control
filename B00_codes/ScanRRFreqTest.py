"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.ScanRRFreq import *  
from B00_codes.PlotPulse import *
import dropbox

####################################################################################################################
reps = 1; ifLooped = (reps != 1); ifFit = 0; laserInit_channel = 3; 
ifAWG=1; SDGnum=1; AWG_channel = 18
num_of_cavity_conditioning = 1; ifInitWvl = 1; timeSleepReadWvl = 0
for i in np.linspace(1,reps,reps):
    # Test for ScanRRFreq
    if True:
        # velNum = 1; vel_current = 62.7; vel_wvl = 637.2; laserRead_channel = 5
        velNum = 2; vel_current = 67; vel_wvl = 636.87; laserRead_channel = 14

        SRSnum = 1; MWPower = -2; MWI_duration = 58; MWFreq  = 2843.87e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
        start = 0; stop = 100; num_sweep_points = 501 # slow scan
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.83; laserRead_channel = 14

        SRSnum = 2; MWPower = -20;   MWI_duration = 44; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11
        start = 80; stop = 81; num_sweep_points = 11 # slow scan
    vpzArray = np.linspace(start, stop, num_sweep_points)

    num_loops                    = int(5e5);   ifInitVpz = 0
    laser_init_delay             = 1e2;        laser_init_duration    = 1e3
    MW_to_read_delay             = 0
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1750, 14: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 2000;       read_laser_duration    = 2000
    AWGbuffer                    = 10;         AWG_output_delay      = 1450      

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
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                'read_laser_duration':    read_laser_duration,
                'MW_to_read_delay':       MW_to_read_delay,   
                'vel_current':            vel_current,           'vel_wvl':             vel_wvl, 'velNum': velNum,
                'ifInitVpz':ifInitVpz,   'num_of_cavity_conditioning': num_of_cavity_conditioning,
                'ifInitWvl':ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'timeSleepReadWvl': timeSleepReadWvl, 'ifAWG': ifAWG, 'AWGbuffer':AWGbuffer, 'AWG_output_delay':AWG_output_delay,
                'SDGnum': SDGnum, 'AWG_channel':AWG_channel
                }

    start = time.time()
    ScanRRFreqObject = ScanRRFreq(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    ScanRRFreqObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = ScanRRFreqObject.getDataFilename()
    ScanRRFreqObject.close()