"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.ScanRRFreq import *  
from B00_codes.PlotPulse import *
import dropbox

####################################################################################################################
reps = 1; ifLooped = (reps != 1); ifFit = 0; laserInit_channel = 3; ifAWG=1; ifIQ=ifAWG
num_of_cavity_conditioning = 1; ifInitWvl = 0; timeSleepReadWvl = 0.5
for i in np.linspace(637, 636, 1):
    # Test for ScanRRFreq
    if False:
        velNum = 1; vel_current = 62.7; vel_wvl = 637.2; laserRead_channel = 5
        SRSnum = 1; MWPower = -5.6; MW_duration = 40; MWFreq  = 2598.1e6   #NV D1
        SDGnum = 1; AWG_channel = 18
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
        start = 14.45; stop = 15.05; num_sweep_points = int(41) # slow scan
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.88; laserRead_channel = 14
        SRSnum = 2; MWPower = -9.6; MW_duration = 40; MWFreq  = 2789.2e6   #NV D2
        SDGnum = 2; AWG_channel = 19
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11
        start = 36.65; stop = 37.25; num_sweep_points = int(40+1) # slow scan
    vpzArray = np.linspace(start, stop, num_sweep_points)

    num_loops                    = int(7e4);   ifInitVpz = 0
    laser_init_delay             = 1e2;        laser_init_duration    = 15e3
    MW_to_read_delay             = 1e2
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1800, 14: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 2e3;        read_laser_duration    = read_duration
    AWG_buffer                   = 10;         AWG_output_delay       = 1450      

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
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MW_duration':        MW_duration,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                'read_laser_duration':    read_laser_duration,
                'MW_to_read_delay':       MW_to_read_delay,   
                'vel_current':            vel_current,           'vel_wvl':             vel_wvl, 'velNum': velNum,
                'ifInitVpz':ifInitVpz,   'num_of_cavity_conditioning': num_of_cavity_conditioning,
                'ifInitWvl':ifInitWvl, 'ifIQ':ifIQ,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'timeSleepReadWvl': timeSleepReadWvl, 'ifAWG': ifAWG, 'AWG_buffer':AWG_buffer, 'AWG_output_delay':AWG_output_delay,
                'SDGnum': SDGnum, 'AWG_channel':AWG_channel
                }

    start = time.time()
    ScanRRFreqObject = ScanRRFreq(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    ScanRRFreqObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')

    dataFilename = ScanRRFreqObject.getDataFilename()
    ScanRRFreqObject.close()