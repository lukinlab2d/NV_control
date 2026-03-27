"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalT2RAWGFast import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
# ConfocalT2RAWGFast - Pulsed ODMR
reps = int(100); ifLooped = (reps!=-1); srate=None#2.5e8
for i in np.linspace(1,reps,reps):
    xArray = np.linspace(0,0.4,11); yArray = np.linspace(-0.6,-0.2,11)
    tausArray = np.linspace(20,3620,21); tausArray2 = np.linspace(20,40020,21)

    SDGnum=1; SRSnum=1; MWPower = -8
    laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

    num_loops                    = int(5e5); settleTime             = 2e9
    laser_init_delay             = 0;        laser_init_duration    = 0
    laser_to_AWG_delay           = 100
    pitime                       = 25;       pitime2                = 34
    BThres                       = -100;     BExt                   = 6948
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    AWG_output_delay             = 1450;     AWGbuffer              = 10
    read_duration                = 300;      DAQ_to_laser_off_delay = 4000
    xref_pitime                  = 0;        pi_incr_factor         = 11/0.2; pi_incr_factor2=48/0.8 #ns/Vx
    read_offset_from_AWG_delay   = 300;      DAQ_error_factor       = 0.98
    
    MWfreqDictFile = 'C:/Users/lukin2dmaterials/data/2025-08-09/#034_ConfocalODMRAWGFast_05-11-24/MWfreqDict.pkl'
    MWfreqDictFilePlus = 'C:/Users/lukin2dmaterials/data/2025-08-09/#034_ConfocalODMRAWGFast_05-11-24/MWfreqDictPlus.pkl'

    settings = {'num_loops':num_loops, 'tausArray':tausArray, 'tausArray2':tausArray2,'srate':srate,
                'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,
                'MWfreqDictFile':MWfreqDictFile,'MWfreqDictFilePlus':MWfreqDictFilePlus,
                'BThres':BThres,'BExt':BExt,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pitime2':pitime2,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                'AWG_output_delay':       AWG_output_delay,      'AWGbuffer': AWGbuffer,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                'AWG_channel':    AWG_channel, 'read_offset_from_AWG_delay': read_offset_from_AWG_delay,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,
                'xArray': xArray, 'yArray':yArray,'settleTime': settleTime,
                'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                'DAQ_error_factor':DAQ_error_factor}

    start = time.time()
    ConfocalT2RAWGFastObject = ConfocalT2RAWGFast(settings=settings, ifPlotPulse=0) 
    ConfocalT2RAWGFastObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    ConfocalT2RAWGFastObject.close()