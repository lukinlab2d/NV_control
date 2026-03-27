"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalT1 import *
import time

def main():
    NO_MS_EQUALS_1 = 0
    Q_FINAL = 1
    THREE_PI_HALF_FINAL = 2
    REF_MINUS_SIG  = 3; pi = np.pi

    ####################################################################################################################
    # ConfocalT1
    reps = int(5); ifLooped = (reps!=-1); srate=None#2.5e8
    for i in np.linspace(1,reps,reps):
        xArray = np.linspace(-0.24,-0.16,3); yArray = np.linspace(-0.56, -0.64, 1)
        tausArray  = np.round(np.logspace(np.log10(40),np.log10(2e6),16),-1)
        tausArray2 = np.round(np.logspace(np.log10(40),np.log10(12e6),16),-1) # off flake

        SDGnum=1; SRSnum=1; MWPower = -8
        laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

        num_loops                    = int(1e5); settleTime             = 2e9
        laser_init_delay             = 0;        laser_init_duration    = 0
        pitime                       = 23;       ifShimon               = 0 
        pitime2                      = 32;       BThres = 1000000; BExt = 6948
        laser_to_AWG_delay           = 5e3;      sig_to_ref_delay_Shimon= 1000
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
        AWG_output_delay             = 1450;     AWG_buffer             = 10
        read_duration                = 300;      DAQ_to_laser_off_delay = 6e3
        xref_pitime                  = -0.24;    pi_incr_factor         = 9/0.12; pi_incr_factor2=49/0.72 #ns/Vx
        MW_to_read_delay             = 100;      AWG_delay              = AWG_output_delay

        MWfreqDictFile='C:/Users/lukin2dmaterials/data/2025-10-05/#010_ConfocalODMRAWGFast_18-51-52/MWfreqDict_9.5K.pkl'
        MWfreqDictFilePlus='C:/Users/lukin2dmaterials/data/2025-10-05/#010_ConfocalODMRAWGFast_18-51-52/MWfreqDictPlus_9.5K.pkl'
        
        settings = {'num_loops':num_loops, 'tausArray':tausArray, 'tausArray2':tausArray2,'srate':srate,
                    'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,
                    'MWfreqDictFile':MWfreqDictFile,'MWfreqDictFilePlus':MWfreqDictFilePlus,
                    'BThres':BThres,'BExt':BExt,'ifShimon':ifShimon, 'sig_to_ref_delay_Shimon':sig_to_ref_delay_Shimon,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pitime2':pitime2,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'AWG_output_delay':       AWG_output_delay,      'AWG_buffer': AWG_buffer,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                    'DAQ_to_laser_off_delay':DAQ_to_laser_off_delay,
                    'AWG_channel':    AWG_channel, 'AWG_delay': AWG_delay, 'MW_to_read_delay': MW_to_read_delay,
                    'xArray': xArray, 'yArray':yArray,'settleTime': settleTime,
                    'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                    }

        start = time.time()
        ConfocalT1Object = ConfocalT1(settings=settings, ifPlotPulse=not(ifLooped)) 
        ConfocalT1Object.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        ConfocalT1Object.close()
        