"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalT2E import *
import time; pi = np.pi

def main():
    # ConfocalT2E
    reps = int(6); ifLooped = (reps!=-1); srate=None#2.5e8
    for i in np.linspace(1,reps,reps):
        xArray = np.linspace(-0.28,-0.20,3); yArray = np.linspace(-0.12, -0.20, 3)
        # tausArray = np.linspace(20,1020,26); tausArray2 = np.linspace(20,30020,13)
        # tausArray = np.round(np.logspace(np.log10(20),np.log10(8e3),26),-1)
        tausArray2 = np.round(np.logspace(np.log10(20),np.log10(60e3),26),-1)
        tausArray = np.concatenate((np.linspace(20,156,35),np.round(np.logspace(np.log10(160),np.log10(4e3),12),-1)))
        # tausArray = np.concatenate((np.linspace(20,200,10),np.round(np.logspace(np.log10(250),np.log10(4e3),12),-1)))
        # tausArray = np.linspace(400,400,1); tausArray2 = tausArray

        SDGnum=1; SRSnum=1; MWPower = -24.2
        laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

        num_loops                    = int(5e5); settleTime             = 2e9
        laser_init_delay             = 0;        laser_init_duration    = 0
        pitime                       = 24;       pi2time                = pitime/2
        pitime2                      = 24;       BThres = 9e9;  BExt    = 6097
        laser_to_AWG_delay           = 5e3;      phi_IQ                 = pi/2
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
        AWG_output_delay             = 1450;     AWG_buffer             = 10
        read_duration                = 300;      DAQ_to_laser_off_delay = 5e3
        xref_pitime                  = -0.28;    pi_incr_factor         = 6/0.12; pi_incr_factor2=49/0.72 #ns/Vx
        MW_to_read_delay             = 100;      AWG_delay              = AWG_output_delay
        T2EcontrastThres             = 10000;    DAQ_error_factor       = 0.98

        MWfreqDictFile='C:/Users/lukin2dmaterials/data/2026-02-10/#023_ConfocalODMRAWGFast_15-51-07/MWfreqDict_16K.pkl'
        MWfreqDictFilePlus='C:/Users/lukin2dmaterials/data/2026-02-10/#023_ConfocalODMRAWGFast_15-51-07/MWfreqDictPlus_16K.pkl'
        T2EDictFile = None#'C:/Users/lukin2dmaterials/data/2025-08-20/#027_ConfocalT2E_13-11-43/T2EDict_42K.pkl'

        settings = {'num_loops':num_loops, 'tausArray':tausArray, 'tausArray2':tausArray2,'srate':srate,
                    'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,'phi_IQ':phi_IQ,'T2EcontrastThres':T2EcontrastThres,
                    'MWfreqDictFile':MWfreqDictFile,'MWfreqDictFilePlus':MWfreqDictFilePlus,'T2EDictFile':T2EDictFile,
                    'BThres':BThres,'BExt':BExt,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pi2time':pi2time,  'pitime2':pitime2,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'AWG_output_delay':       AWG_output_delay,      'AWG_buffer': AWG_buffer,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                    'DAQ_to_laser_off_delay':DAQ_to_laser_off_delay,
                    'AWG_channel':    AWG_channel, 'AWG_delay': AWG_delay, 'MW_to_read_delay': MW_to_read_delay,
                    'xArray': xArray, 'yArray':yArray,'settleTime': settleTime, 'DAQ_error_factor': DAQ_error_factor,
                    'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                    }

        start = time.time()
        ConfocalT2EObject = ConfocalT2E(settings=settings, ifPlotPulse=not(ifLooped)) 
        ConfocalT2EObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        ConfocalT2EObject.close()

if __name__ == "__main__":
    main()