"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalXY8 import *
from B00_codes.ConfocalXY84point import *
import time; pi = np.pi; import pickle

def main(NXY8=-1,numreps=-5,phi_IQ_XY8='a',tau='b'):
    if4point = 1; ifRandomized = 1 # ConfocalXY8
    reps = int(numreps); ifLooped = (reps!=-1); srate=None; mode = 0 # mode = 0: cos mag, mode = 1: sin mag
    for i in np.linspace(1,reps,reps):
        xArray = np.linspace(-0.52, -0.44, 3); yArray = np.linspace(-0.40,-0.48,3)
        # xArray = np.linspace(-0.52, -0.44, 3); yArray = np.linspace(-1.04,-1.12,3)
        tausArray = np.linspace(8,80,37);  tausArray2 = np.linspace(20,2020,51)
        # tausArray = np.linspace(0.1,1.7,26)*pi;  tausArray2 = np.linspace(20,2020,51)
        # tausArray = np.concatenate((np.linspace(20,156,35),np.round(np.logspace(np.log10(160),np.log10(5e2),6),-1)))
        # tausArray = np.array([24,]);  tausArray2 = np.linspace(20,2020,51)

        SDGnum=1; SRSnum=1; MWPower = -24.7
        laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18
        ifCoSweepPhiIQ               = 1;        sweepWhich             = 'tau'
        num_loops                    = int(2e5); settleTime             = 2e9
        numxy8                       = NXY8;     ifIncPower             = 1 # incr pwr istd of pipulse
        laser_init_delay             = 0;        laser_init_duration    = 0
        pitime                       = 24;       pi2time                = pitime/2
        pitime2                      = 24;       BThres = 9e9;  BExt    = 6554
        laser_to_AWG_delay           = 5e3;      #phi_IQ_XY8             = pi/2
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
        AWG_output_delay             = 1450;     AWG_buffer             = 10
        read_duration                = 300;      DAQ_to_laser_off_delay = 5e3
        xref_pitime                  = xArray[0];pi_incr_factor         = 2/0.08; pi_incr_factor2=49/0.72 #ns/Vx
        MW_to_read_delay             = 100;      AWG_delay              = AWG_output_delay; XY8contrastThres=10000
        baseFolder = 'C:/Users/lukin2dmaterials/data/'
        MWfreqDictFile = baseFolder+'2026-03-26/#017_ConfocalODMRAWGFast_14-51-49/MWfreqDict_9.5K.pkl'
        MWfreqDictFilePlus = baseFolder+'2026-03-26/#017_ConfocalODMRAWGFast_14-51-49/MWfreqDictPlus_9.5K.pkl'
        if NXY8==3:
            phi0DictFile = baseFolder+'2026-03-26/#040_ConfocalXY84point_23-16-51/phi0Dict_9.5K_shifted.pkl' #XY8-3
        else:
            phi0DictFile = baseFolder+'2026-03-26/#003_ConfocalXY84point_02-25-31/phi0Dict_9.5K_phizeros.pkl' #XY8-1
            # phi0DictFile = baseFolder+'2026-03-22/#004_ConfocalXY84point_03-09-46/phi0Dict_25K_plus0.025.pkl' #XY8-1
            # phi0DictFile = baseFolder+'2026-03-21/#009_ConfocalXY84point_03-39-08/phi0Dict_25K.pkl' #XY8-1
            # phi0DictFile = baseFolder+'2026-03-14/#015_ConfocalXY84point_03-43-31/phi0Dict_25K.pkl' #XY8-1
            # phi0DictFile = baseFolder+'2026-03-09/#057_ConfocalXY84point_21-36-26/phi0Dict_25K.pkl' #XY8-1
        with open(phi0DictFile, 'rb') as f:
            phi_0_dict = pickle.load(f)
        XY8DictFile = None #baseFolder+'2025-08-20/#027_ConfocalT2E_13-11-43/T2EDict_42K.pkl'
        settings = {'num_loops':num_loops, 'tausArray':tausArray, 'tausArray2':tausArray2,'srate':srate,
                    'mode':  mode,'numxy8':  numxy8, 'ifCoSweepPhiIQ':ifCoSweepPhiIQ, 'phi_0_dict':phi_0_dict,
                    'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,'phi_IQ':phi_IQ_XY8,'XY8contrastThres':XY8contrastThres,
                    'MWfreqDictFile':MWfreqDictFile,'MWfreqDictFilePlus':MWfreqDictFilePlus,'XY8DictFile':XY8DictFile,
                    'BThres':BThres,'BExt':BExt,'tau':tau,'sweepWhich':sweepWhich,'ifRandomized':ifRandomized,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pi2time':pi2time,  'pitime2':pitime2,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'AWG_output_delay':       AWG_output_delay,      'AWG_buffer': AWG_buffer,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                    'DAQ_to_laser_off_delay':DAQ_to_laser_off_delay,
                    'AWG_channel':    AWG_channel, 'AWG_delay': AWG_delay, 'MW_to_read_delay': MW_to_read_delay,
                    'xArray': xArray, 'yArray':yArray,'settleTime': settleTime, 'ifIncPower':ifIncPower,
                    'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                    }

        start = time.time()
        if if4point==0:
            ConfocalXY8Object = ConfocalXY8(settings=settings, ifPlotPulse=not(ifLooped)) 
        else:
            ConfocalXY8Object = ConfocalXY84point(settings=settings, ifPlotPulse=not(ifLooped)) 
        ConfocalXY8Object.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        ConfocalXY8Object.close()

if __name__ == "__main__":
    # taus = np.linspace(8,80,19)
    # for tau in taus:
    main(NXY8=3,numreps=10,phi_IQ_XY8='dakdhskd',tau='tau')