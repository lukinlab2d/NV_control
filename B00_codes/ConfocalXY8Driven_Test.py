"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
import pickle; import dataReader as dr
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalXY8DrivenSDG import *
from B00_codes.ConfocalXY84pointDrivenSDG import *
import time; pi = np.pi
import pickle

def main(NXY8=-1,tau_ns=-5,numreps=-5,amp2='sfsf',phi='sfssfs'):
    if4point = 1 # ConfocalXY8Driven
    reps = int(numreps); ifLooped = (reps!=-1); srate=None; mode = 0 # mode = 0: cos mag, mode = 1: sin mag
    for i in np.linspace(1,reps,reps):
        # xArray = np.linspace(-0.52, -0.44, 3); yArray = np.linspace(-0.40,-0.48,3)
        xArray = np.linspace(-0.52, -0.44, 3); yArray = np.linspace(-1.04,-1.12,3)
        if tau_ns==56:   tausArray = np.linspace(0,0.4,51)
        elif tau_ns>=40: tausArray = np.linspace(0,0.5,51)
        elif tau_ns>=28: tausArray = np.linspace(0,0.8,51)
        elif tau_ns>=16: tausArray = np.linspace(0,1.2,51)
        elif tau_ns>=12: tausArray = np.linspace(0,2.0,51)
        else:            tausArray = np.linspace(0,3.0,51)
        tausArray=tausArray[0:41]; tausArray=tausArray[::2]*1 #41
        if False:
            BACDictFile = (
                'C:/Users/lukin2dmaterials/data/2026-02-01/'
                '#049_ConfocalXY8Driven_21-32-20/XY8-3_BAC_0.4V_9.5K_dict.pkl'
            )
            with open(BACDictFile, 'rb') as f:
                loaded_BACdict = pickle.load(f)
            try:
                Btrue = loaded_BACdict[tau_ns]
            except:
                Btrue = loaded_BACdict[76]
            factor = 8.559066/Btrue

            BACDictFile = (
                'C:/Users/lukin2dmaterials/data/2026-02-13/'
                '#032_ConfocalXY8Driven_10-49-17/XY8-3_BAC_0.4V_17.5K_off_dict_coscoscos.pkl'
            )
            with open(BACDictFile, 'rb') as f:
                loaded_BACdict = pickle.load(f)
            try:
                Btrue = loaded_BACdict[tau_ns]
            except:
                Btrue = loaded_BACdict[78]
            factor2 = 7.91958/Btrue
        else: factor=factor2=1

        SDGnum=1; SRSnum=1; MWPower = -24.2; AWG_channel = 18
        SDGnum2=2; SRSnum2=2; MWPower2 = amp2*factor*factor2; AWG2_channel = 9; hiLoMWPwr_channel=17
        laserInit_channel = 3; laserRead_channel = 3
        sweepWhich                   = 'MWPower2';ifCoSweepPhiIQXY8     = 1
        num_loops                    = int(5e5); settleTime             = 2e9
        numxy8                       = NXY8;     ifIncPower             = 1 # incr pwr istd of pipulse
        laser_init_delay             = 0;        laser_init_duration    = 0
        pitime                       = 24;       pi2time                = pitime/2; tau=tau_ns
        pitime2                      = 24;       BThres = 9e9;  BExt    = 6554
        laser_to_AWG_delay           = 5e3;      phi_IQ_RF              = phi; phi_IQ_XY8 = 'abc'
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
        AWG_output_delay             = 1450;     AWG_buffer             = 10
        read_duration                = 300;      DAQ_to_laser_off_delay = 5e3
        xref_pitime                  = xArray[0];pi_incr_factor         = 2/0.08; pi_incr_factor2=49/0.72 #ns/Vx
        MW_to_read_delay             = 100;      XY8contrastThres       = 10000
        drive_freq                   = 1e9/2/(tau_ns + pitime)
        MWfreqDictFile='C:/Users/lukin2dmaterials/data/2026-03-26/#017_ConfocalODMRAWGFast_14-51-49/MWfreqDict_9.5K.pkl'
        MWfreqDictFilePlus='C:/Users/lukin2dmaterials/data/2026-03-26/#017_ConfocalODMRAWGFast_14-51-49/MWfreqDictPlus_9.5K.pkl'
        phi0DictFile = 'C:/Users/lukin2dmaterials/data/2026-03-11/#021_ConfocalXY84point_15-24-04/phi0Dict_25K.pkl' #XY8-3
        # phi0DictFile = 'C:/Users/lukin2dmaterials/data/2026-03-09/#057_ConfocalXY84point_21-36-26/phi0Dict_25K.pkl' #XY8-1
        
        tausArray2 = np.linspace(1e6,2e6,2)
        with open(phi0DictFile, 'rb') as f:
            phi_0_dict = pickle.load(f)
        XY8DictFile = None#'C:/Users/lukin2dmaterials/data/2025-08-20/#027_ConfocalT2E_13-11-43/T2EDict_42K.pkl'
        settings = {'num_loops':num_loops, 'tausArray':tausArray, 'tausArray2':tausArray2,'srate':srate,
                    'mode':  mode,'numxy8':  numxy8, 'hiLoMWPwr_channel':hiLoMWPwr_channel,
                    'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,'phi_IQ_RF':phi_IQ_RF,'phi_IQ_XY8':phi_IQ_XY8,
                    'XY8contrastThres':XY8contrastThres, 'phi_IQ':phi_IQ_RF,'ifCoSweepPhiIQXY8':ifCoSweepPhiIQXY8,
                    'SRSnum2':SRSnum2, 'MWPower2':MWPower2, 'SDGnum2':SDGnum2,'phi_0_dict':phi_0_dict,
                    'MWfreqDictFile':MWfreqDictFile,'MWfreqDictFilePlus':MWfreqDictFilePlus,'XY8DictFile':XY8DictFile,
                    'BThres':BThres,'BExt':BExt,'sweepWhich':sweepWhich,'drive_freq':drive_freq,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pi2time':pi2time,  'pitime2':pitime2,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'AWG_output_delay':       AWG_output_delay,      'AWG_buffer': AWG_buffer,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                    'DAQ_to_laser_off_delay':DAQ_to_laser_off_delay,
                    'AWG_channel':    AWG_channel,  'MW_to_read_delay': MW_to_read_delay,
                    'AWG2_channel':    AWG2_channel, 'tau':tau,
                    'xArray': xArray, 'yArray':yArray,'settleTime': settleTime, 'ifIncPower':ifIncPower,
                    'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                    }

        start = time.time()
        if if4point==0:
            ConfocalXY8DrivenObject = ConfocalXY8DrivenSDG(settings=settings, ifPlotPulse=not(ifLooped)) 
        else:
            ConfocalXY8DrivenObject = ConfocalXY84pointDrivenSDG(settings=settings, ifPlotPulse=not(ifLooped))
        ConfocalXY8DrivenObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        ConfocalXY8DrivenObject.close()

if __name__ == "__main__":
    taus = np.linspace(8,56,25) # off flake
    for i, tau_ns in enumerate(taus):
        if tau_ns>=0:
            # phi_IQ_RF = phis[i]
            savedDictFile = 'C:/Users/lukin2dmaterials/data/2026-03-07/#080_ConfocalXY84pointDriven_16-09-07/phiOptDict_25K.pkl'
            with open(savedDictFile, 'rb') as f:
                loaded_dict = pickle.load(f)
            phi_IQ_RF = np.round(np.mod(loaded_dict[tau_ns],1),3)*pi
            # for phi_IQ_RF in np.linspace(0.20*pi,0.40*pi,5):
            main(NXY8=3,tau_ns=tau_ns,numreps=1,amp2='dummy',phi=phi_IQ_RF)