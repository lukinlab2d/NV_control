"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalODMRAWGFast import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
# ConfocalODMRAWGFast - Pulsed ODMR
for i in range(1):
    xArray = np.linspace(-0.2,0.6,21); yArray = np.linspace(-0.08,0.72,21); ifDifferential = 0
    
    a1 = 2834e6; b1 = 2854e6; num_sweep_points = 21
    f1 = np.linspace(a1, b1, num_sweep_points)
    a2 = 2912e6; b2 = 2932e6; num_sweep_points = 21
    f2 = np.linspace(a2, b2, num_sweep_points)
    f3 = np.linspace(2795e6,a1-3e6,4); f4 = np.linspace(b2+3e6,2975e6,4)
    f5 = np.linspace(b1+3e6,a2-3e6, 4)
    freqsArray = np.concatenate((f3,f1,f5,f2,f4))

    # freqsArray = np.linspace(2745e6,3015e6,2)

    SDGnum=1; SRSnum=1; MWPower = -2
    laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

    num_loops                    = int(2e5); settleTime             = 2e9
    laser_init_delay             = 0;        laser_init_duration    = 0
    pitime                       = 68;       laser_to_AWG_delay     = 0; pitime2=0
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    AWG_output_delay             = 1450;     AWGbuffer              = 1
    read_duration                = 300;      DAQ_to_laser_off_delay = 400
    xref_pitime                  = -0.20;    pi_incr_factor         = 13; pi_incr_factor2=0 #ns/Vx
    padding                      = 900;      MW_to_DAQ_delay = 0
    padding_green1               = 100;      AWG_delay = 1800
    DAQ_error_factor             = 0.97;     freqsArray2 = None

    if ifDifferential:
        xArray2 = -xArray
        interleaved_array = np.empty(xArray.size + xArray2.size, dtype=xArray.dtype)
        interleaved_array[0::2] = xArray; interleaved_array[1::2] = xArray2
        xArray = interleaved_array

    settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 'freqsArray2':freqsArray2, 
                'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pitime2':pitime2,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                'AWG_output_delay':       AWG_output_delay,      'AWGbuffer': AWGbuffer,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                'AWG_channel':    AWG_channel, 'AWG_delay': AWG_delay, 'MW_to_DAQ_delay': MW_to_DAQ_delay,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,
                'xArray': xArray, 'yArray':yArray,'settleTime': settleTime,
                'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                'padding':padding,'padding_green1':padding_green1,'ifDifferential':ifDifferential,
                'DAQ_error_factor':DAQ_error_factor}

    start = time.time()
    ConfocalODMRAWGFastObject = ConfocalODMRAWGFast(settings=settings, ifPlotPulse=0) 
    ConfocalODMRAWGFastObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    ConfocalODMRAWGFastObject.close()
    

    # for x in xArray:
    #     # smallXArray = xArray[j*10:(j+1)*10] 
    #     # if j==19: smallXArray = xArray[190:201]   
    #     smallXArray = np.linspace(x,x,1)

    # smallXArray = xArray