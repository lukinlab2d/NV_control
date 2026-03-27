"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.RabiDEER import *
import B00_codes.dataReader as dataReader
import time

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 1; ifRandomized=0
for i in range(reps):
    # RabiDEER
    start = 410; stop = 500; num_sweep_points = 10; ifLooped=(reps != 1)
    tausArray = np.linspace(start, stop, num_sweep_points)
    SDGnum=1;  SRSnum=1;  uwPower = 5;  uwFreq = 2925e6; AWG_channel = 18
    SDGnum2=2; SRSnum2=2; uwPower2= 2;  uwFreq2= 2880e6; AWG_channel2 = 9
    laserInit_channel = 3; laserRead_channel = 3; ifFakeRabi=0 #just park at tau=tpi

    num_loops               = int(10e5);  AWGbuffer               = 10
    laser_init_delay        = 0;         laser_init_duration     = 0
    laser_to_AWG_delay      = 0;         AWG_output_delay        = 1450      
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)
    read_duration           = 300;       DAQ_to_laser_off_delay  = 4000
    MWI_duration            = 40;        MW_to_read_delay        = 40
    if True:
        if ifFakeRabi>0:
            start = 1; stop = int(ifFakeRabi); num_sweep_points = int(ifFakeRabi)
            tausArray = np.linspace(start, stop, num_sweep_points)
            num_loops = int(1e5)

        settings = {'num_loops':num_loops, 'tausArray':tausArray, 'MWI_duration':MWI_duration,
                    'uwPower2':uwPower2, 'uwFreq2': uwFreq2,
                    'SRSnum':SRSnum, 'uwPower':uwPower, 'uwFreq': uwFreq, 'SDGnum':SDGnum,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'AWGbuffer':           AWGbuffer,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'AWG_output_delay':    AWG_output_delay,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'AWG_channel':            AWG_channel,           'ifFakeRabi': ifFakeRabi, 
                    'AWG_channel2': AWG_channel2, 'SRSnum2': SRSnum2, 'SDGnum2': SDGnum2,
                    'MW_to_read_delay':MW_to_read_delay, 'ifRandomized':ifRandomized}

        start = time.time()
        RabiDEERObject = RabiDEER(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        RabiDEERObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
    
        RabiDEERObject.close()