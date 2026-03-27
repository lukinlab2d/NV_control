"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.RabiAWG import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 1; ifRandomized=0#; freqs = 
for i in range(int(reps)):
    # RabiAWG
    start = 2; stop = 60; num_sweep_points = 30; ifLooped=(reps != -1)
    tausArray = np.linspace(start, stop, num_sweep_points)
    SDGnum=1; SRSnum=1; uwPower = -24.2+3.5*0; uwFreq = 2692e6; AWG_channel = 18
    # SDGnum=2; SRSnum=2; uwPower = 15; uwFreq = 2704.24e6; AWG_channel = 9
    laserInit_channel = 3; laserRead_channel = 3; ifFakeRabi=0 #just park at tau=tpi

    num_loops               = int(4e5);  AWGbuffer               = 10
    laser_init_delay        = 0;         laser_init_duration     = 0
    laser_to_AWG_delay      = 100;       AWG_output_delay        = 1450      
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)
    read_duration           = 300;       DAQ_to_laser_off_delay  = 2000
    MWI_duration            = 40;        MW_to_read_delay        = 100
    if True:
        if ifFakeRabi>0:
            start = 1; stop = int(ifFakeRabi); num_sweep_points = int(ifFakeRabi)
            tausArray = np.linspace(start, stop, num_sweep_points)
            num_loops = int(1e5)

        settings = {'num_loops':num_loops, 'tausArray':tausArray, 'MWI_duration':MWI_duration,
                    'SRSnum':SRSnum, 'uwPower':uwPower, 'uwFreq': uwFreq, 'SDGnum':SDGnum,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'AWGbuffer':           AWGbuffer,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'AWG_output_delay':    AWG_output_delay,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'AWG_channel':            AWG_channel,           'ifFakeRabi': ifFakeRabi, 
                    'MW_to_read_delay':MW_to_read_delay, 'ifRandomized':ifRandomized}

        start = time.time()
        RabiAWGObject = RabiAWG(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        RabiAWGObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        
        # if not ifLooped: 
        #     dataFilename = RabiAWGObject.getDataFilename()
        #     # dataFilename = 'C:/Users/lukin2dmaterials/data/2024-04-22/#162_RabiAWG_22-51-20/RabiAWGObject_sig_set.dat'
        #     # dataReader.readData(dataFilename)
        #     guess=(0.2, 1000, 0, 0.9, 600)
        #     dataReader.readData(dataFilename, type='RabiDecay', guess=guess, ifFit=0)
        RabiAWGObject.close()