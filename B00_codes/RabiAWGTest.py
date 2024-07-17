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
reps = 1
for i in range(reps):
    # RabiAWG
    start = 2; stop = 122; num_sweep_points = 31; ifLooped =1#(reps != 1)
    tausArray = np.linspace(start, stop, num_sweep_points)
    SDGnum=1; SRSnum=1; uwPower = -2; uwFreq = 2845.62e6; print(uwFreq)
    laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

    num_loops               = int(3e5);  AWGbuffer               = 10
    laser_init_delay        = 0;         laser_init_duration     = 0
    laser_to_AWG_delay      = 0;         AWG_output_delay        = 1450      
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)
    read_duration           = 300;       DAQ_to_laser_off_delay  = 400

    if True:
    
        settings = {'num_loops':num_loops, 'tausArray':tausArray,
                    'SRSnum':SRSnum, 'uwPower':uwPower, 'uwFreq': uwFreq, 'SDGnum':SDGnum,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'AWGbuffer':           AWGbuffer,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'AWG_output_delay':    AWG_output_delay,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'AWG_channel':            AWG_channel}

        start = time.time()
        RabiAWGObject = RabiAWG(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        RabiAWGObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        
        if not ifLooped: 
            dataFilename = RabiAWGObject.getDataFilename()
            # dataFilename = 'C:/Users/lukin2dmaterials/data/2024-04-22/#162_RabiAWG_22-51-20/RabiAWGObject_sig_set.dat'
            # dataReader.readData(dataFilename)
            guess=(0.2, 1000, 0, 0.9, 600)
            dataReader.readData(dataFilename, type='RabiDecay', guess=guess, ifFit=0)
        RabiAWGObject.close()