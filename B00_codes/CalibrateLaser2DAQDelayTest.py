"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.CalibrateLaser2DAQDelay import *

####################################################################################################################
for i in np.linspace(1000,2000,1):
    # CalibrateLaser2DAQDelay
    start = -100; stop = 300; num_sweep_points = 41
    tausArray = np.linspace(start, stop, num_sweep_points)

    laserRead_channel = 8 # 532 is 3, 589 is 6
    
    if True:
        num_loops        = int(0.5e6)
        laser_init_delay = 5e3;     laser_init_duration = 300
        read_duration    = 50;      #laser_to_DAQ_delay = 260

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops,
                    'tausArray': tausArray,
                    'laser_init_delay': laser_init_delay,'laser_init_duration': laser_init_duration,
                    'laserRead_channel': laserRead_channel,
                    'read_duration':    read_duration}
        
        start = time.time()
        CalibrateLaser2DAQDelayObject = CalibrateLaser2DAQDelay(settings=settings, ifPlotPulse=True, 
                                                              ifRandom=False, ifSweepLaserInitDelay=False) # this is implemented as an Instrument
        CalibrateLaser2DAQDelayObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        CalibrateLaser2DAQDelayObject.close()