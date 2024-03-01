"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.RabiForScope import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 1
for i in range(reps):
    # Rabi
    start = 80; stop = 80; num_sweep_points = 31000; ifLooped = (reps != 1)
    tausArray = np.linspace(start, stop, num_sweep_points)
    uwPower = -20; uwFreq = 1000e6
    laserInit_channel = 7; laserRead_channel = 7; 
    # MWswitch_channel = 11; MWI_channel = 12; MWQ_channel = 13
    MWswitch_channel = 2; MWI_channel = 1; MWQ_channel = 0
    
    if True:
        print(uwFreq)

        num_loops               = int(1e6)
        laser_init_delay        = 0;        laser_init_duration = 0
        laser_to_MWI_delay      = 10;       
        laser_to_DAQ_delay      = 10
        read_duration           = 10
        DAQ_to_laser_off_delay  = 10;     MWI_to_switch_delay       = 10 # cannot be between 0 and 10

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'MWI_to_switch_delay':       MWI_to_switch_delay,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                    'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                    }

        start = time.time()
        RabiObject = RabiForScope(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        RabiObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        if not ifLooped: 
            dataFilename = RabiObject.getDataFilename()
            guess=(0.2, 1000, 0, 0.9, 600)
            dataReader.readData(dataFilename, type='RabiDecay', guess=guess, ifFit=1)
        RabiObject.close()
        
        # guess=(0.3, 250, 0, 0.9)
        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-06-14/#008_Rabi_16-57-21/RabiObject_sig_set.dat'
        # dataReader.readData(dataFilename, type='Rabi', guess=guess, typeNorm=REF_MINUS_SIG)




