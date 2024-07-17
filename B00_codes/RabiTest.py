"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.Rabi import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 1
for i in range(reps):
    # Rabi
    start = 12; stop = 412; num_sweep_points = 51; ifLooped = False
    tausArray = np.linspace(start, stop, num_sweep_points)
    SRSnum=2; uwPower = 0; uwFreq = 2869.2e6
    laserInit_channel = 3; laserRead_channel = 3
    MWswitch_channel = 11; MWI_channel = 12; MWQ_channel = 13
    
    if True:
        print(uwFreq)

        num_loops               = int(1e5)
        laser_init_delay        = 0;        laser_init_duration = 0
        laser_to_MWI_delay      = 1000;       
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)  
        read_duration           = 250
        DAQ_to_laser_off_delay  = 20000;     MWI_to_switch_delay       = 10 # cannot be between 0 and 10

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 
                    'SRSnum':SRSnum, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'MWI_to_switch_delay':       MWI_to_switch_delay,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                    'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel
                    }

        start = time.time()
        RabiObject = Rabi(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        RabiObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        if not ifLooped: 
            dataFilename = RabiObject.getDataFilename()
            guess=(0.2, 1000, 0, 0.9, 600)
            dataReader.readData(dataFilename, type='RabiDecay', guess=guess, ifFit=0)
        RabiObject.close()
        
        # guess=(0.3, 250, 0, 0.9)
        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-06-14/#008_Rabi_16-57-21/RabiObject_sig_set.dat'
        # dataReader.readData(dataFilename, type='Rabi', guess=guess, typeNorm=REF_MINUS_SIG)




