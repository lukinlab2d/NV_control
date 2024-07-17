"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalODMRAWG import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
# ConfocalODMRAWG - Pulsed ODMR
for i in range(100):
    xArray = np.linspace(-2,2,41); yArray = np.linspace(0,0,1); xref_pitime = -2  

    start = 2760e6; stop = 2800e6; num_sweep_points = 41
    freqsArray = np.linspace(start, stop, num_sweep_points)
    SDGnum=1; SRSnum=1; MWPower = -2; ifLooped = 1#(reps != 1)
    laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

    num_loops                    = int(1e6); settleTime             = 10e6
    laser_init_delay             = 0;        laser_init_duration    = 0
    pitime                       = 70;       laser_to_AWG_delay     = 0
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    AWG_output_delay             = 1450;     AWGbuffer              = 1
    read_duration                = 300;      DAQ_to_laser_off_delay = 400
    ifSingleGreenRead            = 1;        pi_incr_factor         = 20 #ns/volt_x

    settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 
                'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                'AWG_output_delay':       AWG_output_delay,      'AWGbuffer': AWGbuffer,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                'AWG_channel':    AWG_channel,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay, 'ifSingleGreenRead': ifSingleGreenRead,
                'xArray': xArray, 'yArray':yArray,'settleTime': settleTime,
                'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor}

    start = time.time()
    ConfocalODMRAWGObject = ConfocalODMRAWG(settings=settings, ifPlotPulse=not(ifLooped)) 
    ConfocalODMRAWGObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    ConfocalODMRAWGObject.close()
