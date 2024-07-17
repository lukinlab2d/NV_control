"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalT2R4point import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
# ConfocalT2R4point
for i in range(1):
    xArray = np.linspace(-1.5,0.5,51); yArray = np.linspace(-1,1,51)

    freqsArray = np.linspace(2842e6,2842e6, 1)

    SDGnum=1; SRSnum=1; MWPower = -2
    laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

    num_loops                    = int(0.95e6); settleTime             = 2e9
    laser_init_delay             = 0;        laser_init_duration    = 0
    laser_to_AWG_delay           = 0;        tau                    = 4
    pitime                       = 38;       pitime2                = 38
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    AWG_output_delay             = 1450;     AWGbuffer              = 1
    read_duration                = 300;      DAQ_to_laser_off_delay = 400
    xref_pitime                  = -1.5;     pi_incr_factor         = 11; pi_incr_factor2 = pi_incr_factor#ns/Vx
    read_offset_from_AWG_delay   = 300;      DAQ_error_factor       = 0.97
    

    settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 'tau':tau,
                'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pitime2':pitime2,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                'AWG_output_delay':       AWG_output_delay,      'AWGbuffer': AWGbuffer,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                'AWG_channel':    AWG_channel, 'read_offset_from_AWG_delay': read_offset_from_AWG_delay,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,
                'xArray': xArray, 'yArray':yArray,'settleTime': settleTime,
                'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                'DAQ_error_factor':DAQ_error_factor}

    start = time.time()
    ConfocalT2R4pointObject = ConfocalT2R4point(settings=settings, ifPlotPulse=0) 
    ConfocalT2R4pointObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    ConfocalT2R4pointObject.close()