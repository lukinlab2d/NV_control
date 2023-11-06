"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
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


for SRSnum in np.linspace(1,2,2):
    # Rabi
    start = 10; stop = 970; num_sweep_points = 2; ifLooped = 1
    tausArray = np.linspace(start, stop, num_sweep_points)
    uwPower = -25; uwFreq = 2870e6
    laserInit_channel = 7; laserRead_channel = 7
    MWswitch_channel = 2; MWI_channel = 1; MWQ_channel = 0
    
    if True:
        print(uwFreq)

        num_loops               = int(1e1)
        laser_init_delay        = 0;        laser_init_duration = 0
        laser_to_MWI_delay      = 1000;       
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)  
        read_duration           = 250
        DAQ_to_laser_off_delay  = 1000;     MWI_to_switch_delay       = 10 # cannot be between 0 and 10

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'MWI_to_switch_delay':       MWI_to_switch_delay,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                    'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                    'SRSnum': SRSnum}

        start = time.time()
        RabiObject = Rabi(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        RabiObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        if not ifLooped: dataFilename = RabiObject.getDataFilename()
        # guess=(0.2, 400, 0, 0.9, 600)
        # dataReader.readData(dataFilename, type='RabiDecay', guess=guess, ifFit=1)
        RabiObject.close()
        
        # guess=(0.3, 250, 0, 0.9)
        # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-06-14/#008_Rabi_16-57-21/RabiObject_sig_set.dat'
        # dataReader.readData(dataFilename, type='Rabi', guess=guess, typeNorm=REF_MINUS_SIG)




