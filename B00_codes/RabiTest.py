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
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
from PlotPulse import *
from Rabi import *
import dataReader

####################################################################################################################


for i in np.linspace(1000,2000,1):
    # Rabi
    start = 250; stop = 10; num_sweep_points = 61
    tausArray = np.linspace(start, stop, num_sweep_points)
    uwPower = -30; uwFreq = 2.8705e9
    if True: #uwFreq != 2.868e9:
        print(uwFreq)

        # Test for pulsed ODMR
        num_loops               = int(1e6)
        laser_init_delay_in_ns  = 10;       laser_init_duration_in_ns = 1e3
        laser_to_AFG_delay      = 1000;       
        laser_to_DAQ_delay      = 1400;     read_duration             = 300
        DAQ_to_laser_off_delay  = 2000;     AFG_2_switch_delay        = 10 # cannot be between 0 and 10

        settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'laser_init_delay_in_ns': laser_init_delay_in_ns,'laser_init_duration_in_ns': laser_init_duration_in_ns,
                    'laser_to_AFG_delay':     laser_to_AFG_delay ,   
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'AFG_2_switch_delay':        AFG_2_switch_delay}

        start = time.time()
        RabiObject = Rabi(settings=settings, ifPlotPulse=True) # this is implemented as an Instrument
        RabiObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')

        dataFilename = RabiObject.getDataFilename()
        dataReader.readData(dataFilename, type='Rabi')
        RabiObject.close()
        




