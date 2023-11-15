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
from B00_codes.CalibRedROIonizeRate import *
import B00_codes.dataReader as dataReader

####################################################################################################################

reps = 5
for i in np.linspace(1,reps,reps):
    # CalibRedROIonizeRate
    ifRandomized = 0;  ifPlotPulse = (reps==1)
    ifMeaningfulRef=0; ifRefBright=0; ifRefInitAgain = 0
    laserInit_channel = 7; laserRead_channel = 14; ifInitVpz   = 1; ifInitWvl = 0
    if False: 
        velNum = 1; vel_current = 56.5; vel_wvl = 637.22; vel_vpz_target = 76.7
    else:
        velNum = 2; vel_current = 62; vel_wvl = 636.83; vel_vpz_target = 76.56
    
    SRSnum = 2; MWPower = -17; pi_time = 52; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
    MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11

    start = 6e3; stop = 10e3; num_sweep_points = 1
    tausArray = np.linspace(start, stop, num_sweep_points)    
   
    num_loops                    = int(1e5)
    laser_init_delay             = 2e2;      laser_init_duration       = 10e3
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:     0}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)       
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    DAQ_to_laser_off_delay       = 1e2;      MWI_to_switch_delay       = 10 # cannot be between 0 and 10
    read_duration                = 250;      ref_laser_to_read_delay   = 1e12
    delay_between_reads          = 100
    
    num_reads = int(start/(read_duration + delay_between_reads))
    settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points,
                'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,     'pi_time':             pi_time,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,     'read_duration':       read_duration,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'MWI_to_switch_delay': MWI_to_switch_delay,
                'ifRandomized':           ifRandomized,            'ifRefBright':         ifRefBright,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'ifMeaningfulRef':        ifMeaningfulRef,
                'ref_laser_to_read_delay':ref_laser_to_read_delay, 'ifRefInitAgain':      ifRefInitAgain,
                'delay_between_reads': delay_between_reads,
                'num_reads':              num_reads}

    start = time.time()
    CalibRedROIonizeRateObject = CalibRedROIonizeRate(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    CalibRedROIonizeRateObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = CalibRedROIonizeRateObject.getDataFilename()
    # if ifPlotPulse: dataReader.readData(dataFilename)
    CalibRedROIonizeRateObject.close()
    time.sleep(2)

    # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-23/#059_CalibRedROIonizeRate_17-02-21/CalibRedROIonizeRateObject_sig_set.dat'
    # dataReader.readDataFullData(dataFilename)
        




