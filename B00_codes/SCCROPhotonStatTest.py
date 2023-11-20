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
from B00_codes.SCCROPhotonStat import *
import B00_codes.dataReader as dataReader

####################################################################################################################
reps=12; sweepWhat = 'ti'; whichNV = 2
ifRandomized = 0; ifLooped = (reps != 1); ifMeaningfulRef=True; 
ifInitWvl = 0; ifInitVpz = 0; laserInit_channel = 7; laserIon_channel = 10; laserShelve_channel = 7
for tr in np.linspace(1e5,2.1e6,10):
    # Test for SCCROPhotonStat
    if whichNV == 1: 
        velNum = 1; vel_current = 56.5; vel_wvl = 637.22; laserRead_channel = 5
        vel_vpz_target = 80.85
        SRSnum = 1; MWPower = -18; MWI_duration = 52; MWFreq  = 2747.88e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    else:
        velNum = 2; vel_current = 49; vel_wvl = 636.83; laserRead_channel = 14
        vel_vpz_target = 77.2
        SRSnum = 2; MWPower = -16.5; MWI_duration = 52; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11

    tausArray = np.array((20,50,80,100,130,160,180,200,230,260,280)) # sweep t_i
    
    num_loops                    = int(1e4); 
    laser_init_delay             = 1e2;      laser_init_duration       = 8e3
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14: 900}
    laserRead_to_DAQ_delay       = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150  
    MWI_to_shelve_delay          = 100;      shelve_duration           = 0
    shelve_to_ion_delay          = 0;        ion_duration              = -1  # tweakable
    ion_to_laserRead_delay       = 5e3;      DAQ_duration              = int(tr)      
    laserRead_signal_duration    = DAQ_duration
    
    settings = {'tausArray': tausArray, 'num_loops':num_loops,
                'MWPower':                MWPower,       'MWFreq':  MWFreq,  'SRSnum': SRSnum,
                'vel_current':            vel_current,   'vel_wvl': vel_wvl, 'velNum': velNum,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,     'MWI_duration':        MWI_duration,
                'MWI_to_shelve_delay':    MWI_to_shelve_delay,     'shelve_duration':     shelve_duration,
                'shelve_to_ion_delay':    shelve_to_ion_delay,     'ion_duration':        ion_duration,
                'ion_to_laserRead_delay': ion_to_laserRead_delay,
                'laserRead_to_DAQ_delay': laserRead_to_DAQ_delay,  'DAQ_duration':        DAQ_duration, 
                'laserRead_signal_duration':laserRead_signal_duration,   
                'ifRandomized':           ifRandomized,            'ifMeaningfulRef':     ifMeaningfulRef,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'laserShelve_channel':    laserShelve_channel,     'laserIon_channel':       laserIon_channel,        
                'sweepWhat':              sweepWhat,
                'ifInitVpz':ifInitVpz,   'ifInitWvl':ifInitWvl, 'vel_vpz_target': vel_vpz_target, 
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}

    start = time.time()
    SCCROPhotonStatObject = SCCROPhotonStat(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    SCCROPhotonStatObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = SCCROPhotonStatObject.getDataFilename()
    SCCROPhotonStatObject.close()


