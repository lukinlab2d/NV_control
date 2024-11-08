"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.SCCROPhotonStatBo import *
import B00_codes.dataReader as dataReader

####################################################################################################################
reps=12; sweepWhat = 'ti'; whichNV = 1
ifRandomized = 0; ifLooped = (reps != 1); ifMeaningfulRef=True; 
ifInitWvl = 0; ifInitVpz = 0; laserInit_channel = 7; laserIon_channel = 10; laserShelve_channel = 7
for tr in np.linspace(10e6,30e6,3):
    # Test for SCCROPhotonStatBo
    if whichNV == 1: 
        velNum = 1; vel_current = 50; vel_wvl = 637.22; laserRead_channel = 5
        ifNeedVel = 0; vel_vpz_target = 76.18
        SRSnum = 1; MWPower = -18; MWI_duration = 52; MWFreq  = 2747.88e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    else:
        velNum = 2; vel_current = 49; vel_wvl = 636.83; laserRead_channel = 14
        ifNeedVel = 0; vel_vpz_target = 64.1
        SRSnum = 2; MWPower = -16.5; MWI_duration = 52; MWFreq  = 2838.26e6   #NV D2, 2nd MW path
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11

    tausArray = np.array((20,50,80,100,130,160,180,200,230,260,280,320,360,400,470,550)) # sweep t_i
    
    num_loops                    = int(1e4); 
    laser_init_delay             = 1e2;      laser_init_duration       = 8e3
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14: 900}
    laserRead_to_DAQ_delay       = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150  
    MWI_to_shelve_delay          = 100;      shelve_duration           = 100
    shelve_to_ion_delay          = 600;      ion_duration              = -1  # tweakable
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
                'ifInitVpz':ifInitVpz,   'ifInitWvl':ifInitWvl, 'vel_vpz_target': vel_vpz_target, 'ifNeedVel': ifNeedVel,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}

    start = time.time()
    SCCROPhotonStatBoObject = SCCROPhotonStatBo(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    SCCROPhotonStatBoObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = SCCROPhotonStatBoObject.getDataFilename()
    SCCROPhotonStatBoObject.close()


