"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.CalibRedRRIznRateWithIznPulse import *
import B00_codes.dataReader as dataReader

####################################################################################################################
reps = -1
# for ion in np.array((5e3,7.5e3,10e3)):
for iii in np.linspace(1,5,10):
    # CalibRedRRIznRateWithIznPulse
    ifPlotPulse = (reps==1); ifIonizedRef=1; ifIznWithRR=1; ifIznWithMW = 0; ifMWReadLowDutyCycle = 0
    laserInit_channel = 7; laserIon_channel = 10; ifInitVpz   = 0; ifInitWvl = 0;ifNeedVel = 0
    hiLoMWPwr_channel = 17
    if True: 
        velNum = 1; vel_current = 62.2; vel_wvl = 637.22; vel_vpz_target = 72.3
        ifMWDuringRead = 1; SRSnum = 1; MWPower = -2.7; MWFreq  = 2747.88e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2; laserRead_channel = 5

        ifMW2DuringRead = 1; SRSnum2 = 3; MWPower2 = -6; MWFreq2  = 3007.65e6   #NV D1
        MWI2_channel = 1; MWQ2_channel = 0; MWswitch2_channel = 15
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.83; vel_vpz_target = 76.56
        ifMWDuringRead = 1
        SRSnum  = 2; MWPower = -1; MWFreq   = 2838.26e6   #NV D2, ms-1
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11; laserRead_channel = 14
        ifMW2DuringRead = 1
        SRSnum2 = 4; MWPower2 = 10;   MWFreq2  = 2932.8e6   #NV D2 ms+1
        MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16

    num_loops                    = int(3e3)
    laser_init_delay             = 1e3;       laser_init_duration       = int(2.5e6)
    laserSwitch_delay_directory  = {3:850, 6:1150, 9:1150, 7:900, 5:1750, 10:120, 14:  900}
    RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
    laser_to_ion_delay           = laserSwitch_delay_directory.get(laserInit_channel, 0) + 200
    ion_duration                 = 4200;       ion_to_read_delay         = 2.5e6
    DAQ_to_laser_off_delay       = 1e2;       read_duration             = 7500
    delay_between_reads          = 100;       laserRead_to_MWmix        = RRLaserSwitch_delay
    iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
    total_read_duration          = 3e6
    MWmix_duration_short         = -1;       delay_between_MWmix       = -1

    start = total_read_duration; stop = total_read_duration; num_sweep_points = 1
    tausArray = np.linspace(start, stop, num_sweep_points)  

    num_reads = int(start/(read_duration + delay_between_reads)); ifRandomized = 0 
    nMWread   = int(start/(MWmix_duration_short + delay_between_MWmix))
    settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points,
                'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
                'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
                'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
                'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_ion_delay':     laser_to_ion_delay,      'ion_duration': ion_duration,
                'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'read_duration':    read_duration,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'ion_to_read_delay': ion_to_read_delay,
                'laserRead_to_MWmix':     laserRead_to_MWmix,      'iznLaserSwitch_delay':iznLaserSwitch_delay,
                'ifRandomized':           ifRandomized,            'ifIznWithMW': ifIznWithMW,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'laserIon_channel':       laserIon_channel,
                'ifIonizedRef':           ifIonizedRef,            'ifIznWithRR': ifIznWithRR,
                'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,'hiLoMWPwr_channel':hiLoMWPwr_channel,
                'delay_between_reads':    delay_between_reads,
                'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                'num_reads':              num_reads,               'nMWread': nMWread}

    start = time.time()
    CalibRedRRIznRateWithIznPulseObject = CalibRedRRIznRateWithIznPulse(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    CalibRedRRIznRateWithIznPulseObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = CalibRedRRIznRateWithIznPulseObject.getDataFilename()
    CalibRedRRIznRateWithIznPulseObject.close()
    time.sleep(2)

time.sleep(2)
    




