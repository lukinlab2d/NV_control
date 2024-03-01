"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ODMR_RR import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

####################################################################################################################
reps = 1;  ifLooped = (reps != 1); laserInit_channel = 7
ifNeedVel = 0; ifInitVpz = 0; ifInitWvl = 0

for i in np.linspace(1, reps, reps):
    # Pulsed ODMR_RR
    freqsArray = np.linspace(2922e6, 2942e6, 101)        
    
    if False:
        velNum = 1; vel_current = 56.5; vel_wvl = 637.22
        vel_vpz_target = 69.8; laserRead_channel = 5
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.83
        vel_vpz_target          = 75; laserRead_channel = 14

    if False: 
        SRSnum = 1; MWPower = -18; MWI_duration = 52; MWFreq  = 2747.88e6   #NV D1
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    else:
        SRSnum = 4; MWPower = -9; MWI_duration = 72
        MWI_channel = 0; MWQ_channel = 0; MWswitch_channel = 16
    
    num_loops                    = int(3e4)
    laser_init_delay             = 1e2;        laser_init_duration    = 8e3
    MW_to_read_delay             = 1e2      
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
    read_duration                = 12300;        read_laser_duration    = 12200

    if_tracking = 0; threshold_repumpVpz = 18; threshold_scanVpz = 23; MWFreq = 2747.88e6 
    num_loops_track = 5e3; num_of_cavity_conditioning = 1 # if already cond' for many times, use more avgs. Not that longer
    start = vel_vpz_target - 1.6; stop = vel_vpz_target + 1.6; num_sweep_points = 65; 
    vpzArray = np.linspace(start, stop, num_sweep_points)

    RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
                            'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                            'SRSnum':   SRSnum,
                            'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                            'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        MWI_duration,
                            'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                            'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                            'read_laser_duration':    read_laser_duration,   
                            'MW_to_read_delay':       MW_to_read_delay,   
                            'vel_current':            vel_current,           'vel_wvl':             vel_wvl,
                            'velNum': velNum,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
                            'ifInitWvl':ifInitWvl,
                            'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}
    settings = {'freqsArray': freqsArray, 'num_loops':num_loops, 'MWPower':MWPower, 'SRSnum': SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration': MWI_duration,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_target': vel_vpz_target, 'ifInitVpz':ifInitVpz,
                'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'RRtrackingSettings': RRtrackingSettings, }

    start = time.time()
    ODMR_RRObject = ODMR_RR(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
    ODMR_RRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if ODMR_RRObject.hasTracked:
        vel_vpz_target = ODMR_RRObject.vpz

