"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T2ERR import *
from B00_codes.XY8 import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
pi = np.pi
####################################################################################################################
reps = 5;  ifRandomized = 0; ifLooped = (reps!=1); normalized_style = Q_FINAL
ifAWG=1;   laserInit_channel = 3;      ifScanVpz = 0;    
vel_vpz_start     = 57;  vel_vpz_end       = 62
vel_vpz_step      = 0.1; vel_vpz_step_time = 0.3; 

# T2ERR   
tausArray = np.round(np.logspace(2,np.log10(2e6),101),-1)

if True:
    velNum    = 1; vel_current = 62.7; vel_wvl   = 637.2; laserRead_channel = 5
    ifNeedVel = 0; ifInitVpz   = 0;    ifInitWvl = 0;     vel_vpz_target    = 78.52;   
else:
    velNum = 2; vel_current = 67; vel_wvl = 636.83; laserRead_channel = 14
    vel_vpz_target    = 79.2; ifInitVpz   = 1
if True: 
    SRSnum = 1; MWPower = 0; pi_half = 10; MWFreq  = 2597.9e6   #NV D1
    SDGnum = 1; AWG_channel = 18 
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
else:
    SRSnum = 2; MWPower = -17; pi_half = 26; MWFreq  = 2785.2e6   #NV D2, 2nd MW path
    MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11

# Params
num_loops                    = int(1e5);   phi_IQ                 = pi/2 #rad, angle of last pi/2
laser_init_delay             = 1e2;        laser_init_duration    = 8e3
MW_to_read_delay             = 1e2;        MWI_to_switch_delay    = 10
laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 2e3;        read_laser_duration    = 2e3
AWG_buffer                   = 1;          AWG_output_delay       = 1450  

if_tracking = 0; threshold_repumpVpz = 20; threshold_scanVpz = 17
num_loops_track = 5e3; fracOfTausArray = 0.6; num_of_cavity_conditioning = 1 # if already cond' for many times, use more avgs. Not that longer
start = vel_vpz_target - 1.6; stop = vel_vpz_target + 1.6; num_sweep_points = 65; 
vpzArray = np.linspace(start, stop, num_sweep_points)

for i in np.linspace(1,reps,reps):
    RRtrackingSettings = {'if_tracking': if_tracking, 'threshold_repumpVpz': threshold_repumpVpz, 'threshold_scanVpz': threshold_scanVpz,
                    'fracOfTausArray': fracOfTausArray,
                'vpzArray': vpzArray,    'num_loops':num_loops_track,  'MWPower':MWPower,    'MWFreq': MWFreq,
                    'SRSnum':   SRSnum,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_duration':        2*pi_half,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':       read_duration,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel,
                    'read_laser_duration':    read_laser_duration,   
                    'MW_to_read_delay':       MW_to_read_delay,   
                    'vel_current':            vel_current,           'vel_wvl':             vel_wvl,
                    'velNum': velNum,                'ifInitVpz':ifInitVpz, 'num_of_cavity_conditioning':num_of_cavity_conditioning,
                    'ifInitWvl':ifInitWvl,
                    'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,}
    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,   'MWI_to_switch_delay':       MWI_to_switch_delay,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'normalized_style':       normalized_style,       'pi_half':        pi_half,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'vel_vpz_start': vel_vpz_start, 'vel_vpz_target': vel_vpz_target, 'vel_vpz_end': vel_vpz_end,
                'vel_vpz_step': vel_vpz_step, 'vel_vpz_step_time': vel_vpz_step_time, 'ifScanVpz': ifScanVpz, 
                'ifInitVpz':ifInitVpz, 'ifInitWvl':ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'RRtrackingSettings': RRtrackingSettings, 'ifRandomized': ifRandomized,'ifNeedVel':ifNeedVel, 'phi_IQ':phi_IQ,
                'ifAWG': ifAWG, 'SDGnum': SDGnum, 'AWG_channel':AWG_channel, 'AWG_buffer':AWG_buffer, 'AWG_output_delay':AWG_output_delay  }

    start = time.time()
    print(MWFreq)
    T2ERRObject = T2ERR(settings=settings, ifPlotPulse=(reps == 1)) # this is implemented as an Instrument
    T2ERRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    if if_tracking:
        vel_vpz_target = T2ERRObject.vpz

    T2ERRObject.close()