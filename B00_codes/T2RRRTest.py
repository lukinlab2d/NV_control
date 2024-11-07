"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T2RRR import *
from B00_codes.XY8 import *
import B00_codes.dataReader as dataReader

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2

####################################################################################################################
# T2RRR  
reps = 1;  ifRandomized=0; ifLooped = (reps!=1); ifAWG=1; laserInit_channel = 3;   
ifNeedVel = 0; ifInitVpz = 0; ifInitWvl   = 0
tausArray = np.linspace(4,15004,376)

if False:
    velNum = 1; vel_current = 62.7; vel_wvl = 637.20; laserRead_channel = 5
    SRSnum = 1; MWPower = 0.5; pi_half = 10; MWFreq  = 2597.9e6   #NV D1
    SDGnum = 1; AWG_channel = 18 
    MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
    
else:
    velNum = 2; vel_current = 67; vel_wvl = 636.88; laserRead_channel = 14
    SRSnum = 2; MWPower = -2; pi_half = 10; MWFreq  = 2789.2e6   #NV D2, 2nd MW path
    SDGnum = 2; AWG_channel = 19
    MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11  

# Params
num_loops                    = int(0.7e5)
laser_init_delay             = 1e2;        laser_init_duration    = 8e3
MW_to_read_delay             = 1e2;        MWI_to_switch_delay    = 10
laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1650, 14:900}
laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
laser_to_MW_delay            = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 150
read_duration                = 2e3;        read_laser_duration    = read_duration
AWG_buffer                   = 1;          AWG_output_delay       = 1450     

for i in np.linspace(1,reps,reps):
    
    settings = {'tausArray': tausArray, 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum,
                'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                'laser_to_MW_delay':     laser_to_MW_delay ,   'MWI_to_switch_delay':       MWI_to_switch_delay,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                'MW_to_read_delay':       MW_to_read_delay,      'read_laser_duration':       read_laser_duration,
                'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                'pi_half':        pi_half,
                'vel_current':  vel_current, 'vel_wvl': vel_wvl, 'velNum': velNum,
                'ifInitVpz':ifInitVpz, 'ifInitWvl':ifInitWvl,
                'MWI_channel': MWI_channel,  'MWQ_channel': MWQ_channel,  'MWswitch_channel': MWswitch_channel,
                'ifRandomized': ifRandomized,'ifNeedVel':ifNeedVel,
                'ifAWG': ifAWG, 'SDGnum': SDGnum, 'AWG_channel':AWG_channel, 'AWG_buffer':AWG_buffer, 'AWG_output_delay':AWG_output_delay }

    start = time.time()
    print(MWFreq)
    T2RRRObject = T2RRR(settings=settings, ifPlotPulse=(reps == 1)) # this is implemented as an Instrument
    T2RRRObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    T2RRRObject.close()