"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.CalibRedRRIonizeRate import *

####################################################################################################################
reps = 1
for iii in np.linspace(1,reps,reps):
    # CalibRedRRIonizeRate
    ifRandomized=0;  ifPlotPulse=(reps==1); ifAWG=0; ifMWDuringRead=0; ifMW2DuringRead=0
    laserInit_channel=3; ifInitVpz=0; ifInitWvl=0; ifNeedVel=0
    if False: 
        velNum = 1; vel_current = 62.7; vel_wvl = 637.20; vel_vpz_target = -1; laserRead_channel = 5
        SRSnum = 1;  MWPower = -20; pi_time  = 20; MWFreq   = 2598.1e6   #NV D1 ms-1
        SDGnum  = 1; AWG_channel = 18
        MWI_channel = 1; MWQ_channel = 0; MWswitch_channel = 2
        SRSnum2 = 3; MWPower2 = -14; pi_time2 = 62; MWFreq2 = 3162e6  #NV D1 ms+1
        MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 15
    else:
        velNum = 2; vel_current = 67; vel_wvl = 636.88; vel_vpz_target = -1; laserRead_channel = 14
        SRSnum = 2; MWPower = -30; pi_time  = 20; MWFreq   = 2789.2e6   #NV D2, ms-1
        SDGnum = 2; AWG_channel = 19
        MWI_channel = 12; MWQ_channel = 13; MWswitch_channel = 11; 
        SRSnum2 = 4; MWPower2 = -4; pi_time2 = 74; MWFreq2 = 3037.2e6  #NV D2 ms+1 #power=10
        MWI2_channel = 0; MWQ2_channel = 0; MWswitch2_channel = 16
    
    tr=50e3
    start = tr; stop = tr; num_sweep_points = 1
    tausArray = np.linspace(start, stop, num_sweep_points)    

    num_loops                    = int(2e4)
    laser_init_delay             = 1e2;     laser_init_duration       = 9e4
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1750, 14:   900}
    laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)       
    laser_to_MWI_delay           = laser_to_DAQ_delay_directory.get(laserInit_channel, 0) + 1e3
    MW_to_read_delay             = 0;       DAQ_to_laser_off_delay    = 1e2
    read_duration                = 400
    delay_between_reads          = 300;     laserRead_to_MWmix        = laser_to_DAQ_delay
    AWG_buffer                   = 40;      AWG_output_delay          = 1450  
    
    num_reads = int(start/(read_duration + delay_between_reads))
    settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points,
                'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
                'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
                'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
                'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,     'pi_time':             pi_time,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,     'read_duration':       read_duration,
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'MW_to_read_delay': MW_to_read_delay,
                'laserRead_to_MWmix':     laserRead_to_MWmix,
                'ifRandomized':           ifRandomized,            
                'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                'delay_between_reads': delay_between_reads,
                'num_reads':              num_reads,
                'ifAWG': ifAWG,
                'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,}

    start = time.time()
    CalibRedRRIonizeRateObject = CalibRedRRIonizeRate(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    CalibRedRRIonizeRateObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = CalibRedRRIonizeRateObject.getDataFilename()
    # if ifPlotPulse: dataReader.readData(dataFilename)
    CalibRedRRIonizeRateObject.close()
    time.sleep(2)

    # time.sleep(5)
        




