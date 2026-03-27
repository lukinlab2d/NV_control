"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.T2EDEERNeumann import *
import B00_codes.dataReader as dataReader
import time

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2; pi = np.pi

####################################################################################################################
for ii in range(5):
    reps=int(11); ifLooped=(reps!=1); srate=2.5e8; ifRandomized=1; normalized_style=Q_FINAL
    uwPower  = 5; uwFreq  = 2925e6
    uwPower2 = 2; pitimeDEER = 80 #uwFreq2 = 2910e6; 
    freqs = np.linspace(3052e6,3072e6,reps)
    for i in range(reps):
        uwFreq2 = freqs[i]; print(uwFreq2)
        tausArray = np.linspace(0,3000,31)

        # Params for T2EDEERNeumann
        laserInit_channel            = 3;        laserRead_channel   = 3
        num_loops                    = int(1e6); phi_IQ              = pi/2 #rad, angle of last pi/2
        laser_init_delay             = 0;        laser_init_duration = 0
        pi2time                      = 48;       pitime              = 200
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0) 
        laser_to_AWG_delay           = 0;        read_duration       = 300
        DAQ_to_laser_off_delay       = 4000;     MW_to_read_delay    = 40
        AWG_output_delay             = 1450;     tau_echo            = 4200 

        AWG_channel  = 18; SRSnum  = 1; SDGnum  = 1; AWG_buffer = 10
        AWG_channel2 = 9;  SRSnum2 = 2; SDGnum2 = 2

        if True: 
            if_tracking = 0 # 2 is for the monty setup
            if np.mod(i,5)==4: if_tracking = 1
            laserTrack_channel     = 3;       
            xy_scan_read_time      = 10;       xy_scan_settle_time    = 5;  
            xy_scan_resolution_hor = 20;       xy_scan_resolution_ver = 20
            x_minus_range          = 0.015;    x_plus_range           = 0.015
            y_minus_range          = 0.015;    y_plus_range           = 0.015
            xy_displacement_limit  = 0.015;    num_of_scans           = 2
            time_btwn_trackings    = 10*60 # min*seconds

            xz_scan_resolution_hor = 20;       xz_scan_resolution_ver = 25
            x_minus_range          = 0.015;    x_plus_range           = 0.015
            z_minus_range          = 0.1;      z_plus_range           = 0.1
            xz_displacement_limit  = 0.1; 

        trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
                            'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
                            'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                            'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
                            'xy_displacement_limit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
                            'tracking_period':        1e9,                   'if_tracking':            if_tracking,
                            'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                            'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                            'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                            'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}

        settings = {'num_loops':num_loops, 'tausArray': tausArray, 'uwPower':uwPower, 'uwFreq': uwFreq,
                    'uwPower2':uwPower2, 'uwFreq2': uwFreq2,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration':       laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay ,   'pi2time':    pi2time, 'pitime':  pitime,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay ,   'read_duration':             read_duration,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,'trackingSettings':          trackingSettings,
                    'AWG_output_delay':    AWG_output_delay,   'ifRandomized':              ifRandomized,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':         laserRead_channel,
                    'AWG_channel': AWG_channel, 'SRSnum': SRSnum, 'SDGnum': SDGnum, 'AWG_buffer': AWG_buffer,
                    'AWG_channel2': AWG_channel2, 'SRSnum2': SRSnum2, 'SDGnum2': SDGnum2, 'tau_echo':tau_echo,
                    'phi_IQ': phi_IQ, 'srate':srate,'MW_to_read_delay':MW_to_read_delay,
                    'pitimeDEER':pitimeDEER}
        
        start = time.time()
        T2EDEERNeumannObject = T2EDEERNeumann(settings=settings, ifPlotPulse=not(ifLooped)) # this is implemented as an Instrument
        T2EDEERNeumannObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        T2EDEERNeumannObject.close()