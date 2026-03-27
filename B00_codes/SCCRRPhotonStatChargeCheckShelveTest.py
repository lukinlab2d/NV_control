"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.SCCRRPhotonStatChargeCheckShelve import *
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A
import dropbox, time

####################################################################################################################
# Below are params for the first NV pair
reps = 1; ifPlotPulse=(reps==-1)
for tsh in np.linspace(350,0,1):
    for whichNV in np.linspace(1,2,1):
        # SCCRRPhotonStatChargeCheckShelve
        ifNeedVel=0; ifInitVpz=0; ifInitWvl=0; ifFancySpinInit=0; ifGreenLow=1
        ifAWG=1; ifMWDuringRead=1; ifMW2DuringRead=1; ifMWReadLowDutyCycle=0; ifIQ=ifAWG
        ifRandomized = 0; ifHiloExtra=1; ifIonInit=0; ifIonRR=0
        laserInit_channel = 3; laserIon_channel = 8; hiLoMWPwr_channel = 17
        laserShelve_channel=3; laserGreenLow_channel=19; sweepWhich = 'ti'
        ifTestSig=0; ifRndPhaseNoise=0; AGBW = 25e3; AGfreq = 2.5e6; AGamp = 0.2 # beware of heating!!
        if np.mod(whichNV,2)==1: 
            velNum = 1; vel_current = 62.7; vel_wvl = 637.25; vel_vpz_target=-1; laserRead_channel = 5
            velNum = 2; vel_current = 67.1; vel_wvl = 636.88; vel_vpz_target=-1; laserRead_channel = 14
            SRSnum = 1; MWPower  = 0;  pi_time  = 70; MWFreq  = 3858e6   #NV D1 ms-1 2598.19e6
            SRSnum2= 3; MWPower2 = -2; pi_time2 = -1; MWFreq2 = 1906e6   #NV D1 ms+1 3161.27e6
            MWI_channel  = 1; MWQ_channel  = 0;  MWswitch_channel  = 2; ion_duration= 120 #8e2
            MWI2_channel = 0; MWQ2_channel = 0;  MWswitch2_channel = 15; amp_MW_mix = 1
            SDGnum       = 1; AWG_channel  = 18; srate             = 2.5e8
        else:
            velNum = 2; vel_current = 67.1; vel_wvl = 636.88; vel_vpz_target=-1; laserRead_channel = 14
            SRSnum = 2; MWPower  = -9.5; pi_time  = 44; MWFreq = 2788.70e6   #NV D2, ms-1 -9.5 44ns -8.8
            SRSnum2= 4; MWPower2 = 0.00; pi_time2 = -1; MWFreq2= 3037.20e6   #NV D2 ms+1 #power=10
            MWI_channel  = 12; MWQ_channel  = 13; MWswitch_channel  = 11; ion_duration= 4e3
            MWI2_channel = 0;  MWQ2_channel = 0;  MWswitch2_channel = 16; amp_MW_mix = 0.7
            SDGnum       = 2;  AWG_channel  = 19; srate             = 2.5e8
        
        tausArray = np.array((10,20,30,40,50,60,80,100))#,120,140,160,180))#,200))#,400,1000,2e3,3e3,5e3,7e3,10e3,15e3,20e3)) # sweep ti
        # tausArray = np.array((1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6))#,10e6,12e6,15e6,20e6,25e6,30e6)) # sweep tr
        # tausArray = np.linspace(4,160,40)         # rabi
        # tausArray = np.array((5e4,1e5,2e5,5e5,1e6,2e6,3e6,4e6,5e6))#15e6,20e6,30e6))#40e3,200e3,500e3,1e6,2e6,3e6,4e6,5e6,6.5e6,8e6,10e6))#,10e6,12e6)) # tinit
        # tausArray = np.array((0,10,20,30,40,50,60,80,100,120,140,160,180,200,220,240,260,300,330,360,400,430,460,500)) # tsh
        # tausArray = np.linspace(-600,200,41) # sh2i
        # tausArray = np.linspace(-12,0,7) # amp_MW_mix / MWPower2
       
        num_loops                    = int(2e3);   nSpinInit                 = 0
        laser_init_delay             = 1e3;        laser_init_duration       = int(4e6)
        check_duration               = 2e5
        laserSwitch_delay_directory  = {3:860, 6:1160, 9:1160, 7:900, 5:1700, 10:160, 8:80, 14:800}
        RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
        laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 500
        pi_to_shelve_delay           = 160;        shelve_duration           = 260
        shelve_to_ion_delay          = -300
        ion_to_read_delay            = 3e5;        DAQ_duration              = 4e6
        DAQ_to_laser_off_delay       = 1e2;        laserRead_to_MWmix        = 0
        iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
        greenLaserSwitch_delay       = laserSwitch_delay_directory.get(laserShelve_channel,0)
        AWG_buffer                   = 40;         AWG_output_delay          = 1450
        MWmix_duration_short         = 0;          delay_between_MWmix       = 0
        spinInit_RR_duration         = 0;          spinInit_RR_to_pi_delay   = 0*(RRLaserSwitch_delay+100)
        spinInit_pi_to_RR_delay      = 0;          laser_ionInit_duration    = 10
        hilo_margin_start = 40; hilo_margin_end = 40; hilo_min = 60

        # Create a Dropbox client
        tkFile = 'C:/Users/lukin2dmaterials/data/tk.txt'; lines = []
        with open(tkFile, 'r') as file:
            for line in file:
                if "\n" in line: line = line[0:-1]
                lines.append(line)
        access_token = lines[0]
        app_key = lines[1]
        app_secret = lines[2]
        refresh_token = lines[3]
        #######################################
        dbx = dropbox.Dropbox(oauth2_access_token=access_token,
                            app_key=app_key, 
                            app_secret=app_secret, 
                            oauth2_refresh_token=refresh_token)
        
        settings = {'tausArray': tausArray, 'dbx':dbx, 'greenLaserSwitch_delay':greenLaserSwitch_delay,'ifIonRR':ifIonRR,
                    'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
                    'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
                    'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
                    'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
                    'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
                    'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
                    'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                    'laser_to_pi_delay':      laser_to_pi_delay,       'DAQ_duration': DAQ_duration,
                    'check_duration':         check_duration,
                    'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'pi_time': pi_time, 'pi_time2':pi_time2,
                    'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'ion_to_read_delay': ion_to_read_delay,
                    'laserRead_to_MWmix':     laserRead_to_MWmix,      'iznLaserSwitch_delay':iznLaserSwitch_delay,
                    'ifRandomized':           ifRandomized,            'pi_to_shelve_delay': pi_to_shelve_delay,
                    'shelve_to_ion_delay':    shelve_to_ion_delay,     'shelve_duration': shelve_duration,
                    'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                    'laserIon_channel':       laserIon_channel,        'laserShelve_channel':laserShelve_channel,
                    'ion_duration': ion_duration, 'laserGreenLow_channel':laserGreenLow_channel,
                    'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    'srate': int(srate),
                    'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                    'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                    'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                    'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
                    'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'laser_ionInit_duration':laser_ionInit_duration,
                    'ifRndPhaseNoise':ifRndPhaseNoise, 'AGBW':AGBW, 'AGfreq':AGfreq, 'AGamp':AGamp,'ifTestSig':ifTestSig,
                    'ifAWG': ifAWG, 'ifIQ': ifIQ,'amp_MW_mix':amp_MW_mix, 'ifGreenLow':ifGreenLow, 'ifHiloExtra':ifHiloExtra,'ifIonInit':ifIonInit,
                    'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
                    'hilo_margin_start':hilo_margin_start, 'hilo_margin_end':hilo_margin_end,'hilo_min':hilo_min}
        ####### Random-phase noise ######
        if ifTestSig==1:
            AG = AG33522A()
            AG.disable_PM()
            AG.disable_RFOutput()
            time.sleep(1)
            
            if ifRndPhaseNoise==1:
                AG.set_PMsource()
                AG.set_PMfunction(function='NOIS')
                AG.set_PMdeviation()
                AG.set_noiseBandwidth(bandwidth=AGBW)

            AG.apply(function='SIN', freq=AGfreq, amplitude=AGamp, DCoffset=0)
            if ifRndPhaseNoise==1: AG.enable_PM()
        else:
            AG = AG33522A()
            AG.disable_PM()
            AG.disable_RFOutput()
            AG.disable_PM(channel=2)
            AG.disable_RFOutput(channel=2)
        #################################
        start = time.time()
        SCCRRPhotonStatChargeCheckShelveObject = SCCRRPhotonStatChargeCheckShelve(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
        SCCRRPhotonStatChargeCheckShelveObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        SCCRRPhotonStatChargeCheckShelveObject.close()