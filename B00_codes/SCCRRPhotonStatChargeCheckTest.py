"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.SCCRRPhotonStatChargeCheck import *
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A
import dropbox

####################################################################################################################
# Below are params for the first NV pair
reps = 1; ifPlotPulse=(reps==-1)
for iii in np.linspace(-1e8,-1e8,1):
    for whichNV in np.linspace(1,2,1):
        # SCCRRPhotonStatChargeCheck
        ifNeedVel=0; ifInitVpz=0; ifInitWvl=0; ifFancySpinInit=0; 
        ifAWG=1; ifMWDuringRead=1; ifMW2DuringRead=1; ifMWReadLowDutyCycle=0; ifIQ=ifAWG
        ifRandomized = 0;  sweepWhich = 'MWPower2'; ifHiloExtra=1; ifIonInit=0
        laserInit_channel = 3; laserIon_channel = 10; hiLoMWPwr_channel = 17; 
        ifTestSig=0; ifRndPhaseNoise=0; AGBW = 25e3; AGfreq = 2.5e6; AGamp = 0.2 # beware of heating!!
        if np.mod(whichNV,2)==1: 
            velNum = 1; vel_current = 62.7; vel_wvl = 637.25; vel_vpz_target=-1; laserRead_channel = 5
            SRSnum = 1; MWPower  = 0;  pi_time  = 72; MWFreq  = 3886e6   #NV D1 ms-1 2598.19e6
            SRSnum2= 3; MWPower2 = -4; pi_time2 = -1; MWFreq2 = 1930e6   #NV D1 ms+1 3161.27e6
            MWI_channel  = 1; MWQ_channel  = 0;  MWswitch_channel  = 2; ion_duration= 2e3 #8e2
            MWI2_channel = 0; MWQ2_channel = 0;  MWswitch2_channel = 15; amp_MW_mix = 1
            SDGnum       = 1; AWG_channel  = 18; srate             = 2.5e8
        else:
            velNum = 2; vel_current = 67.1; vel_wvl = 636.88; vel_vpz_target=-1; laserRead_channel = 14
            SRSnum = 2; MWPower  = -9.5; pi_time  = 44; MWFreq = 2788.70e6   #NV D2, ms-1 -9.5 44ns -8.8
            SRSnum2= 4; MWPower2 = 0.00; pi_time2 = -1; MWFreq2= 3037.20e6   #NV D2 ms+1 #power=10
            MWI_channel  = 12; MWQ_channel  = 13; MWswitch_channel  = 11; ion_duration= 4e3
            MWI2_channel = 0;  MWQ2_channel = 0;  MWswitch2_channel = 16; amp_MW_mix = 0.7
            SDGnum       = 2;  AWG_channel  = 19; srate             = 2.5e8
        
        # tausArray = np.array((1e2,5e2,7e2,1e3,1.5e3,2e3,3e3,4e3,6e3,8e3))#,15e3,17.5e3,20e3)) # sweep ti
        # tausArray = np.array((2e5,5e5,1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,10e6))#,12e6,15e6)) # sweep tr
        # tausArray = np.linspace(4,160,40)         # rabi
        # tausArray = np.array((5e5,7e5,1e6,1.25e6,1.5e6,2e6))#,3e6))#,4e6))#,5e6,6.5e6,8e6))#,10e6,12e6)) # tinit
        # tausArray = np.array((50,100,250,500,1e3,2e3,3e3,4e3)) # tionInit
        # tausArray = np.array((5e2,1e3,3e3, 5e3, 1e4, 3e4, 5e4,1e5,3e5,5e5,7e5,9e5,11e5,14e5,17e5,20e5)) #i2r
        tausArray = np.linspace(-8,-2,4) # amp_MW_mix / MWPower2

        num_loops                    = int(3e3);   nSpinInit                 = 0
        laser_init_delay             = 1e2;        laser_init_duration       = int(20e5)
        check_duration               = 2e4
        laserSwitch_delay_directory  = {3:860, 6:1160, 9:1160, 7:900, 5:1660, 10:172, 14:800}
        RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
        laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 160
        pi_to_ion_delay              = 1e2;        
        ion_to_read_delay            = 3e5;        DAQ_duration              = 4e6
        DAQ_to_laser_off_delay       = 1e2;        laserRead_to_MWmix        = 0
        iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
        AWG_buffer                   = 40;         AWG_output_delay          = 1448
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
        
        settings = {'tausArray': tausArray, 'dbx':dbx,
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
                    'ifRandomized':           ifRandomized,            'pi_to_ion_delay': pi_to_ion_delay,
                    'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
                    'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,
                    'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    'srate': int(srate),
                    'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
                    'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
                    'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
                    'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
                    'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'laser_ionInit_duration':laser_ionInit_duration,
                    'ifRndPhaseNoise':ifRndPhaseNoise, 'AGBW':AGBW, 'AGfreq':AGfreq, 'AGamp':AGamp,'ifTestSig':ifTestSig,
                    'ifAWG': ifAWG, 'ifIQ': ifIQ,'amp_MW_mix':amp_MW_mix, 'ifHiloExtra':ifHiloExtra,'ifIonInit':ifIonInit,
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
        SCCRRPhotonStatChargeCheckObject = SCCRRPhotonStatChargeCheck(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
        SCCRRPhotonStatChargeCheckObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        SCCRRPhotonStatChargeCheckObject.close()

        # if ifTestSig==1:
        #     AG.disable_PM()
        #     AG.disable_RFOutput()
        #     time.sleep(1)



###################################################################################################
# # Below are params for the second NV pair
# reps = 1; ifPlotPulse=(reps==-1)
# for whichNV in np.linspace(1,2,2):
#     # SCCRRPhotonStatChargeCheck
#     ifNeedVel=0; ifInitVpz=0; ifInitWvl=0; ifFancySpinInit=0; 
#     ifAWG=1; ifMWDuringRead=1; ifMW2DuringRead=1; ifMWReadLowDutyCycle=0; ifIQ=ifAWG
#     ifRandomized = 0;  sweepWhich = 'tinit'; ifHiloExtra=1; ifIonInit=0
#     laserInit_channel = 3; laserIon_channel = 10; hiLoMWPwr_channel = 17; 
#     ifTestSig=0; ifRndPhaseNoise=0; AGBW = 25e3; AGfreq = 2.5e6; AGamp = 0.2 # beware of heating!!

#     if np.mod(whichNV,2)==1: 
#         velNum = 1; vel_current = 62.7; vel_wvl = 637.22; vel_vpz_target=-1; laserRead_channel = 5
#         ion_duration              = 4e3
#     else:
#         velNum = 2; vel_current = 67.1; vel_wvl = 636.88; vel_vpz_target=-1; laserRead_channel = 14
#         ion_duration              = 1.5e3

#     SRSnum = 1; MWPower  = -6.2; pi_time  = 40; MWFreq  = 2953.76e6   #NV D1 ms-1
#     SRSnum2= 3; MWPower2 = 6.00; pi_time2 = -1; MWFreq2 = 2886.745e6   #NV D1 ms+1
#     MWI_channel  = 1; MWQ_channel  = 0;  MWswitch_channel  = 2
#     MWI2_channel = 0; MWQ2_channel = 0;  MWswitch2_channel = 15; amp_MW_mix = 0.7
#     SDGnum       = 1; AWG_channel  = 18; srate             = 2.5e8
        
#     # tausArray = np.array((5e2,1e3,2e3,3e3,4e3,5e3,6.5e3,8e3,10e3,12.5e3,15e3,17.5e3,20e3,25e3)) # ti
#     # tausArray = np.array((1e6,2e6,3e6,4e6,5e6,6e6,7e6))#,12e6,15e6)) # tr
#     # tausArray = np.linspace(4,164,41)         # rabi
#     tausArray = np.array((1.1e6,1.5e6,2e6))#,3e6))#,4e6))#,5e6,6.5e6,8e6))#,10e6,12e6)) # tinit
#     # tausArray = np.array((50,100,250,500,1e3,2e3,3e3,4e3)) # tionInit
#     # tausArray = np.array((5e2,1e3,3e3, 5e3, 1e4, 3e4, 5e4,1e5,3e5,5e5,7e5,9e5,11e5,14e5,17e5,20e5)) #i2r
#     # tausArray = np.linspace(1,0.7,7) # amp_MW_mix / MWPower2
#     num_loops                    = int(2e3);   nSpinInit                 = 0
#     laser_init_delay             = 1e2;        laser_init_duration       = int(1500e3)
#     check_duration               = 2e4
#     laserSwitch_delay_directory  = {3:860, 6:1160, 9:1160, 7:900, 5:1660, 10:172, 14:800}
#     RRLaserSwitch_delay          = laserSwitch_delay_directory.get(laserRead_channel, 0)       
#     laser_to_pi_delay            = laserSwitch_delay_directory.get(laserInit_channel, 0) + 160
#     pi_to_ion_delay              = 1e2;        
#     ion_to_read_delay            = 3e3;        DAQ_duration              = 3e6
#     DAQ_to_laser_off_delay       = 1e2;        laserRead_to_MWmix        = 0
#     iznLaserSwitch_delay         = laserSwitch_delay_directory.get(laserIon_channel, 0)
#     AWG_buffer                   = 40;         AWG_output_delay          = 1448
#     MWmix_duration_short         = 0;          delay_between_MWmix       = 0
#     spinInit_RR_duration         = 0;          spinInit_RR_to_pi_delay   = 0*(RRLaserSwitch_delay+100)
#     spinInit_pi_to_RR_delay      = 0;          laser_ionInit_duration    = 10
#     hilo_margin_start = 40; hilo_margin_end = 40; hilo_min = 60

#     # Create a Dropbox client
#     tkFile = 'C:/Users/lukin2dmaterials/data/tk.txt'; lines = []
#     with open(tkFile, 'r') as file:
#         for line in file:
#             if "\n" in line: line = line[0:-1]
#             lines.append(line)
#     access_token = lines[0]
#     app_key = lines[1]
#     app_secret = lines[2]
#     refresh_token = lines[3]
#     #######################################
#     dbx = dropbox.Dropbox(oauth2_access_token=access_token,
#                           app_key=app_key, 
#                           app_secret=app_secret, 
#                           oauth2_refresh_token=refresh_token)
    
#     settings = {'tausArray': tausArray, 'dbx':dbx,
#                 'num_loops':num_loops, 'MWPower':MWPower, 'MWFreq': MWFreq, 'SRSnum': SRSnum, 'ifMWDuringRead':ifMWDuringRead,
#                 'MWI_channel': MWI_channel, 'MWQ_channel': MWQ_channel, 'MWswitch_channel': MWswitch_channel,
#                 'MWPower2':MWPower2, 'MWFreq2': MWFreq2, 'SRSnum2': SRSnum2, 'ifMW2DuringRead':ifMW2DuringRead,
#                 'MWI2_channel': MWI2_channel, 'MWQ2_channel': MWQ2_channel, 'MWswitch2_channel': MWswitch2_channel,
#                 'velNum': velNum, 'vel_current': vel_current, 'vel_wvl': vel_wvl, 'vel_vpz_target': vel_vpz_target,
#                 'ifInitVpz': ifInitVpz, 'ifInitWvl': ifInitWvl, 'ifNeedVel': ifNeedVel,
#                 'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
#                 'laser_to_pi_delay':      laser_to_pi_delay,       'DAQ_duration': DAQ_duration,
#                 'check_duration':         check_duration,
#                 'RRLaserSwitch_delay':    RRLaserSwitch_delay,     'pi_time': pi_time, 'pi_time2':pi_time2,
#                 'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'ion_to_read_delay': ion_to_read_delay,
#                 'laserRead_to_MWmix':     laserRead_to_MWmix,      'iznLaserSwitch_delay':iznLaserSwitch_delay,
#                 'ifRandomized':           ifRandomized,            'pi_to_ion_delay': pi_to_ion_delay,
#                 'laserRead_channel':      laserRead_channel,       'laserInit_channel':   laserInit_channel,
#                 'laserIon_channel':       laserIon_channel,        'ion_duration': ion_duration,
#                 'ifMWReadLowDutyCycle':   ifMWReadLowDutyCycle,    'srate': int(srate),
#                 'MWmix_duration_short':   MWmix_duration_short,    'delay_between_MWmix': delay_between_MWmix,
#                 'ifFancySpinInit':        ifFancySpinInit,         'nSpinInit': nSpinInit,
#                 'spinInit_RR_duration':   spinInit_RR_duration,    'spinInit_RR_to_pi_delay': spinInit_RR_to_pi_delay,
#                 'spinInit_pi_to_RR_delay':spinInit_pi_to_RR_delay, 'sweepWhich':sweepWhich,
#                 'hiLoMWPwr_channel':      hiLoMWPwr_channel,       'laser_ionInit_duration':laser_ionInit_duration,
#                 'ifRndPhaseNoise':ifRndPhaseNoise, 'AGBW':AGBW, 'AGfreq':AGfreq, 'AGamp':AGamp,'ifTestSig':ifTestSig,
#                 'ifAWG': ifAWG, 'ifIQ': ifIQ,'amp_MW_mix':amp_MW_mix, 'ifHiloExtra':ifHiloExtra,'ifIonInit':ifIonInit,
#                 'SDGnum': SDGnum,   'AWG_channel':AWG_channel,   'AWG_buffer':AWG_buffer,   'AWG_output_delay':AWG_output_delay,
#                 'hilo_margin_start':hilo_margin_start, 'hilo_margin_end':hilo_margin_end,'hilo_min':hilo_min}
#     ####### Random-phase noise ######
#     if ifTestSig==1:
#         AG = AG33522A()
#         AG.disable_PM()
#         AG.disable_RFOutput()
#         time.sleep(1)
        
#         if ifRndPhaseNoise==1:
#             AG.set_PMsource()
#             AG.set_PMfunction(function='NOIS')
#             AG.set_PMdeviation()
#             AG.set_noiseBandwidth(bandwidth=AGBW)

#         AG.apply(function='SIN', freq=AGfreq, amplitude=AGamp, DCoffset=0)
#         if ifRndPhaseNoise==1: AG.enable_PM()
#     else:
#         AG = AG33522A()
#         AG.disable_PM()
#         AG.disable_RFOutput()
#         AG.disable_PM(channel=2)
#         AG.disable_RFOutput(channel=2)
#     #################################
#     start = time.time()
#     SCCRRPhotonStatChargeCheckObject = SCCRRPhotonStatChargeCheck(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
#     SCCRRPhotonStatChargeCheckObject.runScan()
#     print('Total time = ' + str(time.time() - start) + ' s')
#     SCCRRPhotonStatChargeCheckObject.close()

#     # if ifTestSig==1:
#     #     AG.disable_PM()
#     #     AG.disable_RFOutput()
#     #     time.sleep(1)