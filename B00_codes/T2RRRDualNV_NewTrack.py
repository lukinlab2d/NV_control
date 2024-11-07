"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
from  qcodes.actions import Task as qctask
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot
import numpy as np
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS
from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity
from qcodes_contrib_drivers.drivers.Siglent.SDG6022X import SDG6022X
from copy import deepcopy

import nidaqmx, time
from nidaqmx.constants import *
from qcodes.instrument.base import Instrument
from qcodes.instrument.parameter import Parameter
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
from PIL import Image
from B00_codes.PlotPulse import *  
from B00_codes.Confocal import *
import B00_codes.dataReader as dr
from B00_codes.ScanRRFreq import *  

class T2RRRDualNV_NewTrack(Instrument):

    def __init__(self, name='T2RRRDualNV_NewTrackObject', settings=None, ifPlotPulse=True, **kwargs) -> None:

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.LaserRead2Param =   {'delay_time': 2, 'channel':settings['laserRead2_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        self.MW2Param =         {'delay_time': 2, 'channel':settings['MWI2_channel']}
        self.MWQ2Param =         {'delay_time': 2, 'channel':settings['MWQ2_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        self.AWG2Param =         {'delay_time': 2, 'channel':settings['AWG2_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
    
        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,
                          'AWG': self.AWGParam, 'AWG2': self.AWG2Param,
                          'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                          'LaserRead2': self.LaserRead2Param, 'LaserInit': self.LaserInitParam,
                          'MW_I2': self.MW2Param, 'MW_Q2': self.MWQ2Param, 'MWswitch2': self.MWswitch2Param,
                          'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        ############################################################################################
        # Microwave and velocity 1
        self.tausArray = self.settings['tausArray']
        self.SRSnum = self.settings['SRSnum']; MWPower = self.settings['MWPower']; MWFreq = self.settings['MWFreq']
        self.SDGnum = self.settings['SDGnum']
        self.ifNeedVel1 = self.settings['ifNeedVel1']; self.velNum = self.settings['velNum']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

        # Microwave and velocity 2
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        self.SDGnum2 = self.settings['SDGnum2']
        self.ifNeedVel2 = self.settings['ifNeedVel2']; self.velNum2 = self.settings['velNum2']
        vel_current2 = self.settings['vel_current2']; vel_wvl2 = self.settings['vel_wvl2']; 
        self.vel_vpz_target2 = self.settings['vel_vpz_target2']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)
        ###########################################################################################

        self.add_parameter(
            name = "sig",
            parameter_class  = Signal,
            settings = self.settings,
            measurementObject = self
        )
        self.add_parameter(
            name = "ref",
            parameter_class = Reference,
        )
        self.add_parameter(
            name = "sig2",
            parameter_class = Signal2,
        )
        self.add_parameter(
            name = "ref2",
            parameter_class = Reference2,
        )
        
        self.savedPulseSequencePlots = {}

        # SRS object 1
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(MWFreq) #Hz
        self.srs.set_RFAmplitude(MWPower) #dBm
        if (self.SRSnum != 3) and (self.SRSnum != 4):
            self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        # SRS object 2
        self.srs2 = SRS(SRSnum=self.SRSnum2)
        self.srs2.set_freq(MWFreq2) #Hz
        self.srs2.set_RFAmplitude(MWPower2) #dBm
        if (self.SRSnum2 != 3) and (self.SRSnum2 != 4):
            self.srs2.enableIQmodulation()
        self.srs2.enable_RFOutput()

        # AWG object 1
        ifAWG = self.settings['ifAWG']; self.ifAWG = ifAWG
        if self.ifAWG:
            self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum)
            global AWG; AWG = self.AWG

        # AWG object 2
        if self.ifAWG:
            self.AWG2 = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum2)
            global AWG2; AWG2 = self.AWG2

        # Velocity object
        if self.ifNeedVel1:
            self.vel = Velocity(velNum=self.velNum, 
                                ifInitVpz=self.ifInitVpz, ifInitWvl=self.ifInitWvl,
                                initWvl=vel_wvl)
            if self.ifInitWvl: 
                self.vel.set_track()
                time.sleep(0.5)
                self.vel.set_wvl(vel_wvl)
                time.sleep(1)
                self.vel.set_ready()
                self.vel.set_vpiezo(2)
                self.vel.waitUntilComplete()
                self.vel.set_ready()
                time.sleep(0.7)
            if self.ifInitVpz:
                self.vel.set_vpiezo(2)
                self.vel.waitUntilComplete()
                self.vel.set_ready()
                time.sleep(0.7)

            self.vel.set_current(vel_current)
            
            for i in range(1):
                self.vel.set_vpiezo(self.vel_vpz_target)
                self.vel.waitUntilComplete()
                self.vel.set_ready()
                time.sleep(0.7)
            global vel; vel = self.vel; 

        # Velocity object 2
        if self.ifNeedVel2:
            self.vel2 = Velocity(velNum=self.velNum2, 
                                ifInitVpz=self.ifInitVpz, ifInitWvl=self.ifInitWvl,
                                initWvl=vel_wvl2)
            if self.ifInitWvl: 
                self.vel2.set_track()
                time.sleep(0.5)
                self.vel2.set_wvl(vel_wvl2)
                time.sleep(1)
                self.vel2.set_ready()
                self.vel2.set_vpiezo(2)
                self.vel2.waitUntilComplete()
                self.vel2.set_ready()
                time.sleep(0.7)
            if self.ifInitVpz:
                self.vel2.set_vpiezo(2)
                self.vel2.waitUntilComplete()
                self.vel2.set_ready()
                time.sleep(0.7)

            self.vel2.set_current(vel_current2)
            
            for i in range(1):
                self.vel2.set_vpiezo(self.vel_vpz_target2)
                self.vel2.waitUntilComplete()
                self.vel2.set_ready()
                time.sleep(0.7)
            global vel2; vel2 = self.vel2

        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        global srs; srs = self.srs; global srs2; srs2 = self.srs2
        
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref
        sig2 = self.sig2
        ref2 = self.ref2

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,sig2,ref,ref2,
                                            qctask(sig.plotPulseSequences),
                                            )

        data = loop.get_data_set(name='T2RRRDualNV_NewTrack')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.T2RRRDualNV_NewTrackObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.T2RRRDualNV_NewTrackObject_sig2, name='sig2')
        plot.add(data.T2RRRDualNV_NewTrackObject_ref, name='ref')
        plot.add(data.T2RRRDualNV_NewTrackObject_ref2, name='ref2')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        img = Image.open(dataPlotFile)
        img.show()
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)
        
        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        self.srs2.disable_RFOutput()
        self.srs2.disableModulation()
        if self.ifAWG: 
            self.AWG.turn_off()
            self.AWG2.turn_off()
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T2RRRDualNV_NewTrackObject_sig_set.dat'

    
class Signal(Parameter):
    def __init__(self, name='sig', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.loopCounter = 0
        self.settings = settings
        self.tausArray = self.settings['tausArray']
        self.T2RRRDualNV_NewTrackObject = measurementObject
        self.RRtrackingSettings = self.settings['RRtrackingSettings']
        self.RRtrackingSettings2 = self.settings['RRtrackingSettings2']

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']
        self.ifAWG = self.settings['ifAWG']

        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.vel_vpz_target2 = self.settings['vel_vpz_target2']

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))

        self.ctrtask.stop(); self.ctrtask.close()

        rate = xLineData
        sig = rate[::self.num_reads_per_iter]
        sig2 = rate[1::self.num_reads_per_iter]
        ref = rate[2::self.num_reads_per_iter]
        ref2 = rate[3::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)/(self.read_duration/1e9)/1e3
        global ref_avg;  ref_avg = np.average(ref)/(self.read_duration/1e9)/1e3
        global sig_avg2;  sig_avg2 = np.average(sig2)/(self.read_duration2/1e9)/1e3
        global ref_avg2;  ref_avg2 = np.average(ref2)/(self.read_duration2/1e9)/1e3

        # Line tracking
        if self.settings['laserRead_channel'] == 5 or self.settings['laserRead_channel'] == 14:
            threshold_scanVpz = self.RRtrackingSettings['threshold_scanVpz']
            threshold_scanVpz2 = self.RRtrackingSettings2['threshold_scanVpz']
            #######################################################
            # NV 1
            if self.RRtrackingSettings['if_tracking'] == 1 and np.max((np.max(sig_avg), np.max(ref_avg))) < threshold_scanVpz:
                print()
                time_sleep_after_scan = self.RRtrackingSettings['time_sleep_after_scan']
                wvl_correction = self.RRtrackingSettings['wvl_correction']
                scan_lower_margin = self.RRtrackingSettings['scan_lower_margin']
                scan_upper_margin = self.RRtrackingSettings['scan_upper_margin']
                num_point_scan = self.RRtrackingSettings['num_point_scan']
                print('-----------------Start line tracking for NV1---------------------------')
                lockStatusFile = 'C:\\Users\\lukin2dmaterials\\miniconda3\\envs\\NV_control\\B00_codes\\laser1_lockStatus.txt'
                
                # Turn off laser lock
                try:
                    with open(lockStatusFile, 'r+') as file:
                        file.seek(0)
                        file.write(str(0))
                except OSError as e:
                    print(f"Error: {e}")
                    try:
                        with open(lockStatusFile, 'r+') as file:
                            file.seek(0)
                            file.write(str(0))
                    except OSError as e:
                        print(f"Error: {e}")
                
                # Scan vpz
                lockStatusParam = np.loadtxt(lockStatusFile)
                vpz_old = lockStatusParam[2]
                start = vpz_old - scan_lower_margin; stop = vpz_old + scan_upper_margin; num_sweep_points = num_point_scan
                vpzArray = np.linspace(start, stop, num_sweep_points)
                self.RRtrackingSettings['vpzArray'] = vpzArray
                
                ScanRRFreqObject = ScanRRFreq(settings=self.RRtrackingSettings, ifPlotPulse=0)
                vpz_new, wvl = ScanRRFreqObject.runScanInPulseSequenceMonitoringWM()

                # Turn on laser lock again
                with open(lockStatusFile, 'r+') as file:
                    wvl = wvl + wvl_correction
                    file.seek(0)
                    file.write(str(1)+'\n')
                    file.write(str(wvl)+'\n')
                    file.write(str(vpz_new+0.01))
                    file.truncate()
                
                print("Lock laser 1 to " + str(wvl) + ' THz') # set the Velocity's piezo voltage
                print('-----------------End line tracking for NV1, wait for 50 s---------------------------')
                print()
                time.sleep(time_sleep_after_scan)
            ########################################################################################################################
            # NV 2
            if self.RRtrackingSettings2['if_tracking'] == 1 and np.max((np.max(sig_avg2), np.max(ref_avg2))) < threshold_scanVpz2:
                print()
                time_sleep_after_scan = self.RRtrackingSettings2['time_sleep_after_scan']
                wvl_correction = self.RRtrackingSettings2['wvl_correction']
                scan_lower_margin = self.RRtrackingSettings2['scan_lower_margin']
                scan_upper_margin = self.RRtrackingSettings2['scan_upper_margin']
                num_point_scan = self.RRtrackingSettings2['num_point_scan']
                print('-----------------Start line tracking for NV2---------------------------')
                lockStatusFile = 'C:\\Users\\lukin2dmaterials\\miniconda3\\envs\\NV_control\\B00_codes\\laser2_lockStatus.txt'
                
                # Turn off laser lock
                try:
                    with open(lockStatusFile, 'r+') as file:
                        file.seek(0)
                        file.write(str(0))
                except OSError as e:
                    print(f"Error: {e}")
                    try:
                        with open(lockStatusFile, 'r+') as file:
                            file.seek(0)
                            file.write(str(0))
                    except OSError as e:
                        print(f"Error: {e}")
                
                # Scan vpz
                lockStatusParam = np.loadtxt(lockStatusFile)
                vpz_old = lockStatusParam[2]
                start = vpz_old - scan_lower_margin; stop = vpz_old + scan_upper_margin; num_sweep_points = num_point_scan
                vpzArray = np.linspace(start, stop, num_sweep_points)
                self.RRtrackingSettings2['vpzArray'] = vpzArray
                
                ScanRRFreqObject = ScanRRFreq(settings=self.RRtrackingSettings2, ifPlotPulse=0)
                vpz_new, wvl = ScanRRFreqObject.runScanInPulseSequenceMonitoringWM()

                # Turn on laser lock again
                with open(lockStatusFile, 'r+') as file:
                    wvl = wvl + wvl_correction
                    file.seek(0)
                    file.write(str(1)+'\n')
                    file.write(str(wvl)+'\n')
                    file.write(str(vpz_new + 0.01))
                    file.truncate()
                
                print("Lock laser 2 to " + str(wvl) + ' THz') # set the Velocity's piezo voltage
                print('-----------------End line tracking for NV2, wait for 50 s---------------------------')
                print()
                time.sleep(time_sleep_after_scan)
           
        return sig_avg
    
    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))

        NO_MS_EQUALS_1 = 0
        Q_FINAL = 1
        THREE_PI_HALF_FINAL = 2

        # Pulse parameters
        num_loops               = self.settings['num_loops']
        laser_init_delay        = self.settings['laser_init_delay'];    laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];  laser_to_DAQ_delay  = self.settings['laser_to_DAQ_delay'];  
        MW_duration=pi2time     = self.settings['pi_half'];             MW_duration2=pi2time2       = self.settings['pi_half2']
        read_duration           = self.settings['read_duration'];       read_laser_duration = self.settings['read_laser_duration']
        read_duration2          = self.settings['read_duration2'];      read_laser_duration2= self.settings['read_laser_duration2']
        MW_to_read_delay        = self.settings['MW_to_read_delay']
        shift_btwn_2NV_MW       = self.settings['shift_btwn_2NV_MW'];   shift_btwn_2NV_read = self.settings['shift_btwn_2NV_read']
        laser_to_DAQ_delay2     = self.settings['laser_to_DAQ_delay2']; normalized_style    = self.settings['normalized_style']
        AWG_buffer              = self.settings['AWG_buffer'];          AWG_output_delay    = self.settings['AWG_output_delay']  
        AWG_buffer2             = self.settings['AWG_buffer2'];         AWG_output_delay2   = self.settings['AWG_output_delay2']         
        
        MW_duration_for_AWG = int(2*int((AWG_buffer + 4*MW_duration + tau_ns + 1)/2))
        MW_duration_for_AWG2 = int(2*int((AWG_buffer2 + 4*MW_duration2 + tau_ns + 1)/2))
        
        if tau_ns/2 > 30:
            MWI_to_switch_delay  = self.settings['MWI_to_switch_delay']
        else: 
            MWI_to_switch_delay  = 0

        when_init_end              = laser_init_delay + laser_init_duration

        ################# Signal  ################################################################
        ###### NV1 ######################################
        MW_delay                  = when_init_end + laser_to_MWI_delay              
        global MW_del;      MW_del = MW_delay + AWG_output_delay

        MW2_delay      = MW_delay + pi2time + tau_ns;         MW2_duration = pi2time

        if self.ifAWG:
            when_pulse_end = MW_del + MW_duration_for_AWG
        else:
            when_pulse_end = MW2_delay + MW2_duration

        laser_read_signal_delay    = when_pulse_end   + MW_to_read_delay
        laser_read_signal_duration = read_laser_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration

        ###### NV2 ######################################
        MW_delay2                  = MW_delay + shift_btwn_2NV_MW 
        global MW_del2;     MW_del2 = MW_delay2 + AWG_output_delay2

        MW2_delay2        = MW_delay2 + pi2time2 + tau_ns;         MW2_duration2 = pi2time2

        if self.ifAWG:
            when_pulse_end2 = MW_del2 + MW_duration_for_AWG2
        else:
            when_pulse_end2 = MW2_delay2 + MW2_duration2   
        
        laser_read_signal_delay2    = when_pulse_end2 + MW_to_read_delay + shift_btwn_2NV_read
        laser_read_signal_duration2 = read_laser_duration2
        when_laser_read_signal_end2 = laser_read_signal_delay2 + laser_read_signal_duration2
        read_signal_delay2          = laser_read_signal_delay2 + laser_to_DAQ_delay2
        read_signal_duration2       = read_duration2
        when_read_signal_end2       = read_signal_delay2 + read_signal_duration2

        ################# Reference  ################################################################
        ###### NV1 ######################################
        laser_init_ref_delay        = np.max((when_read_signal_end2,when_read_signal_end)) + laser_init_delay
        when_init_ref_end           = laser_init_ref_delay + laser_init_duration

        MW3_delay         = when_init_ref_end + laser_to_MWI_delay;    MW3_duration = pi2time
        MW4_delay         = MW3_delay + MW3_duration + tau_ns;        MW4_duration = pi2time

        laser_read_ref_delay        = np.max((when_read_signal_end2,when_read_signal_end)) + laser_read_signal_delay
        laser_read_ref_duration     = read_laser_duration
        when_laser_read_ref_end     = laser_read_ref_delay + laser_read_ref_duration
        read_ref_delay              = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration           = read_duration
        when_read_ref_end           = read_ref_delay + read_ref_duration

        sig_to_ref_wait             = laser_read_ref_delay - 2*MW_duration_for_AWG - MW_del

        ###### NV2 ######################################
        MW3_delay2        = when_init_ref_end + laser_to_MWI_delay + shift_btwn_2NV_MW;    MW3_duration2 = pi2time
        MW4_delay2        = MW3_delay2 + MW3_duration2 + tau_ns;                           MW4_duration2 = pi2time

        laser_read_ref_delay2    = np.max((when_read_signal_end2,when_read_signal_end)) + laser_read_signal_delay2
        laser_read_ref_duration2 = read_laser_duration2
        when_laser_read_ref_end2 = laser_read_ref_delay2 + laser_read_ref_duration2
        read_ref_delay2          = laser_read_ref_delay2 + laser_to_DAQ_delay2
        read_ref_duration2       = read_duration2
        when_read_ref_end2       = read_ref_delay2 + read_ref_duration2

        sig_to_ref_wait2         = sig_to_ref_wait #laser_read_ref_delay2 - 2*MW_duration_for_AWG2 - MW_del2

        self.read_duration  = read_signal_duration
        self.read_duration2 = read_signal_duration2
        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")  
        

        if self.ifAWG:
            global ch1plot; global ch2plot
            ch1plot, ch2plot = AWG.send_T2R_seq(pi_2time=int(pi2time), tau = int(tau_ns), 
                                                buffer=int(AWG_buffer), sig_to_ref_wait=int(sig_to_ref_wait))
            global ch1plot2; global ch2plot2
            ch1plot2, ch2plot2 = AWG2.send_T2R_seq(pi_2time=int(pi2time2), tau = int(tau_ns), 
                                                buffer=int(AWG_buffer2), sig_to_ref_wait=int(sig_to_ref_wait2))
          
        ################################################################################################
        # Make pulse sequence (per each freq)
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_delay,             duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_ref_delay,         duration=int(laser_init_duration))] # times are in ns
        
        ###### NV1 ######################################
        pulse_sequence += [spc.Pulse('LaserRead', laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead', laser_read_ref_delay,          duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('Counter',   read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',   read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
        
        if self.ifAWG:
            pulse_sequence += [spc.Pulse('AWG',      MW_delay,      duration=50)]
        else:
            pulse_sequence += [spc.Pulse('MWswitch', MW_delay,    duration=int(MW_duration))]
            pulse_sequence += [spc.Pulse('MWswitch', MW2_delay,   duration=int(MW2_duration))]
        
            pulse_sequence += [spc.Pulse('MWswitch', MW3_delay,   duration=int(MW3_duration))]
            pulse_sequence += [spc.Pulse('MWswitch', MW4_delay,   duration=int(MW4_duration))]

            pulse_sequence += [spc.Pulse('MW_I', MW4_delay-MWI_to_switch_delay, duration=int(MW4_duration + 2*MWI_to_switch_delay))]
            pulse_sequence += [spc.Pulse('MW_Q', MW4_delay-MWI_to_switch_delay, duration=int(MW4_duration + 2*MWI_to_switch_delay))]
        
        ###### NV2 ######################################
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_signal_delay2, duration=int(laser_read_signal_duration2))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_ref_delay2,    duration=int(laser_read_ref_duration2))]
        pulse_sequence += [spc.Pulse('Counter',    read_signal_delay2,       duration=int(read_signal_duration2))] 
        pulse_sequence += [spc.Pulse('Counter',    read_ref_delay2,          duration=int(read_ref_duration2))] 
        if self.ifAWG:
            pulse_sequence += [spc.Pulse('AWG2',      MW_delay2,  duration=50)]
        else:
            pulse_sequence += [spc.Pulse('MWswitch2', MW_delay2,    duration=int(MW_duration2))]
            pulse_sequence += [spc.Pulse('MWswitch2', MW2_delay2,   duration=int(MW2_duration2))]
        
            pulse_sequence += [spc.Pulse('MWswitch2', MW3_delay2,   duration=int(MW3_duration2))]
            pulse_sequence += [spc.Pulse('MWswitch2', MW4_delay2,   duration=int(MW4_duration2))]

            pulse_sequence += [spc.Pulse('MW_I2', MW4_delay2-MWI_to_switch_delay, duration=int(MW4_duration2 + 2*MWI_to_switch_delay))]
            pulse_sequence += [spc.Pulse('MW_Q2', MW4_delay2-MWI_to_switch_delay, duration=int(MW4_duration2 + 2*MWI_to_switch_delay))]
        
        
        ################################################################################################
        self.pulse_sequence = pulse_sequence
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        self.pb.program_pb(pulse_sequence, num_loops=num_loops)

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        self.num_reads = int(num_loops * num_reads_per_iter)
        self.num_reads_per_iter = num_reads_per_iter

        # Pulse width counter. Timebase = signal; gate = PB signal
        self.ctrtask = nidaqmx.Task()
        pulseWidthChan = self.ctrtask.ci_channels.add_ci_pulse_width_chan( # define the pulse width counter
            counter = "cDAQ1Mod1/ctr0",
            name_to_assign_to_channel = "",
            min_val = 0,
            max_val = int(1e8),
            units = TimeUnits.TICKS,
            starting_edge = Edge.RISING,
            )
        self.ctrtask.timing.cfg_implicit_timing(
            sample_mode = AcquisitionType.CONTINUOUS,
            samps_per_chan = 2*self.num_reads # x2 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        print('Set tau to ' + str(tau_ns) + " ns")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1
    
    def plotPulseSequences(self):
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False,
                                            readColor=self.readColor, 
                                            initColor=self.initColor)
                if self.ifAWG:
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, MW_del, label1='ch1-AWG1', label2='ch2-AWG1')
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot2,ch2plot2,MW_del2,
                                                           fig=fig, offset1=17.75, offset2=17.5, label1='ch1-AWG2', label2='ch2-AWG2')
                else:
                    fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T2RRRDualNV_NewTrackObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1

    def close_turnOnAtEnd(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        # pb.turn_on_infinite(channels=channels)
    
class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg

class Signal2(Parameter):
    def __init__(self, name='sig2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_avg2

class Reference2(Parameter):
    def __init__(self, name='ref2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg2