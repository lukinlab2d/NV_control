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

def ns2cycles(time, samp_rate=1e7):
        return int(time/1e9*samp_rate)
    

class RabiRRDualNV(Instrument):

    def __init__(self, name='RabiRRDualNVObject', settings=None, ifPlotPulse=True, **kwargs) -> None:

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
        self.MWI2Param =         {'delay_time': 2, 'channel':settings['MWI2_channel']}
        self.MWQ2Param =         {'delay_time': 2, 'channel':settings['MWQ2_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        self.AWG2Param =         {'delay_time': 2, 'channel':settings['AWG2_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
    
        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,
                          'AWG': self.AWGParam, 'AWG2': self.AWG2Param,
                          'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                          'LaserRead2': self.LaserRead2Param, 'LaserInit': self.LaserInitParam,
                          'MW_I2': self.MWI2Param, 'MW_Q2': self.MWQ2Param, 'MWswitch2': self.MWswitch2Param,
                          'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}; self.RRtrackingSettings = self.settings['RRtrackingSettings']; self.RRtrackingSettings2 = self.settings['RRtrackingSettings2']
        self.metadata.update(self.settings)

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        ############################################################################################
        # Microwave and velocity 1
        self.tausArray = self.settings['tausArray']
        self.SRSnum = self.settings['SRSnum']; MWPower = self.settings['MWPower']; MWFreq = self.settings['MWFreq']
        self.SDGnum = self.settings['SDGnum']; self.ifIQ = self.settings['ifIQ']
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
        ###########################################################################################

        self.add_parameter(
            name = "sig",
            parameter_class  = Signal,
        )
        self.add_parameter(
            name = "ref",
            parameter_class = Reference,
        )
        self.add_parameter(
            name = "sig2",
            parameter_class = Signal2,
            settings = self.settings,
            measurementObject = self
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
        if (self.SRSnum != 3) and (self.SRSnum != 4) and self.ifIQ:
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
        ifAWG2 = self.settings['ifAWG2']; self.ifAWG2 = ifAWG2
        if self.ifAWG2:
            self.AWG2 = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum2)
            global AWG2; AWG2 = self.AWG2

        # Velocity object 1
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
            sig2.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig2,sig,ref2,ref,
                                            qctask(sig2.plotPulseSequences),
                                            )

        data = loop.get_data_set(name='RabiRRDualNV')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.RabiRRDualNVObject_sig2, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig2'
            )
        plot.add(data.RabiRRDualNVObject_sig, name='sig')
        plot.add(data.RabiRRDualNVObject_ref2, name='ref2')
        plot.add(data.RabiRRDualNVObject_ref, name='ref')

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

        # Tracking NV
        self.hasTracked1 = 0; self.hasTracked2 = 0
        if self.settings['laserRead_channel'] == 5 or self.settings['laserRead_channel'] == 14:
                threshold_scanVpz = self.RRtrackingSettings['threshold_scanVpz']
                threshold_scanVpz2 = self.RRtrackingSettings2['threshold_scanVpz']
                
                datafile = self.getDataFilename()
                x_s, sig, ref, sig2, ref2 = dr.readDataNoPlotDual(datafile)
                ref = np.array(ref); ref2 = np.array(ref2)
                
                self.vpz = self.vel_vpz_target; self.vpz2 = self.vel_vpz_target2

                if self.RRtrackingSettings['if_tracking'] == 1 and np.max(ref) < threshold_scanVpz: 
                    print()
                    print('-----------------Start line tracking for NV1---------------------------')
                    ScanRRFreqObject = ScanRRFreq(settings=self.RRtrackingSettings, ifPlotPulse=0)
                    self.vpz = ScanRRFreqObject.runScanInPulseSequence()
                    vel.set_vpiezo(self.vpz)
                    vel.set_ready()
                    print("Set Vpiezo 1 to " + str(np.round(self.vpz,1)) + ' %') # set the Velocity's piezo voltage
                    print('-----------------End line tracking for NV1---------------------------')
                    print()
                    self.timeLastRRtracking = time.time()
                    self.hasTracked1 = 1
                if self.RRtrackingSettings2['if_tracking'] == 1 and np.max(ref2) < threshold_scanVpz2: 
                    print()
                    print('-----------------Start line tracking for NV2---------------------------')
                    ScanRRFreqObject = ScanRRFreq(settings=self.RRtrackingSettings2, ifPlotPulse=0)
                    self.vpz2 = ScanRRFreqObject.runScanInPulseSequence()
                    vel2.set_vpiezo(self.vpz2)
                    vel2.set_ready()
                    print("Set Vpiezo 2 to " + str(np.round(self.vpz2,1)) + ' %') # set the Velocity's piezo voltage
                    print('-----------------End line tracking for NV2---------------------------')
                    print()
                    self.timeLastRRtracking = time.time()
                    self.hasTracked2 = 1
    
        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        self.srs2.disable_RFOutput()
        self.srs2.disableModulation()
        if self.ifAWG: 
            self.AWG.turn_off()
        if self.ifAWG2: 
            self.AWG2.turn_off()
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/RabiRRDualNVObject_sig2_set.dat'

    
class Signal2(Parameter):
    def __init__(self, name='sig2', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.loopCounter = 0
        self.numOfRepumpVpz = 0
        self.numOfRepumpVpz2 = 0
        self.settings = settings
        self.tausArray = self.settings['tausArray']
        self.RabiRRDualNVObject = measurementObject
        self.RRtrackingSettings = self.settings['RRtrackingSettings']
        self.RRtrackingSettings2 = self.settings['RRtrackingSettings2']
        self.timeLastRRtracking = time.time()

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']
        self.ifAWG = self.settings['ifAWG']
        self.ifAWG2 = self.settings['ifAWG2']

        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.vel_vpz_target2 = self.settings['vel_vpz_target2']
        self.ifFakeRabi = self.settings['ifFakeRabi']
        self.timeNow = time.time()

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))
        self.ctrtask.stop(); self.ctrtask.close()

        rate = xLineData
        sig2 = rate[::self.num_reads_per_iter]
        sig = rate[1::self.num_reads_per_iter]
        ref2 = rate[2::self.num_reads_per_iter]
        ref = rate[3::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)/(self.read_duration2/1e9)/1e3
        global ref_avg;  ref_avg = np.average(ref)/(self.read_duration2/1e9)/1e3
        global sig_avg2;  sig_avg2 = np.average(sig2)/(self.read_duration/1e9)/1e3
        global ref_avg2;  ref_avg2 = np.average(ref2)/(self.read_duration/1e9)/1e3

        # Piezo repumping for RO NV 1
        if self.settings['laserRead_channel'] == 5 or self.settings['laserRead_channel'] == 14:
            if self.RRtrackingSettings['if_tracking'] == 1:
                threshold_repumpVpz = self.RRtrackingSettings['threshold_repumpVpz']
                if ref_avg < threshold_repumpVpz:
                    if np.mod(self.numOfRepumpVpz,5) == 0:
                        print()
                        print('-----------------Start resetting Vpiezo for NV1---------------------------')
                        vel.set_vpiezo(50)
                        vel.waitUntilComplete()
                        vel.set_ready()
                        time.sleep(0.7)
                        vel.set_vpiezo(2)
                        vel.waitUntilComplete()
                        vel.set_ready()
                        time.sleep(0.7)
                        for i in range(1):
                            vel.set_vpiezo(self.vel_vpz_target)
                            vel.waitUntilComplete()
                            vel.set_ready()
                            time.sleep(0.7)
                        print('-----------------End resetting Vpiezo for NV1---------------------------')
                        print()
                        self.timeLastRRtracking = time.time()
                    self.numOfRepumpVpz += 1
        
        # Piezo repumping for RO NV 2
        if self.settings['laserRead2_channel'] == 5 or self.settings['laserRead2_channel'] == 14:
            if self.RRtrackingSettings2['if_tracking'] == 1:
                threshold_repumpVpz2 = self.RRtrackingSettings2['threshold_repumpVpz']
                if ref_avg2 < threshold_repumpVpz2:
                    if np.mod(self.numOfRepumpVpz2,5) == 0:
                        print()
                        print('-----------------Start resetting Vpiezo for NV2---------------------------')
                        vel2.set_vpiezo(50)
                        vel2.waitUntilComplete()
                        vel2.set_ready()
                        time.sleep(0.7)
                        vel2.set_vpiezo(2)
                        vel2.waitUntilComplete()
                        vel2.set_ready()
                        time.sleep(0.7)
                        for i in range(1):
                            vel2.set_vpiezo(self.vel_vpz_target2)
                            vel2.waitUntilComplete()
                            vel2.set_ready()
                            time.sleep(0.7)
                        print('-----------------End resetting Vpiezo for NV2---------------------------')
                        print()
                        self.timeLastRRtracking = time.time()
                    self.numOfRepumpVpz2 += 1
                    
        return sig_avg2
    
    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        if self.ifFakeRabi:
            print("Time from last loop = " + str(time.time()-self.timeNow))
        else:
            print("Loop " + str(self.loopCounter))

        # Pulse parameters
        num_loops               = self.settings['num_loops']
        laser_init_delay        = self.settings['laser_init_delay'];    laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];  laser_to_DAQ_delay  = self.settings['laser_to_DAQ_delay']
        read_duration           = self.settings['read_duration'];       read_laser_duration = self.settings['read_laser_duration']
        read_duration2          = self.settings['read_duration2'];      read_laser_duration2= self.settings['read_laser_duration2']
        MW_to_read_delay        = self.settings['MW_to_read_delay']
        shift_btwn_2NV_MW       = self.settings['shift_btwn_2NV_MW'];   shift_btwn_2NV_read = self.settings['shift_btwn_2NV_read']
        laser_to_DAQ_delay2     = self.settings['laser_to_DAQ_delay2']
        AWG_buffer              = self.settings['AWG_buffer'];          AWG_output_delay    = self.settings['AWG_output_delay']  
        AWG_buffer2             = self.settings['AWG_buffer2'];         AWG_output_delay2   = self.settings['AWG_output_delay2']         
        
        if self.ifFakeRabi:
            MWI_duration = self.settings['MWI_duration']
            MWI_duration2 = self.settings['MWI_duration2']
        else:
            MWI_duration        = tau_ns
            MWI_duration2       = tau_ns

        when_init_end              = laser_init_delay + laser_init_duration
        MWI_delay                  = when_init_end    + laser_to_MWI_delay; self.MWI_delay = MWI_delay
        MW_delay_for_AWG           = MWI_delay - AWG_output_delay
        MW_duration_for_AWG        = int(2*int((2*AWG_buffer + MWI_duration + 1)/2)) # to make it even

        if self.ifAWG:
            when_pulse_end = MWI_delay + MW_duration_for_AWG
        else:
            when_pulse_end = MWI_delay + MWI_duration
     
        laser_read_signal_delay    = when_pulse_end   + MW_to_read_delay
        laser_read_signal_duration = read_laser_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration

        MWI_delay2                  = when_init_end   + laser_to_MWI_delay + shift_btwn_2NV_MW; self.MWI_delay2 = MWI_delay2
        MW_delay_for_AWG2           = MWI_delay2 - AWG_output_delay2
        MW_duration_for_AWG2        = int(2*int((2*AWG_buffer2 + MWI_duration2 + 1)/2)) # to make it even

        if self.ifAWG2:
            when_pulse_end2 = MWI_delay2 + MW_duration_for_AWG2
        else:
            when_pulse_end2 = MWI_delay2 + MWI_duration2
        
        laser_read_signal_delay2    = when_pulse_end2 + MW_to_read_delay + shift_btwn_2NV_read
        laser_read_signal_duration2 = read_laser_duration2
        when_laser_read_signal_end2 = laser_read_signal_delay2 + laser_read_signal_duration2
        read_signal_delay2          = laser_read_signal_delay2 + laser_to_DAQ_delay2
        read_signal_duration2       = read_duration2
        when_read_signal_end2       = read_signal_delay2 + read_signal_duration2
        
        laser_init_ref_delay        = np.max((when_read_signal_end2, when_read_signal_end)) + laser_init_delay
        when_init_ref_end           = laser_init_ref_delay + laser_init_duration
        laser_read_ref_delay        = when_init_ref_end + laser_to_MWI_delay + MWI_duration + MW_to_read_delay
        laser_read_ref_duration     = read_laser_duration
        when_laser_read_ref_end     = laser_read_ref_delay + laser_read_ref_duration
        read_ref_delay              = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration           = read_duration
        when_read_ref_end           = read_ref_delay + read_ref_duration

        laser_read_ref_delay2    = when_init_ref_end + laser_to_MWI_delay + shift_btwn_2NV_MW + MWI_duration + MW_to_read_delay + shift_btwn_2NV_read
        laser_read_ref_duration2 = read_laser_duration2
        when_laser_read_ref_end2 = laser_read_ref_delay2 + laser_read_ref_duration2
        read_ref_delay2          = laser_read_ref_delay2 + laser_to_DAQ_delay2
        read_ref_duration2       = read_duration2
        when_read_ref_end2       = read_ref_delay2 + read_ref_duration2

        self.read_duration = read_signal_duration
        self.read_duration2 = read_signal_duration2
        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")    
        
        if self.ifAWG:
            global ch1plot; global ch2plot
            ch1plot, ch2plot = AWG.send_fastRabi_seq(pulse_width=int(MWI_duration), buffer=int(AWG_buffer))
        if self.ifAWG2:
            global ch1plot2; global ch2plot2
            ch1plot2, ch2plot2 = AWG2.send_fastRabi_seq(pulse_width=int(MWI_duration2), buffer=int(AWG_buffer2))

        # Make pulse sequence (per each freq)
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay,    duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,    duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('Counter',      read_signal_delay,       duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',      read_ref_delay,          duration=int(read_ref_duration))] # times are in ns
        if self.ifAWG:
            pulse_sequence += [spc.Pulse('AWG',          MW_delay_for_AWG,    duration=50)]
        else:
            pulse_sequence += [spc.Pulse('MWswitch',     MWI_delay,           duration=int(MWI_duration))]


        pulse_sequence += [spc.Pulse('LaserRead2',    laser_read_signal_delay2, duration=int(laser_read_signal_duration2))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead2',    laser_read_ref_delay2,    duration=int(laser_read_ref_duration2))]
        pulse_sequence += [spc.Pulse('Counter',       read_signal_delay2,       duration=int(read_signal_duration2))] 
        pulse_sequence += [spc.Pulse('Counter',       read_ref_delay2,          duration=int(read_ref_duration2))] 
        if self.ifAWG2:
            pulse_sequence += [spc.Pulse('AWG2',          MW_delay_for_AWG2,    duration=50)]
        else:
            pulse_sequence += [spc.Pulse('MWswitch2',     MWI_delay2,           duration=int(MWI_duration2))]

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

        if not self.ifFakeRabi:
            print('Set tau to ' + str(tau_ns) + " ns")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1

        if self.ifFakeRabi:
            self.timeNow = time.time()
    
    def plotPulseSequences(self):
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False,
                                            readColor=self.readColor, 
                                            initColor=self.initColor)
                if (self.ifAWG) and (self.ifAWG2):
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, self.MWI_delay, label1='chnl1-AWG1', label2='chnl2-AWG1')
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot2,ch2plot2, self.MWI_delay2,
                                                           fig=fig, offset1=17.75, offset2=17.5, label1='chnl1-AWG2', label2='chnl2-AWG2')
                else:
                    if self.ifAWG:
                        fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, self.MWI_delay, label1='chnl1-AWG1', label2='chnl2-AWG1')
                    elif self.ifAWG2:
                        fig = plotPulseObject.makePulsePlotAWG(ch1plot2,ch2plot2, self.MWI_delay2,
                                                           offset1=17.75, offset2=17.5, label1='chnl1-AWG2', label2='chnl2-AWG2')
                    else:
                        fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.RabiRRDualNVObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
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

class Signal(Parameter):
    def __init__(self, name='sig',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_avg

class Reference2(Parameter):
    def __init__(self, name='ref2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg2