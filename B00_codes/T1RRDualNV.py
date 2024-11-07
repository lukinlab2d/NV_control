"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
from copy import deepcopy
from qcodes.actions import Task as qctask
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot
import numpy as np
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS
from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity
from qcodes_contrib_drivers.drivers.Siglent.SDG6022X import SDG6022X

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

class T1RRDualNV(Instrument):

    def __init__(self, name='T1RRDualNVObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz, of the Pulse Blaster
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.LaserRead2Param =  {'delay_time': 2, 'channel':settings['laserRead2_channel']}
        self.LaserIonParam =    {'delay_time': 2, 'channel':settings['laserIon_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        self.MWI2Param =         {'delay_time': 2, 'channel':settings['MWI2_channel']}
        self.MWQ2Param =         {'delay_time': 2, 'channel':settings['MWQ2_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        self.hiLoMWPwrParam =    {'delay_time': 2, 'channel':settings['hiLoMWPwr_channel']}
        self.MWswitch3Param =    {'delay_time': 2, 'channel':settings['MWswitch3_channel']}
        self.MWswitch4Param =    {'delay_time': 2, 'channel':settings['MWswitch4_channel']}
        self.AWG2Param =         {'delay_time': 2, 'channel':settings['AWG2_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 
                          'LaserRead': self.LaserReadParam, 'LaserRead2': self.LaserRead2Param, 'LaserInit': self.LaserInitParam, 'LaserIon': self.LaserIonParam,
                        'AWG': self.AWGParam, 'AWG2': self.AWG2Param,
                        'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'MW_I2': self.MWI2Param, 'MW_Q2': self.MWQ2Param, 'MWswitch2': self.MWswitch2Param, 'hiLoMWPwr': self.hiLoMWPwrParam,
                        'MWswitch3': self.MWswitch3Param, 'MWswitch4': self.MWswitch4Param, 
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.tausArray = self.settings['tausArray']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)

        self.SRSnum = self.settings['SRSnum'];   MWPower = self.settings['MWPower'];   MWFreq = self.settings['MWFreq']
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        self.SRSnum3 = self.settings['SRSnum3']; MWPower3 = self.settings['MWPower3']; MWFreq3 = self.settings['MWFreq3']
        self.SRSnum4 = self.settings['SRSnum4']; MWPower4 = self.settings['MWPower4']; MWFreq4 = self.settings['MWFreq4']
        
        self.ifNeedVel = self.settings['ifNeedVel']; self.velNum = self.settings['velNum']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

        self.SDGnum = self.settings['SDGnum']; self.ifIQ = self.settings['ifIQ']
        self.SDGnum2 = self.settings['SDGnum2']; self.ifIQ2 = self.settings['ifIQ2']

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

        # SRS object 
        if True:
            self.srs = SRS(SRSnum=self.SRSnum)
            self.srs.set_freq(MWFreq) #Hz
            self.srs.set_RFAmplitude(MWPower) #dBm
            if (self.SRSnum != 3) and (self.SRSnum != 4) and self.ifIQ:
                self.srs.enableIQmodulation()
            self.srs.enable_RFOutput()
            global srs; srs = self.srs

            self.srs2 = SRS(SRSnum=self.SRSnum2)
            self.srs2.set_freq(MWFreq2) #Hz
            self.srs2.set_RFAmplitude(MWPower2) #dBm
            if (self.SRSnum2 != 3) and (self.SRSnum2 != 4) and self.ifIQ2:
                self.srs2.enableIQmodulation()
            self.srs2.enable_RFOutput()
            global srs2; srs2 = self.srs2
            
            self.srs3 = SRS(SRSnum=self.SRSnum3)
            self.srs3.set_freq(MWFreq3) #Hz
            self.srs3.set_RFAmplitude(MWPower3) #dBm
            if (self.SRSnum3 != 3) and (self.SRSnum3 != 4):
                self.srs3.enableIQmodulation()
            self.srs3.disable_RFOutput()
            global srs3; srs3 = self.srs3

            self.srs4 = SRS(SRSnum=self.SRSnum4)
            self.srs4.set_freq(MWFreq4) #Hz
            self.srs4.set_RFAmplitude(MWPower4) #dBm
            if (self.SRSnum4 != 3) and (self.SRSnum4 != 4):
                self.srs4.enableIQmodulation()
            self.srs4.disable_RFOutput()
            global srs4; srs4 = self.srs4

        # Velocity object
        if self.ifNeedVel:
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

            for i in range(5):
                self.vel.set_vpiezo(self.vel_vpz_target)
                self.vel.waitUntilComplete()
                self.vel.set_ready()
                time.sleep(0.7)
            global vel; vel = self.vel
        
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

        # Make Pulse Blaster, Counter, SRS global objects
        global pb
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref 

        sig2 = self.sig2 
        ref2 = self.ref2 

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,
                                            sig2,ref2,
                                            qctask(sig.plotPulseSequences),
                                            )
        data = loop.get_data_set(name='T1RRDualNV')
        data.add_metadata(self.settings)
        self.data = data
        
        loop.run()
        print('Data saved to ' + str(data.location) + '/')
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)
        
        self.srs.disable_RFOutput();  self.srs.disableModulation()
        self.srs2.disable_RFOutput(); self.srs2.disableModulation()
        self.srs3.disable_RFOutput(); self.srs3.disableModulation()
        self.srs4.disable_RFOutput(); self.srs4.disableModulation()
        if self.ifAWG: 
            self.AWG.turn_off()
        if self.ifAWG2: 
            self.AWG2.turn_off()

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T1RRDualNVObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.T1RRDualNVObject = measurementObject
        self.loopCounter = 0
        self.ifFancySpinInit = self.settings['ifFancySpinInit']
        self.tausArray = self.settings['tausArray']
        self.ifAWG = self.settings['ifAWG']
        self.ifAWG2 = self.settings['ifAWG2']

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];             
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration    = self.settings['laser_init_duration']
        laser_to_pi_delay       = self.settings['laser_to_pi_delay']
        pi_time                 = self.settings['pi_time'];                pi_time2               = self.settings['pi_time2'] 
        pi_to_ion_delay         = self.settings['pi_to_ion_delay'];        
        RRLaserSwitch_delay     = self.settings['RRLaserSwitch_delay'];    read_duration           = self.settings['read_duration']
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']; RRLaser2Switch_delay   = self.settings['RRLaser2Switch_delay']
        nSpinInit               = self.settings['nSpinInit']
        spinInit_RR_duration    = self.settings['spinInit_RR_duration'];   spinInit_RR_to_pi_delay= self.settings['spinInit_RR_to_pi_delay']
        pi_time3                = self.settings['pi_time3'];               spinInit_pi_to_RR_delay= self.settings['spinInit_pi_to_RR_delay']
        sweepWhich              = self.settings['sweepWhich'];             shift_btwn_2NV_read    = self.settings['shift_btwn_2NV_read']
        AWG_buffer              = self.settings['AWG_buffer'];             AWG_output_delay    = self.settings['AWG_output_delay']  
        AWG_buffer2             = self.settings['AWG_buffer2'];            AWG_output_delay2   = self.settings['AWG_output_delay2']         
        MW_duration_for_AWG = int(2*int((2*AWG_buffer + pi_time + 1)/2))
        MW_duration_for_AWG2 = int(2*int((2*AWG_buffer2 + pi_time2 + 1)/2))

        if sweepWhich == 'ti':
            ion_duration = tau_ns
        elif sweepWhich == 'nSpinInit':
            nSpinInit = tau_ns
        elif sweepWhich == 'pi_to_ion_delay':
            pi_to_ion_delay = int(tau_ns)

        self.num_loops = num_loops

        # Make pulse sequence
        pulse_sequence = []

        when_init_end           = laser_init_delay + laser_init_duration
        ######### The spin init block ###########
        spinInit_delay          = when_init_end + laser_to_pi_delay
        spinInit_duration_each  = spinInit_RR_duration + spinInit_RR_to_pi_delay + pi_time3 + spinInit_pi_to_RR_delay
        spinInit_duration_total = nSpinInit*spinInit_duration_each
        spinInit_lastRR_delay   = spinInit_delay + spinInit_duration_total
        when_spinInit_lastRR_end= spinInit_lastRR_delay + spinInit_RR_duration
        #########################################
        pi_delay                = when_spinInit_lastRR_end + spinInit_RR_to_pi_delay
        global MW_del;   MW_del = pi_delay + AWG_output_delay
        end1 = MW_del + MW_duration_for_AWG
        end2 = pi_delay + pi_time
        if (self.ifAWG==1) or (self.ifAWG2==1):
            when_pi_end = end1
        else:
            when_pi_end = end2
        
        laser_read_signal_delay    = when_pi_end + pi_to_ion_delay
        DAQ_signal_delay           = laser_read_signal_delay + RRLaserSwitch_delay
        DAQ_signal_duration        = read_duration
        when_DAQ_signal_end        = DAQ_signal_delay + DAQ_signal_duration 
        when_laser_read_signal_end = when_DAQ_signal_end + DAQ_to_laser_off_delay
        laser_read_signal_duration = when_laser_read_signal_end - laser_read_signal_delay

        laser_read_signal_NV2_delay = laser_read_signal_delay + shift_btwn_2NV_read
        DAQ_signal_NV2_delay        = laser_read_signal_NV2_delay + RRLaser2Switch_delay
        when_laser_read_signal_NV2_end = laser_read_signal_NV2_delay + laser_read_signal_duration
        when_DAQ_signal_NV2_end        = DAQ_signal_NV2_delay + DAQ_signal_duration

        when_everything_sig_end = np.max((when_laser_read_signal_NV2_end, when_DAQ_signal_NV2_end))
        ###################################################################################################################################
        ###################################################################################################################################
        ###################################################################################################################################
        laser_init_ref_delay = when_everything_sig_end + laser_init_delay
        laser_init_ref_duration = laser_init_duration

        when_init_ref_end       = laser_init_ref_delay + laser_init_ref_duration
        ######### The spin init block ###########
        spinInit_ref_delay           = when_init_ref_end + laser_to_pi_delay
        spinInit_duration_each       = spinInit_RR_duration + spinInit_RR_to_pi_delay + pi_time3 + spinInit_pi_to_RR_delay
        spinInit_duration_total      = nSpinInit*spinInit_duration_each
        spinInit_lastRR_ref_delay    = spinInit_ref_delay + spinInit_duration_total
        when_spinInit_lastRR_ref_end = spinInit_lastRR_ref_delay + spinInit_RR_duration
        #########################################
        pi_ref_delay            = when_spinInit_lastRR_ref_end + spinInit_RR_to_pi_delay
        MW_del_ref              = pi_ref_delay + AWG_output_delay
        if (self.ifAWG==1) or (self.ifAWG2==1):
            when_pi_ref_end = MW_del_ref + MW_duration_for_AWG
        else:
            when_pi_ref_end = pi_ref_delay + pi_time
        
        laser_read_ref_delay    = when_pi_ref_end + pi_to_ion_delay
        DAQ_ref_delay           = laser_read_ref_delay + RRLaserSwitch_delay
        DAQ_ref_duration        = read_duration
        when_DAQ_ref_end        = DAQ_ref_delay + DAQ_ref_duration
        when_laser_read_ref_end = when_DAQ_ref_end + DAQ_to_laser_off_delay
        laser_read_ref_duration = when_laser_read_ref_end - laser_read_ref_delay

        laser_read_ref_NV2_delay = laser_read_ref_delay + shift_btwn_2NV_read
        DAQ_ref_NV2_delay        = laser_read_ref_NV2_delay + RRLaser2Switch_delay

        if self.ifAWG:
            global ch1plot; global ch2plot
            ch1plot, ch2plot = AWG.send_T1_seq(pitime = int(pi_time), buffer=int(AWG_buffer))
        if self.ifAWG2:
            global ch1plot2; global ch2plot2
            ch1plot2, ch2plot2 = AWG2.send_T1_seq(pitime = int(pi_time2), buffer=int(AWG_buffer2))
        
##########################################################################################################
##########################################################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
        
        if self.ifFancySpinInit:
            for i in range(nSpinInit):
                delayRR = spinInit_delay + i*spinInit_duration_each
                delayMW2 = delayRR + spinInit_RR_duration + spinInit_RR_to_pi_delay
                pulse_sequence += [spc.Pulse('LaserRead',  delayRR,            duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch3',  delayMW2,           duration=int(pi_time3))] 
                pulse_sequence += [spc.Pulse('LaserRead2', delayRR,            duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch4',  delayMW2,           duration=int(pi_time3))]
            pulse_sequence += [spc.Pulse('LaserRead',  spinInit_lastRR_delay,  duration=int(spinInit_RR_duration))] 
            pulse_sequence += [spc.Pulse('LaserRead2', spinInit_lastRR_delay,  duration=int(spinInit_RR_duration))] 
        else:
            if self.ifAWG:
                pulse_sequence += [spc.Pulse('AWG',      pi_delay,      duration=pi_time)]
            else:
                pulse_sequence += [spc.Pulse('MWswitch',      pi_delay,            duration=int(pi_time))] 
            if self.ifAWG2:
                pulse_sequence += [spc.Pulse('AWG2',     pi_delay,     duration=pi_time)]
            else:
                pulse_sequence += [spc.Pulse('MWswitch2',     pi_delay,            duration=int(pi_time))]
        
        # Read sig NV1
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] 
        pulse_sequence += [spc.Pulse('Counter',      DAQ_signal_delay,        duration=int(DAQ_signal_duration))] 
        
        # Read sig NV2
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_signal_NV2_delay, duration=int(laser_read_signal_duration))] 
        pulse_sequence += [spc.Pulse('Counter',      DAQ_signal_NV2_delay,        duration=int(DAQ_signal_duration))] 
        
##############################################################################################################################
##############################################################################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay, duration=int(laser_init_ref_duration))] 
        
        if self.ifFancySpinInit:
            for i in range(nSpinInit):
                delayRR = spinInit_ref_delay + i*spinInit_duration_each
                delayMW2 = delayRR + spinInit_RR_duration + spinInit_RR_to_pi_delay
                pulse_sequence += [spc.Pulse('LaserRead',  delayRR,             duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch3',  delayMW2,            duration=int(pi_time3))]
                pulse_sequence += [spc.Pulse('LaserRead2', delayRR,             duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch4',  delayMW2,            duration=int(pi_time3))] 
            pulse_sequence += [spc.Pulse('LaserRead', spinInit_lastRR_ref_delay,duration=int(spinInit_RR_duration))] 
            pulse_sequence += [spc.Pulse('LaserRead2',spinInit_lastRR_ref_delay,duration=int(spinInit_RR_duration))] 
            if not self.ifAWG:
                pulse_sequence += [spc.Pulse('MWswitch',      pi_ref_delay,         duration=int(pi_time))]
            if not self.ifAWG2:
                pulse_sequence += [spc.Pulse('MWswitch2',     pi_ref_delay,         duration=int(pi_time))]
        
        # Read ref NV1
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay, duration=int(laser_read_ref_duration))] 
        pulse_sequence += [spc.Pulse('Counter',      DAQ_ref_delay,        duration=int(DAQ_ref_duration))] 

        # Read ref NV2
        pulse_sequence += [spc.Pulse('LaserRead2',    laser_read_ref_NV2_delay, duration=int(laser_read_ref_duration))] 
        pulse_sequence += [spc.Pulse('Counter',      DAQ_ref_NV2_delay,        duration=int(DAQ_ref_duration))] 
##############################################################################################################################
##############################################################################################################################

        self.read_duration = DAQ_signal_duration
        self.read_duration2 = DAQ_signal_duration
    
        self.pulse_sequence = pulse_sequence
        self.ifPrintTime = True
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False, ifPrintTime=self.ifPrintTime)
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
            samps_per_chan = 5*self.num_reads # x2 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        print('Set tau to ' + str(tau_ns) + " ns")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))
        self.ctrtask.stop(); self.ctrtask.close()
        rate = xLineData

        global sig_data;  sig_data  = rate[::self.num_reads_per_iter]
        global sig_data2; sig_data2 = rate[1::self.num_reads_per_iter]
        global ref_data;  ref_data  = rate[2::self.num_reads_per_iter]
        global ref_data2; ref_data2 = rate[3::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig_data)/(self.read_duration/1e9)/1e3
        global ref_avg;  ref_avg = np.average(ref_data)/(self.read_duration/1e9)/1e3
        global sig_avg2;  sig_avg2 = np.average(sig_data2)/(self.read_duration2/1e9)/1e3
        global ref_avg2;  ref_avg2 = np.average(ref_data2)/(self.read_duration2/1e9)/1e3

        return sig_avg
    
    def plotPulseSequences(self):
        self.readColor = self.settings['LaserRead']['channel']
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False, readColor=self.readColor)
                if (self.ifAWG) and (self.ifAWG2):
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, MW_del, label1='ch1-AWG1', label2='ch2-AWG1')
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot2,ch2plot2, MW_del,
                                                           fig=fig, offset1=17.75, offset2=17.5, label1='ch1-AWG2', label2='ch2-AWG2')
                else:
                    if self.ifAWG:
                        fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, MW_del, label1='ch1-AWG1', label2='ch2-AWG1')
                    elif self.ifAWG2:
                        fig = plotPulseObject.makePulsePlotAWG(ch1plot2,ch2plot2, MW_del,
                                                           offset1=17.75, offset2=17.5, label1='ch1-AWG2', label2='ch2-AWG2')
                    else:
                        fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T1RRDualNVObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
    
    def turn_on_at_end(self):
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
