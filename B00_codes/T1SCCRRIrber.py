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

class T1SCCRRIrber(Instrument):

    def __init__(self, name='T1SCCRRIrberObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz, of the Pulse Blaster
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.LaserIonParam =    {'delay_time': 2, 'channel':settings['laserIon_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        # self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        # self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        # self.MWI2Param =         {'delay_time': 2, 'channel':settings['MWI2_channel']}
        # self.MWQ2Param =         {'delay_time': 2, 'channel':settings['MWQ2_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        self.hiLoMWPwrParam =    {'delay_time': 2, 'channel':settings['hiLoMWPwr_channel']}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 'AWG': self.AWGParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam, 'LaserIon': self.LaserIonParam,
                        'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'MWswitch2': self.MWswitch2Param, 'hiLoMWPwr': self.hiLoMWPwrParam,
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.tausArray = self.settings['tausArray']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)
        self.SDGnum = self.settings['SDGnum']; self.ifIQ = self.settings['ifIQ']; self.srate=self.settings['srate']
        self.SRSnum = self.settings['SRSnum'];   MWPower = self.settings['MWPower'];   MWFreq = self.settings['MWFreq']
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        self.ifMWDuringRead = self.settings['ifMWDuringRead']; self.ifMW2DuringRead = self.settings['ifMW2DuringRead']

        self.ifNeedVel = self.settings['ifNeedVel']; self.velNum = self.settings['velNum']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

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
            name = "sigFullData",
            parameter_class = SigFullData,
        )
        self.add_parameter(
            name = "refFullData",
            parameter_class = RefFullData,
        )
        self.savedPulseSequencePlots = {}

        # SRS object
        if self.ifMWDuringRead:
            self.srs = SRS(SRSnum=self.SRSnum)
            self.srs.set_freq(MWFreq) #Hz
            self.srs.set_RFAmplitude(MWPower) #dBm
            if (self.SRSnum != 3) and (self.SRSnum != 4) and self.ifIQ:
                self.srs.enableIQmodulation()
            self.srs.enable_RFOutput()
            global srs; srs = self.srs

        # SRS object 2
        if self.ifMW2DuringRead:
            self.srs2 = SRS(SRSnum=self.SRSnum2)
            self.srs2.set_freq(MWFreq2) #Hz
            self.srs2.set_RFAmplitude(MWPower2) #dBm
            if (self.SRSnum2 != 3) and (self.SRSnum2 != 4):
                self.srs2.enableIQmodulation()
            self.srs2.enable_RFOutput()
            global srs2; srs2 = self.srs2
        
        # AWG object 1
        ifAWG = self.settings['ifAWG']; self.ifAWG = ifAWG
        if self.ifAWG:
            self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum, srate=self.srate)
            global AWG; AWG = self.AWG

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
        
        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        
        
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigFullData = self.sigFullData
        refFullData = self.refFullData

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            ifIndividualCount = True, 
            numOfPBLoops = self.settings['num_loops'],
            sleepTimeAfterFinishing=0).each(sig,ref, sigFullData, refFullData,
                                            qctask(sig.plotPulseSequences),
                                            )
        data = loop.get_data_set(name='T1SCCRRIrber')
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
        if self.ifAWG: 
            self.AWG.turn_off()

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T1SCCRRIrberObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.T1SCCRRIrberObject = measurementObject
        self.loopCounter = 0
        self.ifMWDuringRead = self.settings['ifMWDuringRead']
        self.ifMW2DuringRead = self.settings['ifMW2DuringRead']
        self.ifAWG = self.settings['ifAWG']
        self.ifHiloExtra = self.settings['ifHiloExtra']
        self.ifMWReadLowDutyCycle = self.settings['ifMWReadLowDutyCycle']
        self.ifFancySpinInit = self.settings['ifFancySpinInit']
        self.tausArray = self.settings['tausArray']

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];             
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration    = self.settings['laser_init_duration']
        laser_to_pi_delay       = self.settings['laser_to_pi_delay'];      pi_time                = self.settings['pi_time']
        pi_to_ion_delay         = self.settings['pi_to_ion_delay'];        ion_duration           = self.settings['ion_duration']
        RRLaserSwitch_delay     = self.settings['RRLaserSwitch_delay'];    DAQ_duration           = self.settings['DAQ_duration']
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']
        laserRead_to_MWmix      = self.settings['laserRead_to_MWmix'];     ion_to_read_delay      = self.settings['ion_to_read_delay']
        iznLaserSwitch_delay    = self.settings['iznLaserSwitch_delay'];   pi_to_hilo_extra_delay = self.settings['pi_to_hilo_extra_delay']
        MWmix_duration_short    = self.settings['MWmix_duration_short'];   delay_between_MWmix    = self.settings['delay_between_MWmix']
        nSpinInit               = self.settings['nSpinInit']
        spinInit_RR_duration    = self.settings['spinInit_RR_duration'];   spinInit_RR_to_pi_delay= self.settings['spinInit_RR_to_pi_delay']
        pi_time2                = self.settings['pi_time2'];               spinInit_pi_to_RR_delay= self.settings['spinInit_pi_to_RR_delay']
        sweepWhich              = self.settings['sweepWhich']
        AWG_buffer              = self.settings['AWG_buffer'];             AWG_output_delay       = self.settings['AWG_output_delay']  
        

        if sweepWhich == 'ti':
            ion_duration = tau_ns
        elif sweepWhich == 'MWmix_duration_short':
            MWmix_duration_short = tau_ns
        elif sweepWhich == 'delay_between_MWmix':
            delay_between_MWmix = tau_ns
        elif sweepWhich == 'nSpinInit':
            nSpinInit = tau_ns
        elif sweepWhich == 'pi_to_ion_delay':
            pi_to_ion_delay = int(tau_ns)

        if self.ifMWReadLowDutyCycle: nMWread = int(DAQ_duration/(MWmix_duration_short + delay_between_MWmix))
        self.num_loops = num_loops

        # Make pulse sequence
        pulse_sequence = []

        when_init_end           = laser_init_delay + laser_init_duration
        ######### The spin init block ###########
        spinInit_delay          = when_init_end + laser_to_pi_delay
        spinInit_duration_each  = spinInit_RR_duration + spinInit_RR_to_pi_delay + pi_time2 + spinInit_pi_to_RR_delay
        spinInit_duration_total = nSpinInit*spinInit_duration_each
        spinInit_lastRR_delay   = spinInit_delay + spinInit_duration_total
        when_spinInit_lastRR_end= spinInit_lastRR_delay + spinInit_RR_duration
        #########################################
        pi_delay                = when_spinInit_lastRR_end + spinInit_RR_to_pi_delay; self.pi_delay = pi_delay
        when_pi_end             = pi_delay + pi_time
        ion_RR_delay            = when_pi_end + pi_to_ion_delay                 
        when_ion_RR_end         = ion_RR_delay + ion_duration
        ion_strong_delay        = ion_RR_delay + RRLaserSwitch_delay - iznLaserSwitch_delay
        when_ion_strong_end     = ion_strong_delay + ion_duration
        
        laser_read_signal_delay    = when_ion_strong_end + ion_to_read_delay
        DAQ_signal_delay           = laser_read_signal_delay + RRLaserSwitch_delay
        DAQ_signal_duration        = DAQ_duration
        when_DAQ_signal_end        = DAQ_signal_delay + DAQ_signal_duration 
        when_laser_read_signal_end = when_DAQ_signal_end + DAQ_to_laser_off_delay
        laser_read_signal_duration = when_laser_read_signal_end - laser_read_signal_delay

        MWmix_delay                = laser_read_signal_delay + laserRead_to_MWmix
        MWmix_duration             = laser_read_signal_duration - laserRead_to_MWmix + RRLaserSwitch_delay
        # hilo_delay                 = when_pi_end + 2e2                        #MWmix_delay - 0.2*ion_to_read_delay
        # hilo_duration              = when_laser_read_signal_end - hilo_delay #MWmix_duration + 0.4*ion_to_read_delay
        hilo_delay                 = laser_read_signal_delay #MWmix_delay - 0.2*ion_to_read_delay
        hilo_duration              = laser_read_signal_duration #MWmix_duration + 0.4*ion_to_read_delay

        hilo_extra_delay           = when_pi_end + pi_to_hilo_extra_delay
        when_hilo_extra_end        = ion_RR_delay - 1e3
        hilo_extra_duration        = when_hilo_extra_end - hilo_extra_delay
        ###################################################################################################################################
        ###################################################################################################################################
        ###################################################################################################################################

        laser_init_ref_delay = (when_laser_read_signal_end + RRLaserSwitch_delay) + laser_init_delay
        laser_init_ref_duration = laser_init_duration

        when_init_ref_end       = laser_init_ref_delay + laser_init_ref_duration
        ######### The spin init block ###########
        spinInit_ref_delay           = when_init_ref_end + laser_to_pi_delay
        spinInit_duration_each       = spinInit_RR_duration + spinInit_RR_to_pi_delay + pi_time2 + spinInit_pi_to_RR_delay
        spinInit_duration_total      = nSpinInit*spinInit_duration_each
        spinInit_lastRR_ref_delay    = spinInit_ref_delay + spinInit_duration_total
        when_spinInit_lastRR_ref_end = spinInit_lastRR_ref_delay + spinInit_RR_duration
        #########################################
        pi_ref_delay            = when_spinInit_lastRR_ref_end + spinInit_RR_to_pi_delay
        when_pi_ref_end         = pi_ref_delay + pi_time
        ion_ref_RR_delay        = when_pi_ref_end + pi_to_ion_delay
        when_ion_ref_RR_end     = ion_ref_RR_delay + ion_duration
        ion_ref_strong_delay    = ion_ref_RR_delay + RRLaserSwitch_delay - iznLaserSwitch_delay
        when_ion_ref_strong_end = ion_ref_strong_delay + ion_duration

        laser_read_ref_delay    = when_ion_ref_strong_end + ion_to_read_delay
        DAQ_ref_delay           = laser_read_ref_delay + RRLaserSwitch_delay
        DAQ_ref_duration        = DAQ_duration
        when_DAQ_ref_end        = DAQ_ref_delay + DAQ_ref_duration
        when_laser_read_ref_end = when_DAQ_ref_end + DAQ_to_laser_off_delay
        laser_read_ref_duration = when_laser_read_ref_end - laser_read_ref_delay

        MWmix_ref_delay         = laser_read_ref_delay + laserRead_to_MWmix
        MWmix_ref_duration      = laser_read_ref_duration - laserRead_to_MWmix + RRLaserSwitch_delay
        # hilo_ref_delay          = when_pi_ref_end + 2e2                     #MWmix_ref_delay - 0.2*ion_to_read_delay
        # hilo_ref_duration       = when_laser_read_ref_end - hilo_ref_delay #MWmix_ref_duration + 0.4*ion_to_read_delay
        hilo_ref_delay          = laser_read_ref_delay#MWmix_ref_delay - 0.2*ion_to_read_delay
        hilo_ref_duration       = laser_read_ref_duration#MWmix_ref_duration + 0.4*ion_to_read_delay

        hilo_extra_ref_delay    = when_pi_ref_end + pi_to_hilo_extra_delay
        when_hilo_extra_ref_end = ion_ref_RR_delay - 1e3

        MW_delay_for_AWG        = pi_delay - AWG_output_delay
        sig_to_ref_wait         = MWmix_ref_delay - MWmix_delay - MWmix_duration
        pi_to_MWmix_wait        = MWmix_delay - pi_delay - pi_time

        if self.loopCounter==0: sleepTime = 12
        else: sleepTime = 6
        if self.ifAWG:
            global ch1plot; global ch2plot
            ch1plot, ch2plot = AWG.send_SCCRRPhotonStatIrber_seq(int(MWmix_duration), int(AWG_buffer),
                            int(sig_to_ref_wait), pitime=int(pi_time), pi_to_MWmix_wait=int(pi_to_MWmix_wait),
                            sleepTime=sleepTime)
        
        
##########################################################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
        
        if self.ifFancySpinInit:
            for i in range(nSpinInit):
                delayRR = spinInit_delay + i*spinInit_duration_each
                delayMW2 = delayRR + spinInit_RR_duration + spinInit_RR_to_pi_delay
                pulse_sequence += [spc.Pulse('LaserRead', delayRR,            duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch2', delayMW2,           duration=int(pi_time2))] 
            pulse_sequence += [spc.Pulse('LaserRead', spinInit_lastRR_delay,  duration=int(spinInit_RR_duration))] 
        else:
            if self.ifAWG:
                pulse_sequence += [spc.Pulse('AWG',      MW_delay_for_AWG, duration=50)]
            else:
                pulse_sequence += [spc.Pulse('MWswitch',     pi_delay,            duration=int(pi_time))] 
        
        pulse_sequence += [spc.Pulse('LaserIon',     ion_strong_delay,        duration=int(ion_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',    ion_RR_delay,            duration=int(ion_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] 
        
        if self.ifMWDuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch',  delay,          duration=int(MWmix_duration_short))] 
            else:
                if self.ifAWG == 0:
                    pulse_sequence += [spc.Pulse('MWswitch', MWmix_delay,         duration=int(MWmix_duration))]
        if self.ifMW2DuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch2',  delay,         duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch2',MWmix_delay,         duration=int(MWmix_duration))]

        if self.ifMWReadLowDutyCycle == 0:
            pulse_sequence += [spc.Pulse('hiLoMWPwr',    hilo_delay,         duration=int(hilo_duration))]
        
        pulse_sequence += [spc.Pulse('Counter',      DAQ_signal_delay,        duration=int(DAQ_signal_duration))] 
        if hilo_extra_duration >= 1e2 and self.ifHiloExtra:
            pulse_sequence += [spc.Pulse('hiLoMWPwr',    hilo_extra_delay,         duration=int(hilo_extra_duration))]
        #########################################################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay, duration=int(laser_init_ref_duration))] 
        
        if self.ifFancySpinInit:
            for i in range(nSpinInit):
                delayRR = spinInit_ref_delay + i*spinInit_duration_each
                delayMW2 = delayRR + spinInit_RR_duration + spinInit_RR_to_pi_delay
                pulse_sequence += [spc.Pulse('LaserRead', delayRR,             duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch2', delayMW2,            duration=int(pi_time2))] 
            pulse_sequence += [spc.Pulse('LaserRead',spinInit_lastRR_ref_delay,duration=int(spinInit_RR_duration))] 
            pulse_sequence += [spc.Pulse('MWswitch',     pi_ref_delay,         duration=int(pi_time))]
        
        pulse_sequence += [spc.Pulse('LaserIon',     ion_ref_strong_delay, duration=int(ion_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',    ion_ref_RR_delay,     duration=int(ion_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay, duration=int(laser_read_ref_duration))] 
        
        if self.ifMWDuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_ref_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch',  delay,       duration=int(MWmix_duration_short))] 
            else:
                if self.ifAWG == 0:
                    pulse_sequence += [spc.Pulse('MWswitch', MWmix_ref_delay,  duration=int(MWmix_ref_duration))]
        if self.ifMW2DuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_ref_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch2', delay,       duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch2',MWmix_ref_delay,  duration=int(MWmix_ref_duration))]
        
        if self.ifMWReadLowDutyCycle == 0:
            pulse_sequence += [spc.Pulse('hiLoMWPwr',    hilo_ref_delay,  duration=int(hilo_ref_duration))]
        
        pulse_sequence += [spc.Pulse('Counter',      DAQ_ref_delay,        duration=int(DAQ_ref_duration))] 
        
        if hilo_extra_duration >= 1e2 and self.ifHiloExtra:
            pulse_sequence += [spc.Pulse('hiLoMWPwr',hilo_extra_ref_delay, duration=int(hilo_extra_duration))]

        self.read_duration = DAQ_signal_duration
    
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

        global sig_data; sig_data = rate[::self.num_reads_per_iter]
        global ref_data; 
        ref_data = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig_data)
        global ref_avg;  ref_avg = np.average(ref_data)

        return sig_avg
    
    def plotPulseSequences(self):
        self.readColor = self.settings['LaserRead']['channel']
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False, readColor=self.readColor)
                if self.ifAWG:
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, self.pi_delay, label1='chnl1-AWG1', label2='chnl2-AWG1')
                else:
                    fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T1SCCRRIrberObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
    
    def turn_on_at_end(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        pb.turn_on_infinite(channels=channels)


class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg

class SigFullData(Parameter):
    def __init__(self, name='sigFullData',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_data
    
class RefFullData(Parameter):
    def __init__(self, name='refFullData',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_data