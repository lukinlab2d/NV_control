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

class CalibRedRRIznRateWithIznPulse(Instrument):

    def __init__(self, name='CalibRedRRIznRateWithIznPulseObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz, of the Pulse Blaster
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.LaserIonParam =    {'delay_time': 2, 'channel':settings['laserIon_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.MWI2Param =         {'delay_time': 2, 'channel':settings['MWI2_channel']}
        self.MWQ2Param =         {'delay_time': 2, 'channel':settings['MWQ2_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        self.hiLoMWPwrParam =    {'delay_time': 2, 'channel':settings['hiLoMWPwr_channel']}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 'AWG': self.AWGParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam, 'LaserIon': self.LaserIonParam,
                        'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'MW_I2': self.MWI2Param, 'MW_Q2': self.MWQ2Param, 'MWswitch2': self.MWswitch2Param,'hiLoMWPwr': self.hiLoMWPwrParam,
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.tausArray = np.linspace(start, stop, num_sweep_points)

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)

        self.SRSnum = self.settings['SRSnum'];   MWPower = self.settings['MWPower'];   MWFreq = self.settings['MWFreq']
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        self.ifMWDuringRead = self.settings['ifMWDuringRead']; self.ifMW2DuringRead = self.settings['ifMW2DuringRead']
        self.SDGnum = self.settings['SDGnum']; self.srate=self.settings['srate']

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
        self.savedPulseSequencePlots = {}

        # SRS object
        if self.ifMWDuringRead:
            self.srs = SRS(SRSnum=self.SRSnum)
            self.srs.set_freq(MWFreq) #Hz
            self.srs.set_RFAmplitude(MWPower) #dBm
            if (self.SRSnum != 3) and (self.SRSnum != 4):
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

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            ifIndividualCount = False, 
            numOfPBLoops = self.settings['num_loops'],
            ifMultReads=1, nread=self.settings['num_reads'],
            sleepTimeAfterFinishing=0).each(
                                            sig, ref, 
                                            qctask(sig.plotPulseSequences),
                                            )
        data = loop.get_data_set(name='CalibRedRRIznRateWithIznPulse')
        data.add_metadata(self.settings)
        self.data = data
        
        loop.run()
        print('Data saved to ' + str(data.location) + '/')
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)
        
        if self.ifMWDuringRead:
            self.srs.disable_RFOutput()
            self.srs.disableModulation()
        if self.ifMW2DuringRead:
            self.srs2.disable_RFOutput()
            self.srs2.disableModulation()
        if self.ifAWG: 
            self.AWG.turn_off()

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/CalibRedRRIznRateWithIznPulseObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.CalibRedRRIznRateWithIznPulseObject = measurementObject
        self.loopCounter = 0
        self.ifIonizedRef = self.settings['ifIonizedRef']
        self.ifMWDuringRead = self.settings['ifMWDuringRead']
        self.ifMW2DuringRead = self.settings['ifMW2DuringRead']
        self.ifAWG = self.settings['ifAWG']
        self.ifPi = self.settings['ifPi']
        self.ifIznWithRR = self.settings['ifIznWithRR']
        self.ifIznWithMW = self.settings['ifIznWithMW']
        self.ifMWReadLowDutyCycle = self.settings['ifMWReadLowDutyCycle']
        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.tausArray = np.linspace(start, stop, num_sweep_points)

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];             
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration     = self.settings['laser_init_duration']
        laser_to_ion_delay      = self.settings['laser_to_ion_delay'];     ion_duration            = self.settings['ion_duration']
        RRLaserSwitch_delay     = self.settings['RRLaserSwitch_delay'];    read_duration           = self.settings['read_duration']
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']
        delay_between_reads     = self.settings['delay_between_reads'];    nread                   = self.settings['num_reads']
        laserRead_to_MWmix      = self.settings['laserRead_to_MWmix'];     ion_to_read_delay       = self.settings['ion_to_read_delay']
        iznLaserSwitch_delay    = self.settings['iznLaserSwitch_delay']
        MWmix_duration_short    = self.settings['MWmix_duration_short'];   delay_between_MWmix     = self.settings['delay_between_MWmix']
        nMWread                 = self.settings['nMWread']
        AWG_buffer              = self.settings['AWG_buffer'];             AWG_output_delay        = self.settings['AWG_output_delay']  
        pi_time                 = self.settings['pi_time'];                laser_to_pi_delay       = self.settings['laser_to_pi_delay']
        pi_to_ion_delay= laser_to_ion_delay-laser_to_pi_delay-pi_time 
        self.num_loops = num_loops

        # Make pulse sequence
        pulse_sequence = []

        when_init_end      = laser_init_delay + laser_init_duration
        ion_RR_delay       = when_init_end + laser_to_ion_delay                 
        when_ion_RR_end    = ion_RR_delay + ion_duration
        ion_strong_delay   = ion_RR_delay + RRLaserSwitch_delay - iznLaserSwitch_delay
        when_ion_strong_end= ion_strong_delay + ion_duration
        
        laser_read_signal_delay    = when_ion_strong_end + ion_to_read_delay
        first_read_signal_delay    = laser_read_signal_delay + RRLaserSwitch_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = first_read_signal_delay + tau_ns
        when_laser_read_signal_end = when_read_signal_end + DAQ_to_laser_off_delay
        laser_read_signal_duration = when_laser_read_signal_end - laser_read_signal_delay

        MWmix_delay                = laser_read_signal_delay + laserRead_to_MWmix; self.MWmix_delay = MWmix_delay
        MWmix_duration             = laser_read_signal_duration - laserRead_to_MWmix + RRLaserSwitch_delay

        MWmix_delay_for_AWG        = MWmix_delay - AWG_output_delay
        MWmix_duration_for_AWG     = int(2*int((2*AWG_buffer + MWmix_duration + 1)/2)) # to make it even
        ###########################################################################
        laser_init_ref_delay = (when_laser_read_signal_end + RRLaserSwitch_delay) + laser_init_delay
        laser_init_ref_duration = laser_init_duration

        when_init_ref_end       = laser_init_ref_delay + laser_init_ref_duration
        pi_ref_delay            = when_init_ref_end + laser_to_pi_delay
        
        ion_ref_RR_delay        = when_init_ref_end + laser_to_ion_delay
        when_ion_ref_RR_end     = ion_ref_RR_delay + ion_duration
        ion_ref_strong_delay    = ion_ref_RR_delay + RRLaserSwitch_delay - iznLaserSwitch_delay
        when_ion_ref_strong_end = ion_ref_strong_delay + ion_duration

        laser_read_ref_delay    = when_ion_ref_strong_end + ion_to_read_delay
        first_read_ref_delay    = laser_read_ref_delay + RRLaserSwitch_delay
        read_ref_duration       = read_duration
        when_read_ref_end       = first_read_ref_delay + tau_ns
        when_laser_read_ref_end = when_read_ref_end + DAQ_to_laser_off_delay
        laser_read_ref_duration = when_laser_read_ref_end - laser_read_ref_delay

        MWmix_ref_delay         = laser_read_ref_delay + laserRead_to_MWmix
        MWmix_ref_duration      = laser_read_ref_duration - laserRead_to_MWmix + RRLaserSwitch_delay

        MWmix_ref_delay_for_AWG        = MWmix_ref_delay - AWG_output_delay
        MWmix_ref_duration_for_AWG     = int(2*int((2*AWG_buffer + MWmix_ref_duration + 1)/2)) # to make it even

        sig_to_ref_wait = MWmix_ref_delay - MWmix_delay - MWmix_duration
        sig_to_pi_wait  = pi_ref_delay - MWmix_delay - MWmix_duration

        if self.loopCounter==0: sleepTime = 12
        else: sleepTime = 6
        if self.ifAWG:
            global ch1plot; global ch2plot
            if self.ifIonizedRef==0:
                ch1plot, ch2plot = AWG.send_fastRabi_seq(pulse_width=int(MWmix_duration), buffer=int(AWG_buffer))
            else:
                ch1plot, ch2plot = AWG.send_CalibRedRRIznRateWithIznPulse_seq(int(MWmix_duration), int(AWG_buffer),
                                    int(sig_to_ref_wait), ifPi=self.ifPi, sig_to_pi_wait=int(sig_to_pi_wait), pi_time=int(pi_time),
                                    sleepTime=sleepTime)
##########################################################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] 
        if self.ifMWDuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch',  delay,          duration=int(MWmix_duration_short))] 
            else:
                if self.ifAWG:
                    pulse_sequence += [spc.Pulse('AWG',      MWmix_delay_for_AWG, duration=50)]
                else:
                    pulse_sequence += [spc.Pulse('MWswitch', MWmix_delay,         duration=int(MWmix_duration))]
        if self.ifMW2DuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch2',  delay,         duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch2', MWmix_delay,        duration=int(MWmix_duration))]
        if self.ifMWReadLowDutyCycle == 0:
            pulse_sequence += [spc.Pulse('hiLoMWPwr',    MWmix_delay,         duration=int(MWmix_duration))]

        for i in range(nread):
            delay = first_read_signal_delay + i*(read_signal_duration + delay_between_reads)
            pulse_sequence += [spc.Pulse('Counter',  delay,             duration=int(read_signal_duration))] 
        
        #####################
        if self.ifIonizedRef:
            if not laser_init_delay == 0:
                pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay, duration=int(laser_init_ref_duration))] 
            pulse_sequence += [spc.Pulse('LaserIon',     ion_ref_strong_delay, duration=int(ion_duration))] 
            if self.ifIznWithRR:
                pulse_sequence += [spc.Pulse('LaserRead',ion_ref_RR_delay,     duration=int(ion_duration))] 
            if self.ifIznWithMW:
                pulse_sequence += [spc.Pulse('MWswitch', ion_ref_strong_delay, duration=int(ion_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch2',ion_ref_strong_delay, duration=int(ion_duration))] 
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
                    pulse_sequence += [spc.Pulse('MWswitch2', MWmix_ref_delay, duration=int(MWmix_ref_duration))]
            if self.ifMWReadLowDutyCycle == 0:
                    pulse_sequence += [spc.Pulse('hiLoMWPwr', MWmix_ref_delay, duration=int(MWmix_ref_duration))]
            
            for i in range(nread):
                delay = first_read_ref_delay + i*(read_ref_duration + delay_between_reads)
                pulse_sequence += [spc.Pulse('Counter',  delay,                duration=int(read_ref_duration))] 
        


        self.read_duration = read_signal_duration
    
        self.pulse_sequence = pulse_sequence
        self.ifPrintTime = True
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False, ifPrintTime=self.ifPrintTime)
        self.pb.program_pb(pulse_sequence, num_loops=num_loops)

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        self.num_reads = int(num_loops * num_reads_per_iter)
        print('Num of reads ', self.num_reads)
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

        print('Set tau to ' + str(tau_ns/1e6) + " ms")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))
        self.ctrtask.stop(); self.ctrtask.close()
        rate = xLineData

        global sig_avg;  global ref_avg
        if self.ifIonizedRef == 0:
            global sig_data; 
            sig_data = np.zeros((self.num_reads_per_iter, self.num_loops))
            for i in range(self.num_reads_per_iter):
                each_read_data = rate[i::self.num_reads_per_iter]
                sig_data[i] = each_read_data
            
            sig_avg = np.average(sig_data, 1)
            ref_avg = np.zeros(len(sig_avg))
        elif self.ifIonizedRef == 1:
            global all_data; 
            all_data = np.zeros((self.num_reads_per_iter, self.num_loops))
            for i in range(self.num_reads_per_iter):
                each_read_data = rate[i::self.num_reads_per_iter]
                all_data[i] = each_read_data
            all_avg = np.average(all_data, 1)
            sig_avg = all_avg[0:int(len(all_avg)/2)]
            ref_avg = all_avg[int(len(all_avg)/2):]

        return sig_avg
    
    def plotPulseSequences(self):
        self.readColor = self.settings['LaserRead']['channel']
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False, readColor=self.readColor)
                if self.ifAWG:
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, self.MWmix_delay, label1='chnl1-AWG1', label2='chnl2-AWG1')
                else:
                    fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.CalibRedRRIznRateWithIznPulseObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
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