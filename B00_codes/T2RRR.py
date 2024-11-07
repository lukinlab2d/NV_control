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


class T2RRR(Instrument):

    def __init__(self, name='T2RRRObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam, 'AWG': self.AWGParam,
                        'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        # Vpiezos, MW power, and MW frequency
        self.tausArray = self.settings['tausArray']
        self.SRSnum = self.settings['SRSnum'];      MWPower = self.settings['MWPower']; MWFreq = self.settings['MWFreq']
        self.SDGnum = self.settings['SDGnum']
        self.velNum = self.settings['velNum']; self.ifNeedVel = self.settings['ifNeedVel']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
       
        self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)

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
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(MWFreq) #Hz
        self.srs.set_RFAmplitude(MWPower) #dBm
        if (self.SRSnum != 3) and (self.SRSnum != 4):
            self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        # AWG object
        ifAWG = self.settings['ifAWG']; self.ifAWG = ifAWG
        if self.ifAWG:
            self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum)
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
            global vel; vel = self.vel

        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        global srs; srs = self.srs
    
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig, ref,
                                            qctask(sig.plotPulseSequences),
                                            ).then(qctask(sig.turn_on_at_end))

        data = loop.get_data_set(name='T2RRR')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.T2RRRObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.T2RRRObject_ref, name='ref')

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
        if self.ifAWG: self.AWG.turn_off()

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T2RRRObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.T2RRRObject = measurementObject
        self.loopCounter = 0
        self.settings = settings
        self.tausArray = self.settings['tausArray']

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']
        self.ifAWG = self.settings['ifAWG']


    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))

        self.ctrtask.stop(); self.ctrtask.close()

        rate = xLineData/(self.read_duration/1e9)/1e3
        sig = rate[::self.num_reads_per_iter]
        ref = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)

        return sig_avg

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))

        NO_MS_EQUALS_1 = 0
        Q_FINAL = 1
        THREE_PI_HALF_FINAL = 2
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];              
        laser_init_delay        = self.settings['laser_init_delay'];     laser_init_duration = self.settings['laser_init_duration']
        laser_to_MW_delay       = self.settings['laser_to_MW_delay'];    pi2time             = self.settings['pi_half']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];   read_duration       = self.settings['read_duration']   
        read_laser_duration     = self.settings['read_laser_duration'];  
        AWG_buffer              = self.settings['AWG_buffer'];           AWG_output_delay    = self.settings['AWG_output_delay']
        MW_duration_for_AWG = int(2*int((AWG_buffer + 2*pi2time + tau_ns + 1)/2))

        if tau_ns/2 > 30:
            MWI_to_switch_delay  = self.settings['MWI_to_switch_delay']
        else: 
            MWI_to_switch_delay  = 0

        ################# Signal  ################################

        when_init_end         = laser_init_delay + laser_init_duration
        MW_delay              = when_init_end + laser_to_MW_delay
        global MW_del; MW_del = MW_delay + AWG_output_delay
                     
        MW2_delay      = MW_delay + pi2time + tau_ns;         MW2_duration = pi2time

        if self.ifAWG:
            when_pulse_end = MW_del + MW_duration_for_AWG
        else:
            when_pulse_end = MW2_delay + MW2_duration     

        laser_read_signal_delay    = when_pulse_end
        laser_read_signal_duration = read_laser_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration

        ################# Reference  ################################
        laser_init_ref_delay = when_read_signal_end + laser_init_delay
        when_init_ref_end    = laser_init_ref_delay + laser_init_duration

        MW3_delay         = when_init_ref_end + laser_to_MW_delay;    MW3_duration = pi2time
        MW4_delay         = MW3_delay + MW3_duration + tau_ns;        MW4_duration = pi2time

        laser_read_ref_delay    = when_read_signal_end + laser_read_signal_delay
        laser_read_ref_duration = read_laser_duration
        read_ref_delay          = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration       = read_duration
        when_read_ref_end       = read_ref_delay + read_ref_duration

        sig_to_ref_wait         = laser_read_ref_delay - 2*MW_duration_for_AWG - MW_del
        self.read_duration      = read_signal_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")
        
        if self.ifAWG:
            global ch1plot; global ch2plot
            ch1plot, ch2plot = AWG.send_T2R_seq(pi_2time=int(pi2time), tau = int(tau_ns), 
                                                buffer=int(AWG_buffer), sig_to_ref_wait=int(sig_to_ref_wait))

        # Make pulse sequence
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit', laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit', laser_init_ref_delay,    duration=int(laser_init_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',     laser_read_signal_delay, duration=int(laser_read_signal_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',     laser_read_ref_delay,    duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('Counter',       read_signal_delay,       duration=int(read_signal_duration))] 
        pulse_sequence += [spc.Pulse('Counter',       read_ref_delay,          duration=int(read_ref_duration))] 

        if self.ifAWG:
            pulse_sequence += [spc.Pulse('AWG',      MW_delay,    duration=50)]
        else:
            pulse_sequence += [spc.Pulse('MWswitch', MW_delay,    duration=int(pi2time))]
            pulse_sequence += [spc.Pulse('MWswitch', MW2_delay,   duration=int(MW2_duration))]
        
            pulse_sequence += [spc.Pulse('MWswitch', MW3_delay,   duration=int(MW3_duration))]
            pulse_sequence += [spc.Pulse('MWswitch', MW4_delay,   duration=int(MW4_duration))]

            pulse_sequence += [spc.Pulse('MW_I', MW4_delay-MWI_to_switch_delay, duration=int(MW4_duration + 2*MWI_to_switch_delay))]
            pulse_sequence += [spc.Pulse('MW_Q', MW4_delay-MWI_to_switch_delay, duration=int(MW4_duration + 2*MWI_to_switch_delay))]
        
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
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, MW_del)
                else:
                    fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T2RRRObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
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