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
import pickle
from qcodes_contrib_drivers.drivers.NationalInstruments.DAQ import *
from qcodes_contrib_drivers.drivers.NationalInstruments.class_file import *

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
    

class ConfocalXY8(Instrument):

    def __init__(self, name='ConfocalXY8Object', settings=None, ifPlotPulse=True, **kwargs) -> None:

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
    
        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,
                          'AWG': self.AWGParam,'PB_type': 'USB',
                          'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']
        self.ifPlotPulse  = ifPlotPulse

        # List of frequencies and power
        self.SRSnum=self.settings['SRSnum']; MWPower = self.settings['MWPower']
        self.SDGnum=self.settings['SDGnum']; self.srate=self.settings['srate']

        self.xArray = self.settings['xArray']; self.yArray = self.settings['yArray']
        
        # SRS object
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(3e9) #Hz
        self.srs.set_RFAmplitude(MWPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        # AWG object
        self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum, srate=self.srate)
        global AWG; AWG = self.AWG

        self.add_parameter(
            name = "sig",
            parameter_class = Signal,
            settings = self.settings,
            measurementObject = self
        )
        self.add_parameter(
            name = "ref",
            parameter_class = Reference,
        )
        self.add_parameter(
            name = "sigStd",
            parameter_class = SigStd,
        )
        self.add_parameter(
            name = "refStd",
            parameter_class = RefStd,
        )
        self.add_parameter(
            name = "y",
            settings = self.settings,
            parameter_class = VoltageY,
        )
        self.add_parameter(
            name = "x",
            settings = self.settings,
            parameter_class = VoltageX,
            measurementObject = self
        )

        # Make Pulse Blaster, Counter, SRS global objects
        global srs; srs = self.srs

        # DAQ object to control galvo
        galvo_card_name = "cDAQ1Mod2"
        galvo_ao_channels = {f'{galvo_card_name}/ao{i}': i for i in range(4)} # dictionary of analog output channels
        self.galvo = DAQAnalogOutputs("ConfocalXY8DAQ", galvo_card_name, galvo_ao_channels)
        global galvo; galvo = self.galvo
    
    def runScan(self):
        y = self.y; x = self.x

        # Loop through x and y
        loopx = Loop(
            x.sweep(keys=self.xArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(x.loop)
        loopy = Loop(
            y.sweep(keys=self.yArray),
            delay=0,
            sleepTimeAfterFinishing=0).each(loopx)

        data = loopy.get_data_set(name='ConfocalXY8')
        data.add_metadata(self.settings)
        self.data = data

        loopy.run()
        print('Data saved to ' + str(data.location) + '/')
        
        if self.ifPlotPulse:
            pulsePlotFilename = data.location + "/pulsePlot.png"
            plotPulseObject = PlotPulse(measurementObject=self, plotFilename=pulsePlotFilename, ifShown=True)
            plotPulseObject.makePulsePlot(ch1plot, ch2plot, MW_del)
        
        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        self.AWG.turn_off()
        self.galvo.close() 
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ConfocalXY8Object_sig_set.dat'


class Signal(Parameter):
    def __init__(self, name='sig', settings=None, measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.loopCounter = 0
        self.ConfocalXY8Object = measurementObject
        self.srate = self.settings['srate']

    def get_raw(self):      
        ctrtask.start()

        pb.start_pulse_seq()
        pb.wait()
        pb.stop_pulse_seq_without_closing()

        data = np.array(ctrtask.read(self.num_reads, timeout=10))

        ctrtask.stop()#; ctrtask.close()
        
        rate = data/(self.read_duration/1e9)/1e3
        sig = rate[::self.num_reads_per_iter]
        ref = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)
        global sig_std;  sig_std = np.std(sig)
        global ref_std;  ref_std = np.std(ref)

        return sig_avg
    
    def set_raw(self, tau_ns):
        print("Loop " + str(self.loopCounter))

        pitime = self.settings['pitime']; pi_incr_factor = self.settings['pi_incr_factor']
        pitime2 = self.settings['pitime2']; pi_incr_factor2 = self.settings['pi_incr_factor2']
        ifIncPower = self.settings['ifIncPower']; MWPower = self.settings['MWPower']

        xref_pitime = self.settings['xref_pitime']
        pi_increment = np.round((x_current-(xref_pitime))*pi_incr_factor)
        pi_increment = int(2*int(pi_increment/2))
        
        pi_increment2 = np.round((x_current-(xref_pitime))*pi_incr_factor2)
        pi_increment2 = int(2*int(pi_increment2/2))

        MWPower_incr = (x_current-(xref_pitime))*pi_incr_factor

        global num_loops; global read_duration
        
        # Pulse parameters
        rest_after_first_pulse=0
        num_loops               = self.settings['num_loops'];              mode = self.settings['mode']
        numxy8                  = self.settings['numxy8']
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration = self.settings['laser_init_duration']
        laser_to_AWG_delay      = self.settings['laser_to_AWG_delay'];     pi2time             = self.settings['pi2time']
        pitime                  = self.settings['pitime'];                 MW_to_read_delay    = self.settings['MW_to_read_delay']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];     read_duration       = self.settings['read_duration']   
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']; AWG_output_delay    = self.settings['AWG_output_delay']
        AWG_buffer              = self.settings['AWG_buffer'];             phi_IQ              = self.settings['phi_IQ']
        
        if ifIncPower==1:
            MWPower = MWPower + MWPower_incr
            srs.set_RFAmplitude(MWPower) #dBm
            pi_increment = 0; pi_increment2 = 0

        pitime = int(pitime + pi_increment);    pitime2  = int(pitime2 + pi_increment2)
        pi2time= int(pi2time + pi_increment/2); pi2time2 = int(pitime2/2)
        if ifOnFlake==0:
            pitime=pitime2; pi2time=pi2time2

        if np.mod(self.loopCounter,len(self.ConfocalXY8Object.tausArray))==0:    
            print('pi_increment2 = ' + str(pi_increment2))
            print('pi_increment = ' + str(pi_increment))
            print('pi_time = ' + str(pitime))
            print('MW power = ' + str(MWPower))
            print('-----o-o-o-o-o-o-----')

        #######################################################
        MW_duration     = int(2*int((AWG_buffer + 2*pi2time + numxy8*tau_ns*8 + numxy8*pitime*8 + 1 + (pitime+rest_after_first_pulse ))/2))
        # the last (pitime+rest_after_first_pulse ) should be updated if we prepare ms=0 vs ms=+-1

        when_init_end = laser_init_delay + laser_init_duration
        
        temp = when_init_end + laser_to_AWG_delay - AWG_output_delay
        if temp > 0:
            AWG_trig_delay = temp
        else:
            AWG_trig_delay = 0    
        
        when_sigMW_end = AWG_output_delay + AWG_trig_delay + MW_duration 
        rounded_when_sigMW_end = np.ceil(when_sigMW_end / 10) * 10

        global MW_del; MW_del = AWG_trig_delay+AWG_output_delay
        
        laser_read_signal_delay    = rounded_when_sigMW_end + MW_to_read_delay
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        laser_read_signal_duration = when_read_signal_end + DAQ_to_laser_off_delay - laser_read_signal_delay
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        
        MW_del_ref      = when_laser_read_signal_end + MW_del
        sig_to_ref_wait = MW_del_ref - when_sigMW_end 

        laser_init_ref_delay = when_laser_read_signal_end + laser_init_delay
        when_init_ref_end    = laser_init_ref_delay + laser_init_duration

        laser_read_ref_delay = when_laser_read_signal_end + laser_read_signal_delay
        read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay;  read_ref_duration    = read_duration; 
        when_read_ref_end    = read_ref_delay + read_ref_duration
        laser_read_ref_duration = when_read_ref_end + DAQ_to_laser_off_delay - laser_read_ref_delay
        self.read_duration = read_signal_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")

        #  Set I and Q for special pi/2 pulse
        I = round(32767*np.cos(phi_IQ))
        Q = round(32767*np.sin(phi_IQ))

        if self.srate is not None:
            if self.loopCounter==0: sleepTime = 10
            else: sleepTime = 1
        else:
            sleepTime = 0

        global ch1plot; global ch2plot
        ch1plot, ch2plot = AWG.send_XY8_seq(mode = mode, numxy8=numxy8, pi_2time=int(pi2time),rest_after_first_pulse=rest_after_first_pulse,
                                            pitime = int(pitime), tau = int(tau_ns),AWGnum=2,
                                              buffer=int(AWG_buffer), sig_to_ref_wait=int(sig_to_ref_wait))

        # Make pulse sequence
        global pulse_sequence; pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,              duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay,          duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,          duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('AWG',          AWG_trig_delay,                duration=50)]
        pulse_sequence += [spc.Pulse('Counter',  read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
           
        self.pulse_sequence = pulse_sequence
        global pb
        pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        pb.program_pb(pulse_sequence, num_loops=num_loops)

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        self.num_reads = int(num_loops * num_reads_per_iter)
        self.num_reads_per_iter = num_reads_per_iter

        # Pulse width counter. Timebase = signal; gate = PB signal
        global ctrtask
        ctrtask = nidaqmx.Task()
        pulseWidthChan = ctrtask.ci_channels.add_ci_pulse_width_chan( # define the pulse width counter
            counter = "cDAQ1Mod1/ctr0",
            name_to_assign_to_channel = "",
            min_val = 0,
            max_val = int(1e8),
            units = TimeUnits.TICKS,
            starting_edge = Edge.RISING,
            )
        ctrtask.timing.cfg_implicit_timing(
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
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False)
                fig = plotPulseObject.makePulsePlot(ch1plot, ch2plot, MW_del)
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.XY8AWGObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
            
    def close(self):
        ctrtask.close()
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        pb.turn_on_infinite(channels=channels)

class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg
    
class SigStd(Parameter):
    def __init__(self, name='sigStd',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_std

class RefStd(Parameter):
    def __init__(self, name='refStd',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_std
    
class VoltageY(Parameter):
    def __init__(self, name='y', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.ConfocalXY8Object = measurementObject

    def set_raw(self,y):
        galvo.voltage_cdaq1mod2ao1(y)
        print('--------------------------------------------------------')
        print('--------------------------------------------------------')
        print('--------------------------------------------------------')
        print('Set y to ' + str(y) + " V")

        global y_current; y_current=y

class VoltageX(Parameter):
    def __init__(self, name='x', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.ConfocalXY8Object = measurementObject
        self.settleTime = settings['settleTime']
        self.tausArray = self.settings['tausArray']
        self.tausArray2 = self.settings['tausArray2']
        self.MWfreqDictFile = self.settings['MWfreqDictFile']
        self.MWfreqDictFilePlus = self.settings['MWfreqDictFilePlus']
        self.XY8DictFile = self.settings['XY8DictFile']
        self.BThres = self.settings['BThres']
        self.XY8contrastThres = self.settings['XY8contrastThres']
        self.BExt = self.settings['BExt']

        with open(self.MWfreqDictFile, 'rb') as f:
            self.MWfreqDict = pickle.load(f)
        prev = None
        for key, value in self.MWfreqDict.items():
            if value >= 0:
                prev = value
            else:
                if prev is not None:
                    self.MWfreqDict[key] = prev  # Replace with last non-negative value
        
        with open(self.MWfreqDictFilePlus, 'rb') as f:
            self.MWfreqDictPlus = pickle.load(f)
        prev = None
        for key, value in self.MWfreqDictPlus.items():
            if value >= 0:
                prev = value
            else:
                if prev is not None:
                    self.MWfreqDictPlus[key] = prev  # Replace with last non-negative value

        if self.XY8DictFile is not None:
            with open(self.XY8DictFile, 'rb') as f:
                self.XY8Dict = pickle.load(f)

        sig = self.ConfocalXY8Object.sig # this is implemented as a Parameter
        ref = self.ConfocalXY8Object.ref # this is implemented as a Parameter
        sigStd = self.ConfocalXY8Object.sigStd
        refStd = self.ConfocalXY8Object.refStd
        self.loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,sigStd, refStd).then(qctask(sig.close))

    def set_raw(self,x):
        sig = self.ConfocalXY8Object.sig # this is implemented as a Parameter
        self.tausArray = self.settings['tausArray']      
        self.tausArray2 = self.settings['tausArray2']       

        global x_current
        x_current = x
        galvo.voltage_cdaq1mod2ao0(x)
        time.sleep(self.settleTime/1e9)
        print('--------------------------------------------------------')
        print('Set x to ' + str(x) + " V")
        print('Current y = '+ str(y_current) + " V")

        MWfreq =      self.MWfreqDict[(np.round(x_current,2), np.round(y_current,2))]
        MWfreqPlus =  self.MWfreqDictPlus[(np.round(x_current,2), np.round(y_current,2))]
        if self.XY8DictFile is not None:
            XY8contrast = self.XY8Dict.get((np.round(x_current, 2), np.round(y_current, 2)), -1000)
        srs.set_freq(MWfreq)
        print("Set SRS freq to " + str(np.round(MWfreq/1e6,3)) + ' MHz')

        Bstray = ((MWfreqPlus-MWfreq)/1e3) / (2*28) - self.BExt
        global ifOnFlake
        
        if self.XY8DictFile is not None:
            if XY8contrast <= self.XY8contrastThres:
                print("XY8contrast = %.4f. On flake" % XY8contrast)
                ifOnFlake=1
                self.ConfocalXY8Object.tausArray = self.tausArray
                self.loop.sweep_values = sig.sweep(keys=self.tausArray)
            else:
                print("XY8contrast = %.4f. Off flake" % XY8contrast)
                ifOnFlake=0
                self.ConfocalXY8Object.tausArray = self.tausArray2
                self.loop.sweep_values = sig.sweep(keys=self.tausArray2)
        else:
            if Bstray <= self.BThres:
                print("Bstray = %.0f $\mu$T. On flake" % Bstray)
                ifOnFlake=1
                self.ConfocalXY8Object.tausArray = self.tausArray
                self.loop.sweep_values = sig.sweep(keys=self.tausArray)
            else:
                print("Bstray = %.0f $\mu$T. Off flake" % Bstray)
                ifOnFlake=0
                self.ConfocalXY8Object.tausArray = self.tausArray2
                self.loop.sweep_values = sig.sweep(keys=self.tausArray2)

        

