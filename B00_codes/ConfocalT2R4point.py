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
    

class ConfocalT2R4point(Instrument):

    def __init__(self, name='ConfocalT2R4pointObject', settings=None, ifPlotPulse=True, **kwargs) -> None:

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
                          'min_pulse_dur': int(5*1e3/self.clock_speed)}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']
        self.ifPlotPulse  = ifPlotPulse

        # List of frequencies and power
        # self.freqsArray = self.settings['freqsArray']
        self.SRSnum=self.settings['SRSnum']; MWPower = self.settings['MWPower']
        self.SDGnum=self.settings['SDGnum']

        self.xArray = self.settings['xArray']; self.yArray = self.settings['yArray']
        
        # SRS object
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(3e9) #Hz
        self.srs.set_RFAmplitude(MWPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        # AWG object
        self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum)
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
            name = "ref2",
            parameter_class = Reference2,
        )
        self.add_parameter(
            name = "ref3",
            parameter_class = Reference3,
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
        self.galvo = DAQAnalogOutputs("ConfocalT2R4pointDAQ", galvo_card_name, galvo_ao_channels)
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

        data = loopy.get_data_set(name='ConfocalT2R4point')
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
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ConfocalT2R4pointObject_sig_set.dat'


class Signal(Parameter):
    def __init__(self, name='sig', settings=None, measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.loopCounter = 0
        self.ConfocalT2R4pointObject = measurementObject

    def get_raw(self):
        self.loopCounter += 1        
        ctrtask.start()

        pb.start_pulse_seq()
        pb.wait()
        pb.stop_pulse_seq_without_closing()

        data = np.array(ctrtask.read(num_reads, timeout=10))

        ctrtask.stop()#; self.ctrtask.close()
        
        rate = data/(read_duration/1e9)/1e3
        sig = rate[::num_reads_per_iter]
        ref = rate[1::num_reads_per_iter]
        ref2 = rate[2::num_reads_per_iter]
        ref3 = rate[3::num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)
        global ref_avg2; ref_avg2 = np.average(ref2)
        global ref_avg3; ref_avg3 = np.average(ref3)

        return sig_avg
    
    def set_raw(self, value):
        srs.set_freq(value)
        print("-------------------------------------")
        print("Loop " + str(self.loopCounter))
        print("Set SRS freq to " + str(np.round(value/1e6,3)) + ' MHz') # set the SRS frequency in reality

        pitime = self.settings['pitime']; pi_incr_factor = self.settings['pi_incr_factor']
        pitime2 = self.settings['pitime2']; pi_incr_factor2 = self.settings['pi_incr_factor2']

        xref_pitime = self.settings['xref_pitime']
        pi_increment = int((current_x-(xref_pitime))*pi_incr_factor)
        pi_increment = int(2*int(pi_increment/2))
        print('pi_increment = ' + str(pi_increment))

        pi_increment2 = int((current_x-(xref_pitime))*pi_incr_factor2)
        pi_increment2 = int(2*int(pi_increment2/2))
        print('pi_increment2 = ' + str(pi_increment2))

        global num_loops; global read_duration

        # Pulse parameters
        tau_ns                  = self.settings['tau'];                     
        pitime                  = pitime + pi_increment;                    pitime2                 = pitime2 + pi_increment2
        num_loops               = self.settings['num_loops'];               laser_init_delay        = self.settings['laser_init_delay']
        laser_init_duration     = self.settings['laser_init_duration'];     laser_to_AWG_delay      = self.settings['laser_to_AWG_delay']
        pi2time                 = int(pitime/2);                            pi2time2                = int(pitime2/2)
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];      read_duration           = self.settings['read_duration']   
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay'];  AWG_output_delay        = self.settings['AWG_output_delay']
        AWGbuffer               = self.settings['AWGbuffer'];               read_offset_from_AWG_delay = self.settings['read_offset_from_AWG_delay']
        DAQ_error_factor        = self.settings['DAQ_error_factor'];        num_loops               = int(num_loops/DAQ_error_factor)
        
        if np.mod(self.loopCounter,2)==1:
            pi2time = pi2time2
        print('Pi half = ' + str(pi2time) + ' ns')

        MW_duration = int(2*int((AWGbuffer + 2*pi2time + tau_ns + 1)/2))

        when_init_end = laser_init_delay + laser_init_duration
        MW_delay     = when_init_end + laser_to_AWG_delay;     when_sigMW_end = MW_delay + AWG_output_delay + MW_duration 
        global MW_del; MW_del = MW_delay+AWG_output_delay

        read_offset = (laser_to_DAQ_delay-read_offset_from_AWG_delay)
        
        laser_read_signal_delay    = when_sigMW_end - read_offset
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        laser_read_signal_duration = when_read_signal_end + DAQ_to_laser_off_delay - laser_read_signal_delay
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        
        laser_read_ref_delay = when_laser_read_signal_end + laser_read_signal_delay
        read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay
        read_ref_duration    = read_duration
        when_read_ref_end    = read_ref_delay + read_ref_duration
        laser_read_ref_duration = when_read_ref_end + DAQ_to_laser_off_delay - laser_read_ref_delay
        when_laser_read_ref_end = laser_read_ref_delay + laser_read_ref_duration

        sig_to_ref_wait = laser_read_signal_duration + (laser_read_signal_delay-2*pi2time-tau_ns)

        laser_read_ref_delay2 = when_laser_read_ref_end + laser_read_signal_delay
        read_ref_delay2       = laser_read_ref_delay2 + laser_to_DAQ_delay
        read_ref_duration2    = read_duration
        when_read_ref_end2    = read_ref_delay2 + read_ref_duration2
        laser_read_ref_duration2 = when_read_ref_end2 + DAQ_to_laser_off_delay - laser_read_ref_delay2
        when_laser_read_ref_end2 = laser_read_ref_delay2 + laser_read_ref_duration2

        laser_read_ref_delay3 = when_laser_read_ref_end2 + laser_read_signal_delay
        read_ref_delay3       = laser_read_ref_delay3 + laser_to_DAQ_delay
        read_ref_duration3    = read_duration
        when_read_ref_end3    = read_ref_delay3 + read_ref_duration3
        laser_read_ref_duration3 = when_read_ref_end3 + DAQ_to_laser_off_delay - laser_read_ref_delay3
        when_laser_read_ref_end3 = laser_read_ref_delay3 + laser_read_ref_duration3

        self.read_duration = read_signal_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")

        global ch1plot; global ch2plot
        ch1plot, ch2plot = AWG.send_T2R4point_seq(pi_2time=int(pi2time), tau = int(tau_ns), buffer=int(AWGbuffer), sig_to_ref_wait=int(sig_to_ref_wait))

        # Make pulse sequence
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,              duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,          duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay2,         duration=int(laser_read_ref_duration2))]
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay3,         duration=int(laser_read_ref_duration3))]
        pulse_sequence += [spc.Pulse('AWG', MW_delay,                     duration=20)]
        pulse_sequence += [spc.Pulse('Counter',  read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',      read_ref_delay2,               duration=int(read_ref_duration2))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',      read_ref_delay3,               duration=int(read_ref_duration3))] # times are in ns
           
        self.pulse_sequence = pulse_sequence

        global num_reads_per_iter; num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        global num_reads; num_reads = int(num_loops * num_reads_per_iter * DAQ_error_factor) # the results look like this: sig-ref-sig-ref-...-sig-ref for num_loops times
        print('DAQ_error_factor = ' + str(DAQ_error_factor))

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
            samps_per_chan = int(2*num_reads) # x5 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        
        # Pulse Blaster object inittialized here to avoid long delay between initialization and first run
        global pb
        pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        pb.program_pb(pulse_sequence, num_loops=num_loops)

        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        self.pb.program_pb(pulse_sequence, num_loops=num_loops)

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

class Reference2(Parameter):
    def __init__(self, name='ref2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg2

class Reference3(Parameter):
    def __init__(self, name='ref3',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg3

class VoltageY(Parameter):
    def __init__(self, name='y', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.ConfocalT2R4pointObject = measurementObject

    def set_raw(self,y):
        galvo.voltage_cdaq1mod2ao1(y)
        print('Set y to ' + str(y) + " V")

        global y_current; y_current=y

class VoltageX(Parameter):
    def __init__(self, name='x', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.ConfocalT2R4pointObject = measurementObject
        self.settleTime = settings['settleTime']
        self.freqsArray = self.settings['freqsArray']

        sig = self.ConfocalT2R4pointObject.sig # this is implemented as a Parameter
        ref = self.ConfocalT2R4pointObject.ref # this is implemented as a Parameter
        ref2 = self.ConfocalT2R4pointObject.ref2
        ref3 = self.ConfocalT2R4pointObject.ref3
        self.loop = Loop(
            sig.sweep(keys=self.freqsArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig, ref, ref2, ref3,).then(qctask(sig.close))

    def set_raw(self,x):
        sig = self.ConfocalT2R4pointObject.sig # this is implemented as a Parameter
        self.freqsArray = self.settings['freqsArray']       
        self.ConfocalT2R4pointObject.freqsArray = self.freqsArray
        self.loop.sweep_values = sig.sweep(keys=self.freqsArray)

        global current_x
        current_x = x
        galvo.voltage_cdaq1mod2ao0(x)
        time.sleep(self.settleTime/1e9)
        print('Set x to ' + str(x) + " V")
        print('Current y = '+ str(y_current) + " V")

        

        
