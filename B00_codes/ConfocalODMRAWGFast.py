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
    

class ConfocalODMRAWGFast(Instrument):

    def __init__(self, name='ConfocalODMRAWGFastObject', settings=None, ifPlotPulse=True, **kwargs) -> None:

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
        self.galvo = DAQAnalogOutputs("ConfocalODMRAWGFastDAQ", galvo_card_name, galvo_ao_channels)
        global galvo; galvo = self.galvo
    
    def runScan(self):
        # sig = self.sig # this is implemented as a Parameter
        # ref = self.ref # this is implemented as a Parameter
        y = self.y; x = self.x

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        # loop = Loop(
        #     sig.sweep(keys=self.freqsArray),
        #     delay = 0,
        #     sleepTimeAfterFinishing=0).each(sig,ref).then(qctask(sig.close))

        # Loop through x and y
        loopx = Loop(
            x.sweep(keys=self.xArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(x.loop)
        loopy = Loop(
            y.sweep(keys=self.yArray),
            delay=0,
            sleepTimeAfterFinishing=0).each(loopx)

        data = loopy.get_data_set(name='ConfocalODMRAWGFast')
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
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ConfocalODMRAWGFastObject_sig_set.dat'


class Signal(Parameter):
    def __init__(self, name='sig', settings=None, measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.loopCounter = 0
        self.ConfocalODMRAWGFastObject = measurementObject

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
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)
        global sig_avg_over_ref_avg; sig_avg_over_ref_avg = sig_avg/ref_avg

        # time.sleep(0.4)

        return sig_avg
    
    def set_raw(self, value):
        srs.set_freq(value)
        print("Loop " + str(self.loopCounter))
        # print("Set SRS freq to " + str(np.round(value/1e6,3)) + ' MHz') # set the SRS frequency in reality

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

class VoltageY(Parameter):
    def __init__(self, name='y', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.ConfocalODMRAWGFastObject = measurementObject

    def set_raw(self,y):
        galvo.voltage_cdaq1mod2ao1(y)
        print('Set y to ' + str(y) + " V")

        global y_current; y_current=y

class VoltageX(Parameter):
    def __init__(self, name='x', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.ConfocalODMRAWGFastObject = measurementObject
        self.settleTime = settings['settleTime']
        self.freqsArray = self.settings['freqsArray']

        sig = self.ConfocalODMRAWGFastObject.sig # this is implemented as a Parameter
        ref = self.ConfocalODMRAWGFastObject.ref # this is implemented as a Parameter
        self.loop = Loop(
            sig.sweep(keys=self.freqsArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref).then(qctask(sig.close))

    def set_raw(self,x):
        sig = self.ConfocalODMRAWGFastObject.sig # this is implemented as a Parameter
        ifDifferential = self.settings['ifDifferential']

        if ifDifferential:
            if x < 0:
                self.freqsArray = self.settings['freqsArray']
                pitime = self.settings['pitime']; pi_incr_factor = self.settings['pi_incr_factor']
                padding = self.settings['padding']
                
            else:
                x = -x
                self.freqsArray = self.settings['freqsArray2']
                pitime = self.settings['pitime2']; pi_incr_factor = self.settings['pi_incr_factor2']
                padding = self.settings['padding'] - (self.settings['pitime2']-self.settings['pitime'])
        else:
            self.freqsArray = self.settings['freqsArray']
            pitime = self.settings['pitime']; pi_incr_factor = self.settings['pi_incr_factor']
            padding = self.settings['padding']
           
        self.ConfocalODMRAWGFastObject.freqsArray = self.freqsArray
        self.loop.sweep_values = sig.sweep(keys=self.freqsArray)

        galvo.voltage_cdaq1mod2ao0(x)
        time.sleep(self.settleTime/1e9)
        print('--------------------------------------------------------')
        print('Set x to ' + str(x) + " V")
        print('Current y = '+ str(y_current) + " V")

        xref_pitime = self.settings['xref_pitime']#; pi_incr_factor = self.settings['pi_incr_factor']
        pi_increment = int((x-(xref_pitime))*pi_incr_factor)
        ## if x <= -1.15:
        ##     # pi_increment = int(((-2)-(xref_pitime))*pi_incr_factor)
        ##     pi_increment = self.settings['pitime2']-self.settings['pitime']
        # if y_current < -0.2:
        #     if x <= 0.0934*y_current + 0.1392:
        #         print('On BN')
        #         pi_increment = self.settings['pitime2']-self.settings['pitime']
        #     else:
        #         print('On FGT')
        # else:
        #     if x <= -0.5837*y_current + 0.0126:
        #         print('On BN')
        #         pi_increment = self.settings['pitime2']-self.settings['pitime']
        #     else:
        #         print('On FGT')
        pi_increment = int(2*int(pi_increment/2))
        # pi_increment = int(10*int(pi_increment/10))
        print('pi_increment = ' + str(pi_increment))

        global num_loops; global read_duration

        # Pulse parameters
        AWGbuffer               = self.settings['AWGbuffer'];               pitime              = pitime + pi_increment
        num_loops               = self.settings['num_loops'];               DAQ_error_factor    = self.settings['DAQ_error_factor']
        num_loops               = int(num_loops/DAQ_error_factor);          AWG_output_delay    = self.settings['AWG_output_delay']
        MW_duration             = int(2*int((2*AWGbuffer + pitime + 1)/2)); MW_to_DAQ_delay     = self.settings['MW_to_DAQ_delay']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];      read_duration       = self.settings['read_duration']   
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay'];  wait_btwn_sig_ref   = DAQ_to_laser_off_delay
        padding                 = padding - pi_increment;                   padding_green1      = self.settings['padding_green1']

        # if pitime < 300 and 950 <= padding and padding <= 1200: padding = 1250
        # elif pitime > 500 and 780 < padding and padding < 850: padding = 850
        print('True padding = ', padding)

        # new delays
        read_signal_duration = read_duration;         read_ref_duration = read_duration
        serious_duration     = 2*read_duration + wait_btwn_sig_ref

        MW_delay = 0; 
        global MW_del; MW_del = MW_delay

        when_pulse_end = MW_delay + MW_duration
        read_signal_delay = when_pulse_end + MW_to_DAQ_delay
        read_ref_delay = read_signal_delay + read_signal_duration + wait_btwn_sig_ref

        total_time = read_signal_delay + serious_duration 

        laser_read_part1_duration = total_time - laser_to_DAQ_delay
        laser_read_part1_delay    = 0
        laser_read_part2_duration = serious_duration - laser_read_part1_duration
        laser_read_part2_delay    = total_time - laser_read_part2_duration + padding
        when_read_part2_ends      = laser_read_part2_delay + laser_read_part2_duration

        if self.settings['AWG_delay'] > (when_read_part2_ends + 20):
            AWG_delay = when_read_part2_ends - (self.settings['AWG_delay']-when_read_part2_ends)
        else:
            AWG_delay = when_read_part2_ends - self.settings['AWG_delay']

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")    

        # Make pulse sequence (per each freq)
        global pulse_sequence; pulse_sequence = []
        pulse_sequence += [spc.Pulse('LaserRead',  laser_read_part1_delay, duration=int(laser_read_part1_duration + padding_green1))] # times are in ns
        pulse_sequence += [spc.Pulse('AWG',        AWG_delay,              duration=20)]
        pulse_sequence += [spc.Pulse('LaserRead',  laser_read_part2_delay, duration=int(laser_read_part2_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',    read_signal_delay,      duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',    read_ref_delay,         duration=int(read_ref_duration))] # times are in ns
        self.pulse_sequence = pulse_sequence

        global num_reads_per_iter; num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        
        global num_reads; num_reads = int(num_loops * num_reads_per_iter * DAQ_error_factor) # the results look like this: sig-ref-sig-ref-...-sig-ref for num_loops times

        ch1, ch2 = AWG.send_ODMR_seq(pulse_width=int(pitime), buffer=int(AWGbuffer))
        global ch1plot; global ch2plot 
        ch1plot = ch1; ch2plot = ch2 

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

