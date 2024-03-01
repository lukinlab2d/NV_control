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

class T1SCC(Instrument):

    def __init__(self, name='T1SCCObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz, of the Pulse Blaster
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.LaserIonParam =    {'delay_time': 2, 'channel':settings['laserIon_channel']}
        self.LaserShelveParam = {'delay_time': 2, 'channel':settings['laserShelve_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':1}
        self.MWQParam =         {'delay_time': 2, 'channel':0}
        self.MWswitchParam =    {'delay_time': 2, 'channel':2}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
        global laserReadChannel; laserReadChannel = self.LaserReadParam['channel']
        global pb

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam, 'LaserIon': self.LaserIonParam, 'LaserShelve': self.LaserShelveParam,
                        'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.tausArray = self.settings['tausArray']
        self.uwPower = self.settings['uwPower']; self.uwFreq = self.settings['uwFreq']

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
        self.add_parameter(
            name = "sigOverRef",
            parameter_class = SigOverRef,
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
        self.srs = SRS()
        self.srs.set_freq(self.uwFreq) #Hz
        self.srs.set_RFAmplitude(self.uwPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigOverRef = self.sigOverRef
        sigFullData = self.sigFullData
        refFullData = self.refFullData

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            ifIndividualCount = True, 
            numOfPBLoops = self.settings['num_loops'],
            sleepTimeAfterFinishing=0).each(
                                            sig, ref,
                                            sigFullData, refFullData,
                                            qctask(sig.plotPulseSequences),
                                            ).then(qctask(sig.turn_on_at_end))

        data = loop.get_data_set(name='T1SCC')
        data.add_metadata(self.settings)
        self.data = data
        
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T1SCCObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.trackingSettings = self.settings['trackingSettings']
        self.T1SCCObject = measurementObject
        self.loopCounter = 0
        self.timeLastTracking = time.time()
        self.tausArray = self.settings['tausArray']

        self.ifMeaningfulRef = self.settings['ifMeaningfulRef']

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))

        self.ctrtask.stop(); self.ctrtask.close()

        # rate = xLineData/(self.read_duration/1e9)/1e3
        rate = xLineData
        global sig_data; sig_data = rate[::self.num_reads_per_iter]
        global ref_data; 
        if self.ifMeaningfulRef:
            ref_data = rate[1::self.num_reads_per_iter]
        else:
            ref_data = sig_data
        global sig_avg;  sig_avg = np.average(sig_data)
        global ref_avg;  ref_avg = np.average(ref_data)
        global sig_avg_over_ref_avg; sig_avg_over_ref_avg = sig_avg/ref_avg

        # NV tracking
        if self.trackingSettings['if_tracking'] == 1:
            # if np.mod(self.loopCounter, self.trackingSettings['tracking_period']) == self.trackingSettings['tracking_period']-1:
            if time.time() - self.timeLastTracking > self.trackingSettings['time_btwn_trackings']:    
                print()
                cfcObject = Confocal(settings=self.trackingSettings, laserChannel=self.settings['laserTrack_channel'])
                cfcObject.optimize_xz()
                time.sleep(1)
                cfcObject.optimize_xy()
                time.sleep(1)
                cfcObject.close()
                self.timeLastTracking = time.time()

        return sig_avg

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster

        print("Loop " + str(self.loopCounter))
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];             
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];     MWI_duration        = self.settings['pi_time']
        MWI_to_shelve_delay     = tau_ns;                                  shelve_duration     = self.settings['shelve_duration']  
        shelve_to_ion_delay     = self.settings['shelve_to_ion_delay'];    ion_duration        = self.settings['ion_duration']
        ion_to_laserRead_delay  = self.settings['ion_to_laserRead_delay']
        laserRead_to_DAQ_delay  = self.settings['laserRead_to_DAQ_delay']; DAQ_duration        = self.settings['DAQ_duration']
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']

        # Make pulse sequence
        pulse_sequence = [];                                                      when_init_end             = laser_init_delay + laser_init_duration
        MWI_delay              = when_init_end          + laser_to_MWI_delay;     when_MWI_end              = MWI_delay        + MWI_duration
        shelve_delay           = when_MWI_end           + MWI_to_shelve_delay;    when_shelve_end           = shelve_delay     + shelve_duration
        ion_delay              = when_shelve_end        + shelve_to_ion_delay;    when_ion_end              = ion_delay        + ion_duration 
        laserRead_signal_delay = when_ion_end           + ion_to_laserRead_delay
        DAQ_signal_delay       = laserRead_signal_delay + laserRead_to_DAQ_delay; when_DAQ_signal_end       = DAQ_signal_delay + DAQ_duration
        zzz                    = 1;                                               when_laserRead_signal_end = when_DAQ_signal_end + DAQ_to_laser_off_delay
        
        laserRead_signal_duration = when_laserRead_signal_end - laserRead_signal_delay
        total_sig_duration        = when_laserRead_signal_end + laserRead_to_DAQ_delay
        
        if not laser_init_delay == 0:
            pulse_sequence     += [spc.Pulse('LaserInit',    laser_init_delay,       duration=int(laser_init_duration))] # times are in ns
        pulse_sequence         += [spc.Pulse('MWswitch',     MWI_delay,              duration=int(MWI_duration))]
        if not shelve_duration == 0:
            pulse_sequence     += [spc.Pulse('LaserShelve',  shelve_delay,           duration=int(shelve_duration))]
        pulse_sequence         += [spc.Pulse('LaserIon',     ion_delay,              duration=int(ion_duration))]
        pulse_sequence         += [spc.Pulse('LaserRead',    laserRead_signal_delay, duration=int(laserRead_signal_duration))] 
        pulse_sequence         += [spc.Pulse('Counter',      DAQ_signal_delay,       duration=int(DAQ_duration))] 
        if self.ifMeaningfulRef:
            pulse_sequence     += [spc.Pulse('LaserInit',    laser_init_delay       + total_sig_duration, duration=int(laser_init_duration))] # times are in ns
            if not shelve_duration == 0:
                pulse_sequence += [spc.Pulse('LaserShelve',  shelve_delay           + total_sig_duration, duration=int(shelve_duration))]
            pulse_sequence     += [spc.Pulse('LaserIon',     ion_delay              + total_sig_duration, duration=int(ion_duration))]
            pulse_sequence     += [spc.Pulse('LaserRead',    laserRead_signal_delay + total_sig_duration, duration=int(laserRead_signal_duration))] 
            pulse_sequence     += [spc.Pulse('Counter',      DAQ_signal_delay       + total_sig_duration, duration=int(DAQ_duration))] 
        
        self.read_duration = DAQ_duration
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
            samps_per_chan = 2*self.num_reads # x2 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        print('Set tau to ' + str(tau_ns/1e6) + " ms")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1
    
    def plotPulseSequences(self):
        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']
        self.ionColor     = self.settings['LaserIon']['channel']
        self.shelveColor  = self.settings['LaserShelve']['channel']
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, 
                                            ifShown=True, ifSave=False, readColor=self.readColor, 
                                            initColor=self.initColor, ionColor=self.ionColor, shelveColor=self.shelveColor)
                fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T1SCCObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
    
    def turn_on_at_end(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        if laserInitChannel == laserReadChannel:
            channels = np.linspace(laserInitChannel,laserInitChannel,1)
        else:
            channels = np.array((laserInitChannel,laserReadChannel))
        pb.turn_on_infinite(channels=channels)


class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg

class SigOverRef(Parameter):
    def __init__(self, name='sigOverRef',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_avg_over_ref_avg
    
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