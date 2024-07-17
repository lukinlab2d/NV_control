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

class T1(Instrument):

    def __init__(self, name='T1Object', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserParam =       {'delay_time': 2, 'channel':3}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':1}
        self.MWQParam =         {'delay_time': 2, 'channel':1}
        self.MWswitchParam =    {'delay_time': 2, 'channel':2}
        global laserChannel; laserChannel = self.LaserParam['channel']

        settings_extra = {'clock_speed': self.clock_speed, 'Laser': self.LaserParam, 'Counter': self.CounterParam, 
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
        self.savedPulseSequencePlots = {}

        # SRS object
        self.srs = SRS(SRSnum=2)
        self.srs.set_freq(self.uwFreq) #Hz
        self.srs.set_RFAmplitude(self.uwPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigOverRef = self.sigOverRef

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig, ref, sigOverRef,
                                            qctask(sig.plotPulseSequences),
                                            ).then(qctask(sig.turn_on_at_end))

        data = loop.get_data_set(name='T1')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.T1Object_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.T1Object_ref, name='ref')
        # plot = QtPlot(
        #     data.T1Object_sigOverRef, # this is implemented as a Parameter
        #     figsize = (1200, 600),
        #     interval = 1,
        #     name = 'sig/ref'
        #     )

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

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T1Object_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.trackingSettings = self.settings['trackingSettings']
        self.T1Object = measurementObject
        self.loopCounter = 0
        self.timeLastTracking = time.time()
        self.tausArray = self.settings['tausArray']

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
        num_loops               = self.settings['num_loops'];              ifShimon                = self.settings['ifShimon']
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration     = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];     MWI_duration            = self.settings['pi_time']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];     read_duration           = self.settings['read_duration']   
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']; sig_to_ref_delay_Shimon = self.settings['sig_to_ref_delay_Shimon']

        if True:
            MWI_to_switch_delay  = self.settings['MWI_to_switch_delay']

        # Make pulse sequence
        pulse_sequence = []

        if ifShimon == 0:
            when_init_end   = laser_init_delay+laser_init_duration
            MWI_delay = when_init_end+laser_to_MWI_delay;                 when_pulse_end = MWI_delay+MWI_duration
            
            laser_read_signal_delay    = when_pulse_end + tau_ns
            read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay;   read_signal_duration = read_duration
            when_read_signal_end       = read_signal_delay + read_signal_duration
            laser_read_signal_duration = when_read_signal_end + DAQ_to_laser_off_delay - laser_read_signal_delay; when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
            
            laser_read_ref_delay = when_laser_read_signal_end + laser_to_MWI_delay + MWI_duration + tau_ns
            read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay 
            read_ref_duration    = read_duration; when_read_ref_end = read_ref_delay + read_ref_duration
            laser_read_ref_duration = when_read_ref_end + DAQ_to_laser_off_delay - laser_read_ref_delay

            if not laser_init_delay == 0:
                pulse_sequence += [spc.Pulse('Laser',laser_init_delay,              duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('Laser',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('Laser',    laser_read_ref_delay,          duration=int(laser_read_ref_duration))]
            pulse_sequence += [spc.Pulse('MWswitch', MWI_delay,                     duration=int(MWI_duration))]
            pulse_sequence += [spc.Pulse('Counter',  read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
        elif ifShimon == 1:
            when_init_end           = laser_init_delay + laser_init_duration
            laser_read_signal_delay = when_init_end + tau_ns
            read_signal_delay       = laser_read_signal_delay + laser_to_DAQ_delay;   read_signal_duration = read_duration
            when_read_signal_end    = read_signal_delay + read_signal_duration

            read_ref_delay          = when_read_signal_end + sig_to_ref_delay_Shimon;    read_ref_duration = read_duration
            when_read_ref_end       = read_ref_delay + read_ref_duration

            laser_read_signal_duration = when_read_ref_end + DAQ_to_laser_off_delay - laser_read_signal_delay

            pulse_sequence += [spc.Pulse('Laser',    laser_init_delay,              duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('Laser',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('Counter',  read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
        
        self.read_duration = read_signal_duration
        
        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")
        
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
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False)
                fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T1Object.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1

    def turn_on_at_end(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserChannel,laserChannel,1)
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