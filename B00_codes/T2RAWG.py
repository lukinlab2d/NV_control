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
from B00_codes.PlotPulseNew import *    
from B00_codes.Confocal import *

class T2RAWG(Instrument):

    def __init__(self, name='T2RAWGObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 #MHz
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
        global pb

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,
                        'AWG': self.AWGParam, 'PB_type': 'USB',
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.tausArray = self.settings['tausArray']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)

        self.SRSnum=self.settings['SRSnum']; self.uwPower = self.settings['uwPower']; self.uwFreq = self.settings['uwFreq']
        self.SDGnum = self.settings['SDGnum']

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
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(self.uwFreq) #Hz
        self.srs.set_RFAmplitude(self.uwPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()
    
        # AWG object
        self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum)
        global AWG; AWG = self.AWG

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

        data = loop.get_data_set(name='T2RAWG')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.T2RAWGObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.T2RAWGObject_ref, name='ref')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        # img = Image.open(dataPlotFile)
        # img.show()

        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        self.AWG.turn_off()
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T2RAWGObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.trackingSettings = self.settings['trackingSettings']
        self.T2RAWGObject = measurementObject
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
                cfcObject = Confocal(settings=self.trackingSettings, laserChannel=self.settings['laserRead_channel'])
                # cfcObject.optimize_xy()
                # time.sleep(1)
                cfcObject.optimize_xz()
                time.sleep(1)
                cfcObject.optimize_xy()
                time.sleep(1)
                cfcObject.close()
                self.timeLastTracking = time.time()
        elif self.trackingSettings['if_tracking'] == 2:
            if time.time() - self.timeLastTracking > self.trackingSettings['time_btwn_trackings']: 
                print()
                cfcObject = Confocal(settings=self.trackingSettings, laserChannel=self.settings['laserRead_channel'])
                x1, y1, z = cfcObject.optimize_xy(direction=1)
                time.sleep(1)
                x2, y2, z = cfcObject.optimize_xy(direction=-1)
                time.sleep(1)
                cfcObject.set_coordinate_fnc((x1+x2)/2, (y1+y2)/2, z)
                cfcObject.close()
                self.timeLastTracking = time.time()
                
        return sig_avg

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))
        
        # Pulse parameters
        num_loops               = self.settings['num_loops']
        laser_init_delay        = self.settings['laser_init_delay'];        laser_init_duration = self.settings['laser_init_duration']
        laser_to_MW_delay       = self.settings['laser_to_MW_delay'];       pi2time        = self.settings['pi2time']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];      read_duration       = self.settings['read_duration']   
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay'];  AWG_output_delay = self.settings['AWG_output_delay']
        AWG_buffer = self.settings['AWG_buffer'];                           #read_offset_from_AWG_delay = self.settings['read_offset_from_AWG_delay']
        MW_duration = int(2*int((AWG_buffer + 2*pi2time + tau_ns + 1)/2))
        
        when_init_end  = laser_init_delay + laser_init_duration
        MW_delay       = when_init_end + laser_to_MW_delay
        when_sigMW_end = MW_delay + AWG_output_delay + MW_duration 
        global MW_del; MW_del = MW_delay + AWG_output_delay

        # read_offset = (laser_to_DAQ_delay-read_offset_from_AWG_delay)
        
        laser_read_signal_delay    = when_sigMW_end
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

        sig_to_ref_wait = laser_read_ref_delay - 2*MW_duration - MW_del
        self.read_duration = read_signal_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")

        global ch1plot; global ch2plot
        ch1plot, ch2plot = AWG.send_T2R_seq(pi_2time=int(pi2time), tau = int(tau_ns), 
                                            buffer=int(AWG_buffer), sig_to_ref_wait=int(sig_to_ref_wait))

        # Make pulse sequence
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,              duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,          duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('AWG', MW_delay,                     duration=20)]
        pulse_sequence += [spc.Pulse('Counter',  read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
           
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
                fig = plotPulseObject.makePulsePlot(ch1plot, ch2plot, MW_del)
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T2RAWGObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
    
    # def turn_on_mid_sweep(self):
    #     pb = TurnOnLaser.turnOnLaser(channel=laserChannel)

    def turn_on_at_end(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        # pb.turn_on_infinite(channels=channels)


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