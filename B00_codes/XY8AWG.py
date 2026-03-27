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

class XY8AWG(Instrument):

    def __init__(self, name='XY8AWGObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
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
        self.SDGnum = self.settings['SDGnum']; self.srate=self.settings['srate']

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
            name = "sigStd",
            parameter_class = SigStd,
        )
        self.add_parameter(
            name = "refStd",
            parameter_class = RefStd,
        )
        self.savedPulseSequencePlots = {}

        # SRS object
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(self.uwFreq) #Hz
        self.srs.set_RFAmplitude(self.uwPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()
    
        # AWG object
        self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum, srate=self.srate)
        global AWG; AWG = self.AWG

    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigStd = self.sigStd
        refStd = self.refStd

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig, ref, sigStd, refStd,
                                            qctask(sig.plotPulseSequences),
                                            ).then(qctask(sig.turn_on_at_end))

        data = loop.get_data_set(name='XY8AWG')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.XY8AWGObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.XY8AWGObject_ref, name='ref')

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
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/XY8AWGObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.trackingSettings = self.settings['trackingSettings']
        self.XY8AWGObject = measurementObject
        self.loopCounter = 0
        self.timeLastTracking = time.time()
        self.tausArray = self.settings['tausArray']
        self.srate = self.settings['srate']

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
        global sig_std;  sig_std = np.std(sig)
        global ref_std;  ref_std = np.std(ref)

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
        rest_after_first_pulse=0
        num_loops               = self.settings['num_loops']
        mode                    = self.settings['mode'];                   numxy8              = self.settings['numxy8']
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration = self.settings['laser_init_duration']
        laser_to_AWG_delay      = self.settings['laser_to_AWG_delay'];     tau                 = self.settings['tau']
        pi2time                 = self.settings['pi2time'];                pitime              = self.settings['pitime']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];     read_duration       = self.settings['read_duration']   
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']; AWG_output_delay    = self.settings['AWG_output_delay']
        AWG_buffer              = self.settings['AWG_buffer'];             sweepWhich          = self.settings['sweepWhich']
        
        if sweepWhich=='N':
            numxy8 = tau_ns
            tau_ns = tau
        
        # sig_to_ref_wait = (read_duration + laser_to_DAQ_delay + DAQ_to_laser_off_delay) + laser_to_AWG_delay
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
        
        laser_read_signal_delay    = rounded_when_sigMW_end
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay;   read_signal_duration = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        laser_read_signal_duration = when_read_signal_end + DAQ_to_laser_off_delay - laser_read_signal_delay
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        
        MW_del_ref      = when_laser_read_signal_end + MW_del
        # print(when_laser_read_signal_end)
        # print(MW_del)
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
        
        if self.srate is not None:
            if self.loopCounter==0: sleepTime = 10
            else: sleepTime = 0.75
        else:
            sleepTime = 0

        global ch1plot; global ch2plot
        ch1plot, ch2plot = AWG.send_XY8_seq(mode = mode, numxy8=numxy8, pi_2time=int(pi2time),rest_after_first_pulse=rest_after_first_pulse,
                                            pitime = int(pitime), tau = int(tau_ns),AWGnum=2,
                                              buffer=int(AWG_buffer), sig_to_ref_wait=int(sig_to_ref_wait),
                                              sleepTime=sleepTime)

        # Make pulse sequence
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,              duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay,          duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,          duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('AWG',          AWG_trig_delay,                     duration=50)]
        pulse_sequence += [spc.Pulse('Counter',  read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
        
        # print(laser_read_signal_delay)
        # print(when_sigMW_end)
        # print(MW_duration)

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

        if sweepWhich == 'tau':
            print('Set tau to ' + str(tau_ns) + " ns")
        elif sweepWhich == 'N':
            print('Set N to ' + str(numxy8))
        if not self.settings['ifPlotPulse']: self.loopCounter += 1
    
    def plotPulseSequences(self):
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False)
                fig = plotPulseObject.makePulsePlot(ch1plot, ch2plot, MW_del)
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.XY8AWGObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
    
    # def turn_on_mid_sweep(self):
    #     pb = TurnOnLaser.turnOnLaser(channel=laserChannel)

    def turn_on_at_end(self):
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