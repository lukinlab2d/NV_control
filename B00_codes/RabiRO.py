"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
"""
from  qcodes.actions import Task as qctask
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot
import numpy as np
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS
from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity
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
from B00_codes.ScanROFreq import *  

from qcodes.dataset import do1d, do2d, dond, LinSweep, LogSweep, ArraySweep
from qcodes.utils.dataset.doNd import plot
from qcodes.dataset.sqlite.database import initialise_or_create_database_at
from qcodes.dataset.experiment_container import load_or_create_experiment
from qcodes.tests.instrument_mocks import DummyInstrument, DummyInstrumentWithMeasurement
from qcodes.dataset.measurements import Measurement
from qcodes.dataset.plotting import plot_dataset

def ns2cycles(time, samp_rate=1e7):
        return int(time/1e9*samp_rate)
    

class RabiRO(Instrument):

    def __init__(self, name='RabiROObject', settings=None, ifPlotPulse=True, **kwargs) -> None:

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
    
        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,
                          'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                          'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        # Vpiezos, MW power, and MW frequency
        self.tausArray = self.settings['tausArray']
        self.SRSnum = self.settings['SRSnum'];      MWPower = self.settings['MWPower']; MWFreq = self.settings['MWFreq']
        
        self.velNum = self.settings['velNum']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
        self.vel_vpz_start = self.settings['vel_vpz_start']; self.vel_vpz_step = self.settings['vel_vpz_step']
        self.vel_vpz_step_time = self.settings['vel_vpz_step_time']; self.vel_vpz_target = self.settings['vel_vpz_target']
        self.ifScanVpz = self.settings['ifScanVpz']

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
        self.srs.set_freq(MWFreq) #Hz
        self.srs.set_RFAmplitude(MWPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        # Velocity object
        self.vel = Velocity(velNum=self.velNum)
        self.vel.set_track()
        time.sleep(0.5)
        self.vel.set_wvl(vel_wvl)
        time.sleep(1)
        self.vel.set_ready()
        self.vel.set_current(vel_current)
        self.vel.set_vpiezo(1.1)
        self.vel.waitUntilComplete()
        self.vel.set_ready()
        time.sleep(2)
        if self.ifScanVpz:
            self.vel.set_vpiezo(self.vel_vpz_start)
            self.vel.waitUntilComplete()
            self.vel.set_ready()
            time.sleep(1)

            self.set_vpz()
            time.sleep(3)
        else:
            self.vel.set_vpiezo(self.vel_vpz_target)
            self.vel.waitUntilComplete()
            self.vel.set_ready()
            time.sleep(1)

        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        global srs; srs = self.srs
        global vel; vel = self.vel
    
    def set_vpz(self):
        nstep = int((self.vel_vpz_target-self.vel_vpz_start)/self.vel_vpz_step + 1)
        for vpz in np.linspace(self.vel_vpz_start, self.vel_vpz_target, nstep):
            self.vel.set_vpiezo(vpz)
            self.vel.waitUntilComplete()
            self.vel.set_ready()
            time.sleep(self.vel_vpz_step_time)
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigOverRef = self.sigOverRef

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,sigOverRef,
                                            qctask(sig.plotPulseSequences),
                                            ).then(qctask(sig.sweepVelToEnd))

        data = loop.get_data_set(name='RabiRO')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.RabiROObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.RabiROObject_ref, name='ref')

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
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/RabiROObject_sig_set.dat'

    
class Signal(Parameter):
    def __init__(self, name='sig', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.loopCounter = 0
        self.settings = settings
        self.tausArray = self.settings['tausArray']
        self.RabiROObject = measurementObject
        self.ROtrackingSettings = self.settings['ROtrackingSettings']
        self.timeLastROTracking = time.time()

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        self.vel_vpz_start = self.settings['vel_vpz_start']; self.vel_vpz_step = self.settings['vel_vpz_step']; self.vel_vpz_end = self.settings['vel_vpz_end']
        self.vel_vpz_step_time = self.settings['vel_vpz_step_time']; self.vel_vpz_target = self.settings['vel_vpz_target']

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

        # Line tracking for RO
        if self.settings['laserRead_channel'] == 5:
            if self.ROtrackingSettings['if_tracking'] == 1:
                threshold = self.ROtrackingSettings['threshold']
                if ref_avg < threshold: 
                    print()
                    print('-----------------Start line tracking---------------------------')
                    ScanROFreqObject = ScanROFreq(settings=self.ROtrackingSettings, ifPlotPulse=0)
                    vpz = ScanROFreqObject.runScanInPulseSequence()
                    vel.set_vpiezo(vpz)
                    vel.set_ready()
                    print("Set Vpiezo to " + str(np.round(vpz,1)) + ' %') # set the Velocity's piezo voltage
                    print('-----------------End line tracking---------------------------')
                    print()
                    self.timeLastROTracking = time.time()

        return sig_avg
    
    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))

        # Pulse parameters
        num_loops               = self.settings['num_loops']
        laser_init_delay        = self.settings['laser_init_delay'];    laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];  MWI_duration        = tau_ns
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];  read_duration       = self.settings['read_duration']   
        read_laser_duration     = self.settings['read_laser_duration']; MW_to_read_delay    = self.settings['MW_to_read_delay']
        
        when_init_end              = laser_init_delay + laser_init_duration
        MWI_delay                  = when_init_end    + laser_to_MWI_delay;  when_pulse_end = MWI_delay + MWI_duration
        
        laser_read_signal_delay    = when_pulse_end   + MW_to_read_delay;           laser_read_signal_duration = read_laser_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay;  read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        
        laser_init_ref_delay = when_read_signal_end + laser_init_delay
        when_init_ref_end    = laser_init_ref_delay + laser_init_duration
        laser_read_ref_delay = when_init_ref_end + laser_to_MWI_delay + MWI_duration + MW_to_read_delay
        read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration    = read_duration; when_read_ref_end = read_ref_delay + read_ref_duration
        laser_read_ref_duration = read_laser_duration
        self.read_duration = read_signal_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")    

        # Make pulse sequence (per each freq)
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay,    duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,    duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('MWswitch',     MWI_delay,               duration=int(MWI_duration))]
        pulse_sequence += [spc.Pulse('Counter',      read_signal_delay,       duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',      read_ref_delay,          duration=int(read_ref_duration))] # times are in ns
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
                fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.RabiROObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1

    def close_turnOnAtEnd(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        pb.turn_on_infinite(channels=channels)
    
    def sweepVelToEnd(self):
        if (self.settings['ifScanVpz']):
            nstep = int((self.vel_vpz_end-self.vel_vpz_target)/self.vel_vpz_step + 1)
            for vpz in np.linspace(self.vel_vpz_target, self.vel_vpz_end, nstep):
                vel.set_vpiezo(vpz)
                vel.waitUntilComplete()
                vel.set_ready()
                time.sleep(self.vel_vpz_step_time)
        else:
            return
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