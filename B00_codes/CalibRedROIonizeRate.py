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
from copy import deepcopy
from qcodes.actions import Task as qctask
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot
import numpy as np
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS
from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity

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

class CalibRedROIonizeRate(Instrument):

    def __init__(self, name='CalibRedROIonizeRateObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz, of the Pulse Blaster
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

        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.tausArray = np.linspace(start, stop, num_sweep_points)

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)

        self.SRSnum = self.settings['SRSnum'];      MWPower = self.settings['MWPower']; MWFreq = self.settings['MWFreq']
        
        self.velNum = self.settings['velNum']
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
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(MWFreq) #Hz
        self.srs.set_RFAmplitude(MWPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        # Velocity object
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
        
        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        global srs; srs = self.srs
        global vel; vel = self.vel
    
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
                                            sig, #ref, 
                                            qctask(sig.plotPulseSequences),
                                            )
        data = loop.get_data_set(name='CalibRedROIonizeRate')
        data.add_metadata(self.settings)
        self.data = data
        
        # plot = QtPlot(
        #     data.CalibRedROIonizeRateObject_sig, # this is implemented as a Parameter
        #     figsize = (1200, 600),
        #     interval = 1,
        #     name = 'sig'
        #     )
        # plot.add(data.CalibRedROIonizeRateObject_ref, name='ref')

        # loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        # dataPlotFilename = data.location + "/dataPlot.png"
        # dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        # img = Image.open(dataPlotFile)
        # img.show()
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/CalibRedROIonizeRateObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.CalibRedROIonizeRateObject = measurementObject
        self.loopCounter = 0
        self.ifMeaningfulRef = self.settings['ifMeaningfulRef']
        self.ifRefBright = self.settings['ifRefBright']
        self.ifRefInitAgain = self.settings['ifRefInitAgain']
        start = self.settings['start']; stop = self.settings['stop']; num_sweep_points = self.settings['num_sweep_points']
        self.tausArray = np.linspace(start, stop, num_sweep_points)

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster

        # if np.mod(self.loopCounter,100) == 0:
        #     print("Loop " + str(self.loopCounter))
        print("Loop " + str(self.loopCounter))
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];             
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration     = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];     MWI_duration            = self.settings['pi_time']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];     read_duration           = self.settings['read_duration']
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']; ref_laser_to_read_delay = self.settings['ref_laser_to_read_delay']
        delay_between_reads     = self.settings['delay_between_reads'];    nread                  = self.settings['num_reads']

        self.num_loops = num_loops

        # Make pulse sequence
        pulse_sequence = []

        when_init_end   = laser_init_delay +  laser_init_duration
        MWI_delay = when_init_end + laser_to_MWI_delay;                 when_pulse_end = MWI_delay+MWI_duration
        
        laser_read_signal_delay    = when_pulse_end
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay + tau_ns;   read_signal_duration = read_duration
        first_read_signal_delay    = laser_read_signal_delay + laser_to_DAQ_delay
        when_read_signal_end       = read_signal_delay
        laser_read_signal_duration = when_read_signal_end + DAQ_to_laser_off_delay - laser_read_signal_delay
        when_laser_read_signal_end = when_read_signal_end + DAQ_to_laser_off_delay
        
        laser_read_ref_delay = (when_laser_read_signal_end + laser_to_DAQ_delay) + (laser_to_MWI_delay + MWI_duration) + when_init_end 
        laser_init_again_delay = (when_laser_read_signal_end + laser_to_DAQ_delay) + laser_init_delay
        laser_init_again_duration = laser_init_duration

        if self.ifMeaningfulRef:
            read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay + ref_laser_to_read_delay
            read_ref_duration    = read_duration
        else:
            read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay
            read_ref_duration = 30
        when_read_ref_end = read_ref_delay + read_ref_duration
        laser_read_ref_duration = when_read_ref_end + DAQ_to_laser_off_delay - laser_read_ref_delay

        # nread = int((laser_read_signal_duration - DAQ_to_laser_off_delay - laser_to_DAQ_delay) / (read_signal_duration + delay_between_reads))
    
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,              duration=int(laser_init_duration))] # times are in ns
        # print('Laser delay in us:')
        # print(laser_read_signal_delay/1e3)
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
        # if self.ifRefBright:
        #     pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,       duration=int(laser_read_ref_duration))]
        # if self.ifRefInitAgain:
        #     pulse_sequence += [spc.Pulse('LaserInit',  laser_init_again_delay,       duration=int(laser_init_again_duration))] # times are in ns
        for i in range(nread):
            delay = first_read_signal_delay + i*(read_signal_duration + delay_between_reads)
            # print('DAQ delay in us:')
            # print(delay/1e3)
            pulse_sequence += [spc.Pulse('Counter',  delay,             duration=int(read_signal_duration))] # times are in ns
        
        # pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns
        
        self.read_duration = read_signal_duration
    
        
        
        self.pulse_sequence = pulse_sequence

        # if np.mod(self.loopCounter,100) == 0:
        #     self.ifPrintTime = True
        # else: self.ifPrintTime = False
        self.ifPrintTime = True
        
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False, ifPrintTime=self.ifPrintTime)
        self.pb.program_pb(pulse_sequence, num_loops=num_loops)

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        self.num_reads = int(num_loops * num_reads_per_iter)
        print(self.num_reads)
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

        # rate = xLineData/(self.read_duration/1e9)/1e3
        rate = xLineData
        global sig_data; 
        sig_data = np.zeros((self.num_reads_per_iter, self.num_loops))
        for i in range(self.num_reads_per_iter):
            each_read_data = rate[i::self.num_reads_per_iter]
            sig_data[i] = each_read_data
        
        global sig_avg;  sig_avg = np.average(sig_data, 1)

        return sig_avg
    
    def plotPulseSequences(self):
        self.readColor = self.settings['LaserRead']['channel']
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False, readColor=self.readColor)
                fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.CalibRedROIonizeRateObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
    
    def turn_on_at_end(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        pb.turn_on_infinite(channels=channels)


class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return 0

class SigOverRef(Parameter):
    def __init__(self, name='sigOverRef',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return 0
    
class SigFullData(Parameter):
    def __init__(self, name='sigFullData',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_data
    
class RefFullData(Parameter):
    def __init__(self, name='refFullData',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return 0