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
from PlotPulse import *    
from Confocal import *

class XY8(Instrument):

    def __init__(self, name='XY8Object', settings=None, ifPlotPulse=True, tausArray=None, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':1}
        self.MWQParam =         {'delay_time': 2, 'channel':0}
        self.MWswitchParam =    {'delay_time': 2, 'channel':2}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
        global pb

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,
                        'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.tausArray = self.settings['tausArray']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)

        self.uwPower = self.settings['uwPower']; self.uwFreq = self.settings['uwFreq']

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
        self.srs = SRS()
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

        data = loop.get_data_set(name='XY8')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.XY8Object_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.XY8Object_ref, name='ref')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        img = Image.open(dataPlotFile)
        img.show()

        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/XY8Object_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.trackingSettings = self.settings['trackingSettings']
        self.XY8Object = measurementObject
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

        NO_MS_EQUALS_1 = 0
        Q_FINAL = 1
        THREE_PI_HALF_FINAL = 2
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];               ifStartInY          = self.settings['ifStartInY']
        laser_init_delay        = self.settings['laser_init_delay'];        laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay     = self.settings['laser_to_MWI_delay'];     
        piHalf                  = self.settings['pi_half'];                 pi = 2*piHalf
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];      read_duration       = self.settings['read_duration']   
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay'];  normalized_style    = self.settings['normalized_style']

        if tau_ns/2 > 30:
            MWI_to_switch_delay  = self.settings['MWI_to_switch_delay']
        else: 
            MWI_to_switch_delay  = 0
        
        when_init_end = laser_init_delay + laser_init_duration
        MW_delays     = []
        MW_delays.append(when_init_end + laser_to_MWI_delay)         
        MW_delays.append(MW_delays[0] + piHalf + tau_ns/2)
        for i in range(2,9):
            MW_delays.append(MW_delays[i-1] + pi + tau_ns)
        MW_delays.append(MW_delays[8] + pi + tau_ns/2)

        when_pulse_end = MW_delays[9] + piHalf;             pulse_duration = 8*tau_ns + 9*pi

        laser_read_signal_delay    = when_pulse_end
        read_signal_delay          = when_pulse_end + laser_to_DAQ_delay;   read_signal_duration = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        laser_read_signal_duration = when_read_signal_end + DAQ_to_laser_off_delay - laser_read_signal_delay
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration

        MW_delays.append(when_laser_read_signal_end + laser_to_MWI_delay)
        MW_delays.append(MW_delays[10] + piHalf + tau_ns/2)
        for i in range(12,19):
            MW_delays.append(MW_delays[i-1] + pi + tau_ns)   
        MW_delays.append(MW_delays[18] + pi + tau_ns/2)
        
        laser_read_ref_delay = when_laser_read_signal_end + laser_to_MWI_delay + pulse_duration
        read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration    = read_duration; when_read_ref_end = read_ref_delay + read_ref_duration
        laser_read_ref_duration = when_read_ref_end + DAQ_to_laser_off_delay - laser_read_ref_delay
        self.read_duration = read_signal_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")
        
        if ifStartInY:
            MWI_list1 = np.array((0,2,4,5,7))
        else: 
            MWI_list1 = np.array((2,4,5,7))

        MWI_list2 = np.array((12,14,15,17))

        # Make pulse sequence ?-XYXYYXYX-?
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_delay,        duration=int(laser_init_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',        laser_read_signal_delay, duration=int(laser_read_signal_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',        laser_read_ref_delay,    duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('MWswitch',     MW_delays[0],            duration=int(piHalf))]
        for i in range(1,9):
            pulse_sequence += [spc.Pulse('MWswitch', MW_delays[i],            duration=int(pi))]
        pulse_sequence += [spc.Pulse('MWswitch',     MW_delays[9],            duration=int(piHalf))]
        for i in MWI_list1:
            if i == 0:
                pulse_sequence += [spc.Pulse('MW_I',     MW_delays[i]-MWI_to_switch_delay, duration=int(piHalf + 2*MWI_to_switch_delay))]
            else:
                pulse_sequence += [spc.Pulse('MW_I',     MW_delays[i]-MWI_to_switch_delay, duration=int(pi + 2*MWI_to_switch_delay))]
        pulse_sequence += [spc.Pulse('Counter',      read_signal_delay,       duration=int(read_signal_duration))] 
        pulse_sequence += [spc.Pulse('Counter',      read_ref_delay,          duration=int(read_ref_duration))] 

        if not normalized_style == NO_MS_EQUALS_1:
            pulse_sequence += [spc.Pulse('MWswitch', MW_delays[10],           duration=int(piHalf))]
            for i in range(11,19):
                pulse_sequence += [spc.Pulse('MWswitch',MW_delays[i],         duration=int(pi))]
            
            if normalized_style == Q_FINAL:
                MWI19_duration = piHalf
                pulse_sequence += [spc.Pulse('MW_Q', MW_delays[19]-MWI_to_switch_delay, duration=int(MWI19_duration + 2*MWI_to_switch_delay))]
                if not ifStartInY:
                    pulse_sequence += [spc.Pulse('MW_I', MW_delays[19]-MWI_to_switch_delay, duration=int(MWI19_duration + 2*MWI_to_switch_delay))]
            elif normalized_style == THREE_PI_HALF_FINAL:
                MWI19_duration = 3*piHalf

            pulse_sequence += [spc.Pulse('MWswitch', MW_delays[19],           duration=int(MWI19_duration))]
            for i in MWI_list2:
                pulse_sequence += [spc.Pulse('MW_I', MW_delays[i]-MWI_to_switch_delay, duration=int(pi + 2*MWI_to_switch_delay))]
        
        self.pulse_sequence = pulse_sequence
        # print(pulse_sequence)
        
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
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False, plotFilename='C:/Users/lukin2dmaterials/data/2022-12-14/aaa.png')
                fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.XY8Object.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
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

class SigOverRef(Parameter):
    def __init__(self, name='sigOverRef',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_avg_over_ref_avg