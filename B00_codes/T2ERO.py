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


class T2ERO(Instrument):

    def __init__(self, name='T2EROObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
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
        self.settings = {**settings, **settings_extra}; self.ROtrackingSettings = self.settings['ROtrackingSettings']
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
        self.ifScanVpz = self.settings['ifScanVpz']; self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

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

        if self.ifScanVpz:
            self.vel.set_vpiezo(self.vel_vpz_start)
            self.vel.waitUntilComplete()
            self.vel.set_ready()
            time.sleep(1)

            self.set_vpz()
            time.sleep(3)
        else:
            for i in range(1):
                self.vel.set_vpiezo(self.vel_vpz_target)
                self.vel.waitUntilComplete()
                self.vel.set_ready()
                time.sleep(0.7)

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

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig, ref,
                                            qctask(sig.plotPulseSequences),
                                            ).then(qctask(sig.sweepVelToEnd))

        data = loop.get_data_set(name='T2ERO')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.T2EROObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.T2EROObject_ref, name='ref')

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

        # tracking at the end if needed
        self.hasTracked = 0
        if self.settings['laserRead_channel'] == 5 or self.settings['laserRead_channel'] == 14:
            if self.ROtrackingSettings['if_tracking'] == 1:
                threshold_scanVpz = self.ROtrackingSettings['threshold_scanVpz']
                
                datafile = self.getDataFilename()
                x_s, sig, ref = dr.readDataNoPlot(datafile)
                sig = np.array(sig); ref = np.array(ref)
                self.vpz = self.vel_vpz_target

                if np.max((np.max(ref), np.max(sig))) < threshold_scanVpz: 
                    print()
                    print('-----------------Start line tracking---------------------------')
                    ScanROFreqObject = ScanROFreq(settings=self.ROtrackingSettings, ifPlotPulse=0)
                    self.vpz = ScanROFreqObject.runScanInPulseSequence()
                    vel.set_vpiezo(self.vpz)
                    vel.set_ready()
                    print("Set Vpiezo to " + str(np.round(self.vpz,1)) + ' %') # set the Velocity's piezo voltage
                    print('-----------------End line tracking---------------------------')
                    print()
                    self.timeLastROTracking = time.time()
                    self.hasTracked = 1
        
        self.srs.disable_RFOutput()
        self.srs.disableModulation()

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T2EROObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.T2EROObject = measurementObject
        self.loopCounter = 0
        self.numOfRepumpVpz = 0
        self.settings = settings
        self.tausArray = self.settings['tausArray']
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

        # Line tracking or piezo repumping for RO
        if self.settings['laserRead_channel'] == 5 or self.settings['laserRead_channel'] == 14:
            if self.ROtrackingSettings['if_tracking'] == 1:
                threshold_repumpVpz = self.ROtrackingSettings['threshold_repumpVpz']
                if sig_avg < threshold_repumpVpz:
                    if np.mod(self.numOfRepumpVpz,6) == 0:
                        print()
                        print('-----------------Start resetting Vpiezo---------------------------')
                        vel.set_vpiezo(50)
                        vel.waitUntilComplete()
                        vel.set_ready()
                        time.sleep(0.7)
                        vel.set_vpiezo(2)
                        vel.waitUntilComplete()
                        vel.set_ready()
                        time.sleep(0.7)
                        for i in range(1):
                            vel.set_vpiezo(self.vel_vpz_target)
                            vel.waitUntilComplete()
                            vel.set_ready()
                            time.sleep(0.7)
                        print('-----------------End resetting Vpiezo---------------------------')
                        print()
                        self.timeLastROTracking = time.time()
                    self.numOfRepumpVpz += 1

        return sig_avg

    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))

        NO_MS_EQUALS_1 = 0
        Q_FINAL = 1
        THREE_PI_HALF_FINAL = 2
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];              
        laser_init_delay        = self.settings['laser_init_delay'];        laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];      MWI_duration        = self.settings['pi_half']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];      read_duration       = self.settings['read_duration']   
        read_laser_duration     = self.settings['read_laser_duration'];     normalized_style    = self.settings['normalized_style']
        MW_to_read_delay        = self.settings['MW_to_read_delay']

        if tau_ns/2 > 30:
            MWI_to_switch_delay  = self.settings['MWI_to_switch_delay']
        else: 
            MWI_to_switch_delay  = 0
        
        when_init_end = laser_init_delay + laser_init_duration
        MWI_delay     = when_init_end + laser_to_MWI_delay;                 
        MWI2_delay    = MWI_delay + MWI_duration + tau_ns/2;    MWI2_duration = 2*MWI_duration
        MWI3_delay    = MWI2_delay + MWI2_duration + tau_ns/2;  MWI3_duration = MWI_duration    
        when_pulse_end= MWI3_delay + MWI3_duration;             pulse_duration= tau_ns + 4*MWI_duration

        laser_read_signal_delay    = when_pulse_end + MW_to_read_delay;            laser_read_signal_duration = read_laser_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay; read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration

        laser_init_ref_delay = when_read_signal_end + laser_init_delay
        when_init_ref_end    = laser_init_ref_delay + laser_init_duration

        MWI4_delay         = when_init_ref_end + laser_to_MWI_delay;    MWI4_duration = MWI_duration
        MWI5_delay         = MWI4_delay + MWI4_duration + tau_ns/2;     MWI5_duration = 2*MWI_duration
        MWI6_delay         = MWI5_delay + MWI5_duration + tau_ns/2;        
        if normalized_style == Q_FINAL:                                 MWI6_duration = MWI_duration
        elif normalized_style == THREE_PI_HALF_FINAL:                   MWI6_duration = 3*MWI_duration
        else:                                                           MWI6_duration = 0
        when_pulse_ref_end = MWI6_delay + MWI6_duration

        laser_read_ref_delay = when_pulse_ref_end + MW_to_read_delay;      laser_read_ref_duration = read_laser_duration
        read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration    = read_duration;                              when_read_ref_end = read_ref_delay + read_ref_duration
        self.read_duration = read_signal_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")

        # Make pulse sequence
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_delay,             duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_ref_delay,         duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay,       duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,          duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('MWswitch', MWI_delay,                     duration=int(MWI_duration))]
        pulse_sequence += [spc.Pulse('MWswitch', MWI2_delay,                    duration=int(MWI2_duration))]
        pulse_sequence += [spc.Pulse('MWswitch', MWI3_delay,                    duration=int(MWI3_duration))]
        pulse_sequence += [spc.Pulse('Counter',  read_signal_delay,             duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',  read_ref_delay,                duration=int(read_ref_duration))] # times are in ns

        if not normalized_style == NO_MS_EQUALS_1:
            pulse_sequence += [spc.Pulse('MWswitch', MWI4_delay,                     duration=int(MWI4_duration))]
            pulse_sequence += [spc.Pulse('MWswitch', MWI5_delay,                     duration=int(MWI5_duration))]

            if normalized_style == Q_FINAL:
                MWI6_duration = MWI_duration
                pulse_sequence += [spc.Pulse('MW_I', MWI6_delay-MWI_to_switch_delay, duration=int(MWI6_duration + 2*MWI_to_switch_delay))]
                pulse_sequence += [spc.Pulse('MW_Q', MWI6_delay-MWI_to_switch_delay, duration=int(MWI6_duration + 2*MWI_to_switch_delay))] # times are in ns
            elif normalized_style == THREE_PI_HALF_FINAL:
                MWI6_duration = 3*MWI_duration
            pulse_sequence += [spc.Pulse('MWswitch', MWI6_delay,                     duration=int(MWI6_duration))]
        
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
                self.T2EROObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1

    def turn_on_at_end(self):
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