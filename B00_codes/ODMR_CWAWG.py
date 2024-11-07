"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
from  qcodes.actions import Task as qctask
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
from B00_codes.PlotPulseNew import * ###
from B00_codes.Confocal import *


class ODMR_CWAWG(Instrument):

    def __init__(self, name='ODMRAWGObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        self.ifPlotPulse=ifPlotPulse

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

        # List of frequencies and power
        self.freqsArray = self.settings['freqsArray']
        uwPower = self.settings['uwPower']
        self.SRSnum=self.settings['SRSnum']
        self.SDGnum=self.settings['SDGnum']

        # Pulse lengths
        AWGbuffer             = self.settings['AWGbuffer']
        num_loops             = self.settings['num_loops'];               wait_btwn_sig_ref    = self.settings['wait_btwn_sig_ref']
        AWG_output_delay      = self.settings['AWG_output_delay'];        MW_duration          = self.settings['MW_duration']    
        read_signal_delay     = AWG_output_delay;                         read_signal_duration = MW_duration + 2*AWGbuffer + self.settings['MW_off_to_read_signal_off']
        when_read_signal_ends = read_signal_delay + read_signal_duration; read_ref_delay       = when_read_signal_ends + wait_btwn_sig_ref;  
        read_ref_duration     = read_signal_duration;                                      
        when_read_ref_ends    = read_ref_delay + read_ref_duration;       laser_delay          = self.settings['laser_delay'];                
        laser_duration        = when_read_ref_ends - laser_delay + 100
        global MW_del; MW_del = AWG_output_delay ###

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")    

        # Make pulse sequence (per each freq)
        pulse_sequence = []
        pulse_sequence += [spc.Pulse('LaserInit',  laser_delay,  duration=int(laser_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('AWG', 0, duration=20)]
        pulse_sequence += [spc.Pulse('Counter', read_signal_delay,   duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter', read_ref_delay,  duration=int(read_ref_duration))] # times are in ns
        self.pulse_sequence = pulse_sequence
        
        # SRS object
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(3e9) #Hz
        self.srs.set_RFAmplitude(uwPower) #dBm
        self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        # AWG object
        self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum)
        global AWG; AWG = self.AWG
        ch1, ch2 = AWG.send_ODMR_seq(pulse_width=int(MW_duration), buffer=int(AWGbuffer)) ###
        global ch1plot; global ch2plot ###
        ch1plot = ch1; ch2plot = ch2 ###

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        num_reads = int(num_loops * num_reads_per_iter) # the results look like this: sig-ref-sig-ref-...-sig-ref for num_loops times

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
            samps_per_chan = int(1.2*num_reads) # x2 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        self.add_parameter(
            name = "sig",
            parameter_class = Signal,
            num_reads = num_reads,
            num_loops = num_loops,
            num_reads_per_iter = num_reads_per_iter,
            read_duration = read_signal_duration,
            pulse_sequence = pulse_sequence,
            settings = self.settings
        )
        self.add_parameter(
            name = "ref",
            parameter_class = Reference,
        )
        self.add_parameter(
            name = "sigOverRef",
            parameter_class = SigOverRef,
        )

        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        global ctrtask; ctrtask = self.ctrtask
        global srs; srs = self.srs
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigOverRef = self.sigOverRef

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.freqsArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,sigOverRef).then(qctask(sig.close))

        data = loop.get_data_set(name='ODMR_CWAWG')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.ODMRAWGObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.ODMRAWGObject_ref, name='ref')
        # plot = QtPlot(
        #     data.ODMRAWGObject_sigOverRef, # this is implemented as a Parameter
        #     figsize = (1200, 600),
        #     interval = 1,
        #     name = 'sig/ref'
        #     )

        loop.with_bg_task(plot.update)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        # img = Image.open(dataPlotFile)
        # img.show()
        
        if self.ifPlotPulse:
            pulsePlotFilename = data.location + "/pulsePlot.png"
            plotPulseObject = PlotPulse(measurementObject=self, plotFilename=pulsePlotFilename, ifShown=True)
            plotPulseObject.makePulsePlot(ch1plot, ch2plot, MW_del) ###
        
        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        self.AWG.turn_off()
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ODMRAWGObject_sig_set.dat'

    
class Signal(Parameter):
    def __init__(self, num_reads: int, num_reads_per_iter: int, num_loops: int,
                 read_duration: int, name='sig', pulse_sequence=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.num_reads = num_reads
        self.num_reads_per_iter = num_reads_per_iter
        self.num_loops = num_loops
        self.read_duration = read_duration
        self.loopCounter = 0
        self.pulse_sequence = pulse_sequence
        self.settings = settings
        self.trackingSettings = self.settings['trackingSettings']

    def get_raw(self):
        self.loopCounter += 1
        ctrtask.start()

        # Pulse Blaster object inittialized here to avoid long delay between initialization and first run
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        self.pb.program_pb(self.pulse_sequence, num_loops=self.num_loops)
        pb = self.pb

        pb.start_pulse_seq()
        pb.wait()
        pb.stop_pulse_seq()
        pb.close()

        xLineData = np.array(ctrtask.read(self.num_reads, timeout=10))

        ctrtask.stop()
        
        rate = xLineData/(self.read_duration/1e9)/1e3
        sig = rate[::self.num_reads_per_iter]
        ref = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global ref_avg;  ref_avg = np.average(ref)
        global sig_avg_over_ref_avg; sig_avg_over_ref_avg = sig_avg/ref_avg

        # NV tracking
        if self.trackingSettings['if_tracking'] == 1:
            if np.mod(self.loopCounter, self.trackingSettings['tracking_period']) == self.trackingSettings['tracking_period']-1:
                print()
                cfcObject = Confocal(settings=self.trackingSettings)
                cfcObject.optimize_xy()
                time.sleep(1)
                cfcObject.optimize_xz()
                time.sleep(1)
                cfcObject.optimize_xy()
                time.sleep(1)
                cfcObject.close()

        return sig_avg
    
    def set_raw(self, value):
        srs.set_freq(value)
        print("Loop " + str(self.loopCounter))
        print("Set SRS freq to " + str(np.round(value/1e6,3)) + ' MHz') # set the SRS frequency in reality

    def close(self):
        ctrtask.close()
        # pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        # channels = np.linspace(laserInitChannel,laserInitChannel,1)
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