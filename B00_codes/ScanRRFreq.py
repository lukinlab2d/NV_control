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

from qcodes.dataset import do1d, do2d, dond, LinSweep, LogSweep, ArraySweep
from qcodes.utils.dataset.doNd import plot
from qcodes.dataset.sqlite.database import initialise_or_create_database_at
from qcodes.dataset.experiment_container import load_or_create_experiment
from qcodes.tests.instrument_mocks import DummyInstrument, DummyInstrumentWithMeasurement
from qcodes.dataset.measurements import Measurement
from qcodes.dataset.plotting import plot_dataset

from RohdeSchwarz.SMC import *

def ns2cycles(time, samp_rate=1e7):
        return int(time/1e9*samp_rate)
    

class ScanRRFreq(Instrument):

    def __init__(self, name='ScanRRFreqObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        self.ifPlotPulse=ifPlotPulse

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.if2sources   = settings['if2sources']

        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        if self.if2sources:
            self.AWGParam2 =         {'delay_time': 2, 'channel':settings['AWG_channel2']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
    
        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam, 'AWG': self.AWGParam,
                          'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                          'min_pulse_dur': int(5*1e3/self.clock_speed)}
        if self.if2sources:
            settings_extra['AWG2'] = self.AWGParam2
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        # Vpiezos, MW power, and MW frequency
        self.vpzArray = self.settings['vpzArray']
        self.SRSnum = self.settings['SRSnum'];   MWPower = self.settings['MWPower'];   MWFreq = self.settings['MWFreq']
        if self.if2sources:
            self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
            self.SDGnum2 = self.settings['SDGnum2']
        self.velNum = self.settings['velNum']; self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl'] 
        self.SDGnum = self.settings['SDGnum']
        self.ifIQ = self.settings['ifIQ']

        # SRS object
        self.srs = SRS(SRSnum=self.SRSnum)
        self.srs.set_freq(MWFreq) #Hz
        self.srs.set_RFAmplitude(MWPower) #dBm
        if self.ifIQ:
            self.srs.enableIQmodulation()
        self.srs.enable_RFOutput()

        if self.if2sources:
            self.srs2 = SRS(SRSnum=self.SRSnum2)
            self.srs2.set_freq(MWFreq2) #Hz
            self.srs2.set_RFAmplitude(MWPower2) #dBm
            if self.ifIQ:
                self.srs2.enableIQmodulation()
            self.srs2.enable_RFOutput()

        # AWG object
        ifAWG = self.settings['ifAWG']; self.ifAWG = ifAWG
        if self.ifAWG:
            self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum)
            global AWG; AWG = self.AWG
            if self.if2sources:
                self.AWG2 = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum2)
                global AWG2; AWG2 = self.AWG2

        # Pulse parameters
        num_loops               = self.settings['num_loops'];           AWG_buffer           = self.settings['AWG_buffer']
        laser_init_delay        = self.settings['laser_init_delay'];    laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];  MW_duration        = self.settings['MW_duration']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];  read_duration       = self.settings['read_duration']   
        read_laser_duration     = self.settings['read_laser_duration']; MW_to_read_delay    = self.settings['MW_to_read_delay']
        AWG_output_delay        = self.settings['AWG_output_delay'];   
        if self.if2sources:
            MW_duration2        = self.settings['MW_duration2']
            shift_btwn_2NV_MW   = self.settings['shift_btwn_2NV_MW']
            AWG_output_delay2   = AWG_output_delay
        
        when_init_end              = laser_init_delay + laser_init_duration
        MWI_delay                  = when_init_end    + laser_to_MWI_delay; self.MWI_delay = MWI_delay
        MW_delay_for_AWG           = MWI_delay-AWG_output_delay
        MW_duration_for_AWG        = int(2*int((2*AWG_buffer + MW_duration + 1)/2)) # to make it even

        if ifAWG:
            when_pulse_end = MWI_delay + MW_duration_for_AWG
        else:
            when_pulse_end = MWI_delay + MW_duration
        
        if self.if2sources:
            MWI_delay2            = when_init_end  + shift_btwn_2NV_MW; self.MWI_delay2 = MWI_delay2
            MW_delay_for_AWG2     = MWI_delay2 - AWG_output_delay2
            MW_duration_for_AWG2  = int(2*int((2*AWG_buffer + MW_duration2 + 1)/2)) # to make it even

            if ifAWG:
                when_pulse_end2 = MWI_delay2 + MW_duration_for_AWG2
            else:
                when_pulse_end2 = MWI_delay2 + MW_duration2
        else: when_pulse_end2=0
        
        laser_read_signal_delay    = np.max((when_pulse_end,when_pulse_end2))   + MW_to_read_delay  
        laser_read_signal_duration = read_laser_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        
        laser_init_ref_delay = when_laser_read_signal_end + laser_init_delay
        when_init_ref_end    = laser_init_ref_delay + laser_init_duration
        # laser_read_ref_delay = when_init_ref_end + laser_to_MWI_delay + MW_duration + MW_to_read_delay
        laser_read_ref_delay = laser_read_signal_delay + (laser_init_ref_delay-laser_init_delay)
        read_ref_delay       = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration    = read_duration;       when_read_ref_end = read_ref_delay + read_ref_duration
        laser_read_ref_duration = read_laser_duration

        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")  

        if ifAWG:
            global ch1plot; global ch2plot
            ch1plot, ch2plot = AWG.send_fastRabi_seq(pulse_width=int(MW_duration), buffer=int(AWG_buffer)) 
        if self.if2sources:
            global ch1plot2; global ch2plot2
            ch1plot2, ch2plot2 = AWG2.send_fastRabi_seq(pulse_width=int(MW_duration2), buffer=int(AWG_buffer))  

        # Make pulse sequence (per each freq)
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay,    duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,    duration=int(laser_read_ref_duration))]
        if ifAWG:
            pulse_sequence += [spc.Pulse('AWG',          MW_delay_for_AWG,    duration=50)]
        if self.if2sources:
            pulse_sequence += [spc.Pulse('AWG2',         MW_delay_for_AWG2,    duration=50)]
        else:
            pulse_sequence += [spc.Pulse('MWswitch',     MWI_delay,           duration=int(MW_duration))]
        pulse_sequence += [spc.Pulse('Counter',      read_signal_delay,       duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',      read_ref_delay,          duration=int(read_ref_duration))] # times are in ns
        self.pulse_sequence = pulse_sequence

        # Velocity objects
        self.vel = Velocity(velNum=self.velNum, ifInitVpz=self.ifInitVpz, ifInitWvl=self.ifInitWvl, initWvl=vel_wvl)
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
            samps_per_chan = 2*num_reads # x2 to make sure buffer doesn't overflow
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
        global vel; vel = self.vel
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        self.add_parameter(
            name = "wvlFromWM",
            parameter_class = WvlFromWM,
            settings = self.settings,
        )
        wvlFromWM = self.wvlFromWM

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.vpzArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,wvlFromWM,
                                            ).then(qctask(sig.close))

        data = loop.get_data_set(name='ScanRRFreq')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.ScanRRFreqObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.ScanRRFreqObject_ref, name='ref')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        # img = Image.open(dataPlotFile)
        # img.show()
        
        if self.ifPlotPulse:
            pulsePlotFilename = data.location + "/pulsePlot.png"
            plotPulseObject = PlotPulse(measurementObject=self, plotFilename=pulsePlotFilename, ifShown=True,
                                        readColor=self.readColor, initColor=self.initColor)
            if self.ifAWG:
                fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, self.MWI_delay, label1='chnl1-AWG1', label2='chnl2-AWG1')
                if self.if2sources:
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot2, ch2plot2, self.MWI_delay2,
                                                           fig=fig, offset1=17.75, offset2=17.5, label1='chnl1-AWG2', label2='chnl2-AWG2')
            else:
                plotPulseObject.makePulsePlot()
        
        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        if self.ifAWG: self.AWG.turn_off()

        if self.if2sources:
            self.srs2.disable_RFOutput()
            self.srs2.disableModulation()
            if self.ifAWG: self.AWG2.turn_off()

    def runScanInPulseSequence(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigOverRef = self.sigOverRef

        loop = Loop(
            sig.sweep(keys=self.vpzArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,sigOverRef,
                                            ).then(qctask(sig.close))

        data = loop.get_data_set(name='ScanRRFreq')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.ScanRRFreqObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.ScanRRFreqObject_ref, name='ref')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        plot.save(filename=dataPlotFilename, type='data')

        # Read in data and find the piezo voltage to maximize contrast
        datafile = self.getDataFilename()
        x_s, sig, ref = dr.readDataNoPlot(datafile)
        sig = np.array(sig); ref = np.array(ref)
        contrast = np.abs(sig-ref)

        return x_s[np.argmax(contrast)]
    
    def runScanInPulseSequenceMonitoringWM(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        self.add_parameter(
            name = "wvlFromWM",
            parameter_class = WvlFromWM,
            settings = self.settings,
        )
        wvlFromWM = self.wvlFromWM

        loop = Loop(
            sig.sweep(keys=self.vpzArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,ref,wvlFromWM,
                                            ).then(qctask(sig.close))

        data = loop.get_data_set(name='ScanRRFreq')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.ScanRRFreqObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.ScanRRFreqObject_ref, name='ref')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        plot.save(filename=dataPlotFilename, type='data')

        # Read in data and find the piezo voltage to maximize ref
        datafile = self.getDataFilename()
        x_s, sig, ref, wvl = dr.readDataNoPlotWM(datafile)
        sig = np.array(sig); ref = np.array(ref); wvl = np.array(wvl)

        return x_s[np.argmax(ref)], wvl[np.argmax(ref)]
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ScanRRFreqObject_sig_set.dat'

    
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

        return sig_avg
    
    def set_raw(self, vpz):
        for i in range(self.settings['num_of_cavity_conditioning']):
            vel.set_vpiezo(np.round(vpz,6))
            vel.waitUntilComplete()
            vel.set_ready()
            time.sleep(0.3)
        print("Loop " + str(self.loopCounter))
        print("Set Vpiezo to " + str(np.round(vpz,1)) + ' %') # set the Velocity's piezo voltage

    def close_turnOnAtEnd(self):
        ctrtask.close()
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        pb.turn_on_infinite(channels=channels)

    def close(self):
        ctrtask.close()
       
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
    
class WvlFromWM(Parameter):
    def __init__(self, name='wvlFromWM', settings=None,**kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.dbx = self.settings['dbx']
        self.velNum = self.settings['velNum']
        if 'timeSleepReadWvl' in self.settings:
            self.timeSleepReadWvl = self.settings['timeSleepReadWvl']
        else:
            self.timeSleepReadWvl = 0.6

    def get_raw(self):
        time.sleep(self.timeSleepReadWvl)
        if self.velNum==1:
            dropbox_path = '/wlm_laser1.txt'
        # elif self.velNum==2:
        #     dropbox_path = '/wlm_laser2.txt'
            _, response = self.dbx.files_download(dropbox_path)
            file_content = response.content

            # Convert bytes to string
            data_str = file_content.decode('utf-8')

            # Split the string into lines
            lines = data_str.split('\r\n')

            # Convert each line to a float
            data_float = [float(line) for line in lines if line]
            data = np.array(data_float)
            lastData = np.round(data[-1],6)
        elif self.velNum==2:
            filepath = 'C:/Users/lukin2dmaterials/pylabnet/b00_wlm//wlm_laser2.txt'

            data = []
            with open(filepath, 'r') as file:
                for line in file:
                    data.append(float(line.strip()))
            data = np.array(data)
            lastData = np.round(data[-1],6)

        return lastData