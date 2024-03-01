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

from RohdeSchwarz.SMC import *

def ns2cycles(time, samp_rate=1e7):
        return int(time/1e9*samp_rate)
    

class ScanRRFreqPumpFree(Instrument):

    def __init__(self, name='ScanRRFreqPumpFreeObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        self.ifPlotPulse=ifPlotPulse

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz, of the Pulse Blaster
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        self.hiLoMWPwrParam =    {'delay_time': 2, 'channel':settings['hiLoMWPwr_channel']}
        self.MWswitch3Param =    {'delay_time': 2, 'channel':settings['MWswitch3_channel']}
        self.MWswitch4Param =    {'delay_time': 2, 'channel':settings['MWswitch4_channel']}

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 
                          'LaserRead': self.LaserReadParam,
                        'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'MWswitch2': self.MWswitch2Param, 'hiLoMWPwr': self.hiLoMWPwrParam,
                        'MWswitch3': self.MWswitch3Param, 'MWswitch4': self.MWswitch4Param, 
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.vpzArray = self.settings['vpzArray']
        self.SRSnum = self.settings['SRSnum'];   MWPower = self.settings['MWPower'];   MWFreq = self.settings['MWFreq']
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        self.SRSnum3 = self.settings['SRSnum3']; MWPower3 = self.settings['MWPower3']; MWFreq3 = self.settings['MWFreq3']
        self.SRSnum4 = self.settings['SRSnum4']; MWPower4 = self.settings['MWPower4']; MWFreq4 = self.settings['MWFreq4']

        self.ifNeedVel = self.settings['ifNeedVel']; self.velNum = self.settings['velNum']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

        self.readColor    = self.settings['LaserRead']['channel']

        
        # Pulse parameters
        num_loops            = self.settings['num_loops']
        read_laser_delay     = self.settings['read_laser_delay'];    read_laser_duration = self.settings['read_laser_duration']
        laser_to_DAQ_delay   = self.settings['laser_to_DAQ_delay'];  DAQ_duration        = self.settings['DAQ_duration']   
        
        DAQ_delay = read_laser_delay + laser_to_DAQ_delay
        MW_delay  = read_laser_delay;                                MW_duration         = read_laser_duration

        # Make pulse sequence (per each freq)
        pulse_sequence = []
        pulse_sequence += [spc.Pulse('LaserRead',    read_laser_delay,        duration=int(read_laser_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('MWswitch',     MW_delay,                duration=int(MW_duration))]
        pulse_sequence += [spc.Pulse('MWswitch2',    MW_delay,                duration=int(MW_duration))]
        pulse_sequence += [spc.Pulse('MWswitch3',    MW_delay,                duration=int(MW_duration))]
        pulse_sequence += [spc.Pulse('MWswitch4',    MW_delay,                duration=int(MW_duration))]
        pulse_sequence += [spc.Pulse('hiLoMWPwr',    MW_delay,                duration=int(MW_duration))]
        pulse_sequence += [spc.Pulse('Counter',      DAQ_delay,               duration=int(DAQ_duration))] # times are in ns
        self.pulse_sequence = pulse_sequence
        
        # SRS object
        if True:
            self.srs = SRS(SRSnum=self.SRSnum)
            self.srs.set_freq(MWFreq) #Hz
            self.srs.set_RFAmplitude(MWPower) #dBm
            if (self.SRSnum != 3) and (self.SRSnum != 4):
                self.srs.enableIQmodulation()
            self.srs.enable_RFOutput()
            global srs; srs = self.srs

            self.srs2 = SRS(SRSnum=self.SRSnum2)
            self.srs2.set_freq(MWFreq2) #Hz
            self.srs2.set_RFAmplitude(MWPower2) #dBm
            if (self.SRSnum2 != 3) and (self.SRSnum2 != 4):
                self.srs2.enableIQmodulation()
            self.srs2.enable_RFOutput()
            global srs2; srs2 = self.srs2
            
            self.srs3 = SRS(SRSnum=self.SRSnum3)
            self.srs3.set_freq(MWFreq3) #Hz
            self.srs3.set_RFAmplitude(MWPower3) #dBm
            if (self.SRSnum3 != 3) and (self.SRSnum3 != 4):
                self.srs3.enableIQmodulation()
            self.srs3.enable_RFOutput()
            global srs3; srs3 = self.srs3

            self.srs4 = SRS(SRSnum=self.SRSnum4)
            self.srs4.set_freq(MWFreq4) #Hz
            self.srs4.set_RFAmplitude(MWPower4) #dBm
            if (self.SRSnum4 != 3) and (self.SRSnum4 != 4):
                self.srs4.enableIQmodulation()
            self.srs4.enable_RFOutput()
            global srs4; srs4 = self.srs4

        # Velocity object
        if self.ifNeedVel:
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
            global vel; vel = self.vel
        

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
            DAQ_duration = DAQ_duration,
            pulse_sequence = pulse_sequence,
            settings = self.settings
        )
        
        # Make Pulse Blaster, Counter, SRS global objects
        global pb
        global ctrtask; ctrtask = self.ctrtask
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
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
            sleepTimeAfterFinishing=0).each(sig,wvlFromWM,
                                            ).then(qctask(sig.close))

        data = loop.get_data_set(name='ScanRRFreqPumpFree')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.ScanRRFreqPumpFreeObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        img = Image.open(dataPlotFile)
        img.show()
        
        if self.ifPlotPulse:
            pulsePlotFilename = data.location + "/pulsePlot.png"
            plotPulseObject = PlotPulse(measurementObject=self, plotFilename=pulsePlotFilename, ifShown=True,
                                        readColor=self.readColor)
            plotPulseObject.makePulsePlot()
        
        self.srs.disable_RFOutput()
        self.srs.disableModulation()

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ScanRRFreqPumpFreeObject_sig_set.dat'

    
class Signal(Parameter):
    def __init__(self, num_reads: int, num_reads_per_iter: int, num_loops: int,
                 DAQ_duration: int, name='sig', pulse_sequence=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.num_reads = num_reads
        self.num_reads_per_iter = num_reads_per_iter
        self.num_loops = num_loops
        self.DAQ_duration = DAQ_duration
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
        
        rate = xLineData/(self.DAQ_duration/1e9)/1e3
        sig = rate
        global sig_avg;  sig_avg = np.average(sig)

        return sig_avg
    
    def set_raw(self, vpz):
        for i in range(self.settings['num_of_cavity_conditioning']):
            vel.set_vpiezo(np.round(vpz,6))
            vel.waitUntilComplete()
            vel.set_ready()
            time.sleep(0.3)
        print("Loop " + str(self.loopCounter))
        print("Set Vpiezo to " + str(np.round(vpz,1)) + ' %') # set the Velocity's piezo voltage

    def close(self):
        ctrtask.close()
       
    
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
        elif self.velNum==2:
            dropbox_path = '/wlm_laser2.txt'
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

        return lastData