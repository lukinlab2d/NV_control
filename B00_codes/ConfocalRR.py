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
from copy import deepcopy

from qcodes_contrib_drivers.drivers.NationalInstruments.DAQ import *
from qcodes_contrib_drivers.drivers.NationalInstruments.class_file import *

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
from B00_codes.ScanRRFreq import *  

def ns2cycles(time, samp_rate=1e7):
        return int(time/1e9*samp_rate)
    

class ConfocalRR(Instrument):

    def __init__(self, name='ConfocalRRObject', settings=None, ifPlotPulse=True, **kwargs) -> None:

        # clock speed is in MHz - is 'status' needed in the dictionary?
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.LaserRead2Param =   {'delay_time': 2, 'channel':settings['laserRead2_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.AWGParam =         {'delay_time': 2, 'channel':settings['AWG_channel']}
        self.MWI2Param =         {'delay_time': 2, 'channel':settings['MWI2_channel']}
        self.MWQ2Param =         {'delay_time': 2, 'channel':settings['MWQ2_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
    
        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,'AWG': self.AWGParam,
                          'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                          'LaserRead2': self.LaserRead2Param, 'LaserInit': self.LaserInitParam,
                          'MW_I2': self.MWI2Param, 'MW_Q2': self.MWQ2Param, 'MWswitch2': self.MWswitch2Param,
                          'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)
    
        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        ############################################################################################
        # Microwave and velocity 1
        self.xArray = self.settings['xArray']; self.yArray = self.settings['yArray']
        self.ifNeedSRS = self.settings['ifNeedSRS']
        self.SRSnum = self.settings['SRSnum']; MWPower = self.settings['MWPower']; MWFreq = self.settings['MWFreq']
        self.SDGnum = self.settings['SDGnum']
        self.ifNeedVel1 = self.settings['ifNeedVel1']; self.velNum = self.settings['velNum']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

        # Microwave and velocity 2
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        
        self.ifNeedVel2 = self.settings['ifNeedVel2']; self.velNum2 = self.settings['velNum2']
        vel_current2 = self.settings['vel_current2']; vel_wvl2 = self.settings['vel_wvl2']; 
        self.vel_vpz_target2 = self.settings['vel_vpz_target2']
        ###########################################################################################

        self.add_parameter(
            name = "sig2",
            parameter_class = Signal2,
            settings = self.settings,
            measurementObject = self
        )
        self.add_parameter(
            name = "sig",
            parameter_class  = Signal,
        )
        self.add_parameter(
            name = "y",
            parameter_class = VoltageY,
        )
        
        self.savedPulseSequencePlots = {}

        # SRS object
        if self.ifNeedSRS:
            self.srs = SRS(SRSnum=self.SRSnum)
            self.srs.set_freq(MWFreq) #Hz
            self.srs.set_RFAmplitude(MWPower) #dBm
            if (self.SRSnum != 3) and (self.SRSnum != 4):
                self.srs.enableIQmodulation()
            self.srs.enable_RFOutput()

            # SRS object 2
            self.srs2 = SRS(SRSnum=self.SRSnum2)
            self.srs2.set_freq(MWFreq2) #Hz
            self.srs2.set_RFAmplitude(MWPower2) #dBm
            if (self.SRSnum2 != 3) and (self.SRSnum2 != 4):
                self.srs2.enableIQmodulation()
            self.srs2.enable_RFOutput()

            global srs; srs = self.srs; global srs2; srs2 = self.srs2
        
        # AWG object
        ifAWG = self.settings['ifAWG']; self.ifAWG = ifAWG
        if self.ifAWG:
            self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum)
            global AWG; AWG = self.AWG

        # Velocity object
        if self.ifNeedVel1:
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
            
            for i in range(1):
                self.vel.set_vpiezo(self.vel_vpz_target)
                self.vel.waitUntilComplete()
                self.vel.set_ready()
                time.sleep(0.7)
            global vel; vel = self.vel; 

        # Velocity object 2
        if self.ifNeedVel2:
            self.vel2 = Velocity(velNum=self.velNum2, 
                                ifInitVpz=self.ifInitVpz, ifInitWvl=self.ifInitWvl,
                                initWvl=vel_wvl2)
            if self.ifInitWvl: 
                self.vel2.set_track()
                time.sleep(0.5)
                self.vel2.set_wvl(vel_wvl2)
                time.sleep(1)
                self.vel2.set_ready()
                self.vel2.set_vpiezo(2)
                self.vel2.waitUntilComplete()
                self.vel2.set_ready()
                time.sleep(0.7)
            if self.ifInitVpz:
                self.vel2.set_vpiezo(2)
                self.vel2.waitUntilComplete()
                self.vel2.set_ready()
                time.sleep(0.7)

            self.vel2.set_current(vel_current2)
            
            for i in range(1):
                self.vel2.set_vpiezo(self.vel_vpz_target2)
                self.vel2.waitUntilComplete()
                self.vel2.set_ready()
                time.sleep(0.7)
            global vel2; vel2 = self.vel2

        # Make Pulse Blaster, Counter, SRS global objects
        global pb

        # DAQ object to control galvo
        galvo_card_name = "cDAQ1Mod2"
        galvo_ao_channels = {f'{galvo_card_name}/ao{i}': i for i in range(4)} # dictionary of analog output channels
        self.galvo = DAQAnalogOutputs("ConfocalRRDAQ", galvo_card_name, galvo_ao_channels)
        global galvo; galvo = self.galvo
    
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        sig2 = self.sig2
        y = self.y

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loopx = Loop(
            sig2.sweep(keys=self.xArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig2,sig,
                                            qctask(sig2.plotPulseSequences),
                                            )
        loopy = Loop(
            y.sweep(keys=self.yArray),
            delay=0,
            sleepTimeAfterFinishing=0).each(loopx)

        data = loopy.get_data_set(name='ConfocalRR')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.ConfocalRRObject_sig2, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig2'
            )
        # plot.add(data.ConfocalRRObject_sig, name='sig')

        loopy.with_bg_task(plot.update, bg_final_task=None)
        loopy.run()
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

        if self.ifNeedSRS:
            self.srs.disable_RFOutput()
            self.srs.disableModulation()
            self.srs2.disable_RFOutput()
            self.srs2.disableModulation()
        if self.ifAWG: self.AWG.turn_off()
        self.galvo.close() 
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/ConfocalRRObject_sig2_set.dat'

    
class Signal2(Parameter):
    def __init__(self, name='sig2', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.loopCounter = 0
        self.numOfRepumpVpz = 0
        self.numOfRepumpVpz2 = 0
        self.settings = settings
        self.xArray = self.settings['xArray']
        self.ConfocalRRObject = measurementObject
        self.settleTime = self.settings['settleTime']

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']
        self.ifAWG = self.settings['ifAWG']

        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.vel_vpz_target2 = self.settings['vel_vpz_target2']

        # Pulse parameters
        num_loops               = self.settings['num_loops']
        laser_init_delay        = self.settings['laser_init_delay'];    laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];  MWI_duration        = self.settings['MWI_duration']
        laser_to_DAQ_delay      = self.settings['laser_to_DAQ_delay'];  read_duration       = self.settings['read_duration']   
        read_laser_duration     = self.settings['read_laser_duration']; MW_to_read_delay    = self.settings['MW_to_read_delay']
        shift_btwn_2NV_MW       = self.settings['shift_btwn_2NV_MW'];   shift_btwn_2NV_read = self.settings['shift_btwn_2NV_read']
        laser_to_DAQ_delay2     = self.settings['laser_to_DAQ_delay2']
        AWG_buffer              = self.settings['AWG_buffer'];          AWG_output_delay    = self.settings['AWG_output_delay']
        MW_duration_for_AWG = int(2*int((AWG_buffer + MWI_duration + 1)/2))

        self.MW_duration_for_AWG = MW_duration_for_AWG; self.AWG_buffer = AWG_buffer

        when_init_end              = laser_init_delay + laser_init_duration
        MWI_delay                  = when_init_end    + laser_to_MWI_delay
        global MW_del;      MW_del = MWI_delay + AWG_output_delay

        if self.ifAWG:
            when_pulse_end         = MWI_delay + AWG_output_delay + MW_duration_for_AWG
        else:
            when_pulse_end         = MWI_delay + MWI_duration
                
        laser_read_signal_delay    = when_pulse_end   + MW_to_read_delay
        laser_read_signal_duration = read_laser_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration

        MWI_delay2                  = when_init_end   + laser_to_MWI_delay + shift_btwn_2NV_MW
        when_pulse_end2             = MWI_delay2      + MWI_duration
        laser_read_signal_delay2    = when_pulse_end2 + MW_to_read_delay + shift_btwn_2NV_read
        laser_read_signal_duration2 = read_laser_duration
        when_laser_read_signal_end2 = laser_read_signal_delay2 + laser_read_signal_duration2
        read_signal_delay2          = laser_read_signal_delay2 + laser_to_DAQ_delay2
        read_signal_duration2       = read_duration
        when_read_signal_end2       = read_signal_delay2 + read_signal_duration2
        
        # laser_init_ref_delay        = np.max((when_read_signal_end2, when_read_signal_end)) + laser_init_delay
        # when_init_ref_end           = laser_init_ref_delay + laser_init_duration
        # laser_read_ref_delay        = when_init_ref_end + laser_to_MWI_delay + MWI_duration + MW_to_read_delay
        # laser_read_ref_duration     = read_laser_duration
        # when_laser_read_ref_end     = laser_read_ref_delay + laser_read_ref_duration
        # read_ref_delay              = laser_read_ref_delay + laser_to_DAQ_delay;  
        # read_ref_duration           = read_duration
        # when_read_ref_end           = read_ref_delay + read_ref_duration

        # laser_read_ref_delay2    = when_init_ref_end + laser_to_MWI_delay + shift_btwn_2NV_MW + MWI_duration + MW_to_read_delay + shift_btwn_2NV_read
        # laser_read_ref_duration2 = read_laser_duration
        # when_laser_read_ref_end2 = laser_read_ref_delay2 + laser_read_ref_duration2
        # read_ref_delay2          = laser_read_ref_delay2 + laser_to_DAQ_delay2
        # read_ref_duration2       = read_duration
        # when_read_ref_end2       = read_ref_delay2 + read_ref_duration2

        self.read_duration = read_signal_duration
        # if read_signal_duration != read_ref_duration:
        #     raise Exception("Duration of reading signal and reference must be the same")    
        

        # Make pulse sequence (per each freq)
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
            # pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay,    duration=int(laser_init_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] # times are in ns
        # pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay,    duration=int(laser_read_ref_duration))]
        # if self.ifAWG:
        #     pulse_sequence += [spc.Pulse('AWG',          MWI_delay,    duration=20)]
        # else:
        #     pulse_sequence += [spc.Pulse('MWswitch',     MWI_delay,           duration=int(MWI_duration))]
        pulse_sequence += [spc.Pulse('Counter',      read_signal_delay,       duration=int(read_signal_duration))] # times are in ns
        # pulse_sequence += [spc.Pulse('Counter',      read_ref_delay,          duration=int(read_ref_duration))] # times are in ns
        
        pulse_sequence += [spc.Pulse('LaserRead2',    laser_read_signal_delay2, duration=int(laser_read_signal_duration2))] # times are in ns
        # pulse_sequence += [spc.Pulse('LaserRead2',    laser_read_ref_delay2,    duration=int(laser_read_ref_duration2))]
        # pulse_sequence += [spc.Pulse('MWswitch2',     MWI_delay2,               duration=int(MWI_duration))]
        pulse_sequence += [spc.Pulse('Counter',       read_signal_delay2,       duration=int(read_signal_duration2))] 
        # pulse_sequence += [spc.Pulse('Counter',       read_ref_delay2,          duration=int(read_ref_duration2))] 
        
        self.pulse_sequence = pulse_sequence; self.num_loops = num_loops

        num_reads_per_iter = 0
        for pulse in pulse_sequence:
            if pulse.channel_id == 'Counter': num_reads_per_iter +=1
        self.num_reads = int(num_loops * num_reads_per_iter)
        self.num_reads_per_iter = num_reads_per_iter

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))
        self.ctrtask.stop(); self.ctrtask.close()

        rate = xLineData/(self.read_duration/1e9)/1e3
        sig2 = rate[::self.num_reads_per_iter]
        sig = rate[1::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)
        global sig_avg2;  sig_avg2 = np.average(sig2)
                    
        return sig_avg2
    
    # def set_coordinate_fnc(x,y,z):
    #     galvo.voltage_cdaq1mod2ao0(x) # this is for the x-mirror
    #     galvo.voltage_cdaq1mod2ao1(y) # this is for the y-mirror
    #     galvo.voltage_cdaq1mod2ao2(z) # this is for the z-piezo
    
    def set_raw(self, x):
        print("Loop " + str(self.loopCounter))

        # Set galvo
        galvo.voltage_cdaq1mod2ao0(x)
        time.sleep(self.settleTime/1e9)

        # Program Pulse Blaster
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False)
        self.pb.program_pb(self.pulse_sequence, num_loops=self.num_loops)

        if self.ifAWG:
            global ch1plot; global ch2plot
            ch1plot, ch2plot = AWG.send_fastRabi_seq(pulse_width=int(self.MW_duration_for_AWG), buffer=int(self.AWG_buffer))

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

        print('Set x to ' + str(x) + " V")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1
    
    def plotPulseSequences(self):
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.xArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False,
                                            readColor=self.readColor, 
                                            initColor=self.initColor)
                if self.ifAWG:
                    fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, MW_del)
                else:
                    fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.xArray)-1: # only save first and last pulse sequence
                self.ConfocalRRObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1

    def close_turnOnAtEnd(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        # pb.turn_on_infinite(channels=channels)
    

class Signal(Parameter):
    def __init__(self, name='sig',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_avg
    
class VoltageY(Parameter):
    def __init__(self, name='y', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.ConfocalRRObject = measurementObject

    def set_raw(self,y):
        galvo.voltage_cdaq1mod2ao1(y)
        print('Set y to ' + str(y) + " V")