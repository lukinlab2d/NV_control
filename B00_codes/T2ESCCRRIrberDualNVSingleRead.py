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
from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity

import nidaqmx, time, dropbox
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
from B00_codes.ScanRRFreq import *  

class T2ESCCRRIrberDualNVSingleRead(Instrument):

    def __init__(self, name='T2ESCCRRIrberDualNVSingleReadObject', settings=None, ifPlotPulse=True, **kwargs) -> None:
        
        super().__init__(name, **kwargs)
        self.clock_speed = 500 # MHz, of the Pulse Blaster
        self.LaserInitParam =   {'delay_time': 2, 'channel':settings['laserInit_channel']}
        self.LaserReadParam =   {'delay_time': 2, 'channel':settings['laserRead_channel']}
        self.LaserRead2Param =  {'delay_time': 2, 'channel':settings['laserRead2_channel']}
        self.LaserIonParam =    {'delay_time': 2, 'channel':settings['laserIon_channel']}
        self.CounterParam =     {'delay_time': 2, 'channel':4}
        self.MWIParam =         {'delay_time': 2, 'channel':settings['MWI_channel']}
        self.MWQParam =         {'delay_time': 2, 'channel':settings['MWQ_channel']}
        self.MWswitchParam =    {'delay_time': 2, 'channel':settings['MWswitch_channel']}
        self.MWI2Param =         {'delay_time': 2, 'channel':settings['MWI2_channel']}
        self.MWQ2Param =         {'delay_time': 2, 'channel':settings['MWQ2_channel']}
        self.MWswitch2Param =    {'delay_time': 2, 'channel':settings['MWswitch2_channel']}
        self.hiLoMWPwrParam =    {'delay_time': 2, 'channel':settings['hiLoMWPwr_channel']}
        self.MWswitch3Param =    {'delay_time': 2, 'channel':settings['MWswitch3_channel']}
        self.MWswitch4Param =    {'delay_time': 2, 'channel':settings['MWswitch4_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']

        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam, 
                          'LaserRead': self.LaserReadParam, 'LaserRead2': self.LaserRead2Param, 'LaserInit': self.LaserInitParam, 'LaserIon': self.LaserIonParam,
                        'MW_I': self.MWIParam, 'MW_Q': self.MWQParam, 'MWswitch': self.MWswitchParam,'PB_type': 'USB',
                        'MW_I2': self.MWI2Param, 'MW_Q2': self.MWQ2Param, 'MWswitch2': self.MWswitch2Param, 'hiLoMWPwr': self.hiLoMWPwrParam,
                        'MWswitch3': self.MWswitch3Param, 'MWswitch4': self.MWswitch4Param, 
                        'min_pulse_dur': int(5*1e3/self.clock_speed), 'ifPlotPulse': ifPlotPulse}
        self.settings = {**settings, **settings_extra}
        self.metadata.update(self.settings)

        self.tausArray = self.settings['tausArray']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)

        self.SRSnum = self.settings['SRSnum'];   MWPower = self.settings['MWPower'];   MWFreq = self.settings['MWFreq']
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        self.SRSnum3 = self.settings['SRSnum3']; MWPower3 = self.settings['MWPower3']; MWFreq3 = self.settings['MWFreq3']
        self.SRSnum4 = self.settings['SRSnum4']; MWPower4 = self.settings['MWPower4']; MWFreq4 = self.settings['MWFreq4']
        self.ifMWDuringRead = self.settings['ifMWDuringRead']; self.ifMW2DuringRead = self.settings['ifMW2DuringRead']

        self.ifNeedVel = self.settings['ifNeedVel']; self.velNum = self.settings['velNum']
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
        self.add_parameter(
            name = "sigFullData",
            parameter_class = SigFullData,
        )
        self.add_parameter(
            name = "refFullData",
            parameter_class = RefFullData,
        )
        
        self.savedPulseSequencePlots = {}

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
            if self.ifMW2DuringRead: self.srs3.enable_RFOutput()
            global srs3; srs3 = self.srs3

            self.srs4 = SRS(SRSnum=self.SRSnum4)
            self.srs4.set_freq(MWFreq4) #Hz
            self.srs4.set_RFAmplitude(MWPower4) #dBm
            if (self.SRSnum4 != 3) and (self.SRSnum4 != 4):
                self.srs4.enableIQmodulation()
            if self.ifMW2DuringRead: self.srs4.enable_RFOutput()
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

            for i in range(5):
                self.vel.set_vpiezo(self.vel_vpz_target)
                self.vel.waitUntilComplete()
                self.vel.set_ready()
                time.sleep(0.7)
            global vel; vel = self.vel
        
        # Make Pulse Blaster, Counter, SRS global objects
        global pb
              
    
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref # this is implemented as a Parameter
        sigFullData = self.sigFullData
        refFullData = self.refFullData

        # For each iteration, sweep tau (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            ifIndividualCount = True, 
            numOfPBLoops = self.settings['num_loops'],
            sleepTimeAfterFinishing=0).each(sig,ref, sigFullData, refFullData,
                                            qctask(sig.plotPulseSequences),
                                            )
        data = loop.get_data_set(name='T2ESCCRRIrberDualNVSingleRead')
        data.add_metadata(self.settings)
        self.data = data
        
        loop.run()
        print('Data saved to ' + str(data.location) + '/')
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)
        
        self.srs.disable_RFOutput();  self.srs.disableModulation()
        self.srs2.disable_RFOutput(); self.srs2.disableModulation()
        self.srs3.disable_RFOutput(); self.srs3.disableModulation()
        self.srs4.disable_RFOutput(); self.srs4.disableModulation()

    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/T2ESCCRRIrberDualNVSingleReadObject_sig_set.dat'
    
class Signal(Parameter):
    def __init__(self, settings=None, name='sig', measurementObject=None, **kwargs):
        super().__init__(name, **kwargs)
        self.settings = settings
        self.T2ESCCRRIrberDualNVSingleReadObject = measurementObject
        self.loopCounter = 0
        self.ifMWDuringRead  = self.settings['ifMWDuringRead']
        self.ifMW2DuringRead = self.settings['ifMW2DuringRead']
        self.ifMWReadLowDutyCycle = self.settings['ifMWReadLowDutyCycle']
        self.ifFancySpinInit = self.settings['ifFancySpinInit']
        self.ifAntiCorrel    = self.settings['ifAntiCorrel']
        self.ifSinDetect     = self.settings['ifSinDetect']
        self.ifJustRef_CorrACorr = self.settings['ifJustRef_CorrACorr']
        self.tausArray = self.settings['tausArray']
        self.RRtrackingSettings  = self.settings['RRtrackingSettings']
        self.RRtrackingSettings2 = self.settings['RRtrackingSettings2']

    def set_raw(self, tau_ns):
        NO_MS_EQUALS_1 = 0
        Q_FINAL = 1
        THREE_PI_HALF_FINAL = 2

        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))
        
        # Pulse parameters
        num_loops               = self.settings['num_loops'];             
        laser_init_delay        = self.settings['laser_init_delay'];       laser_init_duration    = self.settings['laser_init_duration']
        laser_to_pi_delay       = self.settings['laser_to_pi_delay'];      pi_time                = self.settings['pi_time']
        pi_to_ion_delay         = self.settings['pi_to_ion_delay'];        ion_duration           = self.settings['ion_duration']
        RRLaserSwitch_delay     = self.settings['RRLaserSwitch_delay'];    DAQ_duration           = self.settings['DAQ_duration']
        DAQ_to_laser_off_delay  = self.settings['DAQ_to_laser_off_delay']; ion_duration2          = self.settings['ion_duration2']
        laserRead_to_MWmix      = self.settings['laserRead_to_MWmix'];     ion_to_read_delay      = self.settings['ion_to_read_delay']
        iznLaserSwitch_delay    = self.settings['iznLaserSwitch_delay'];   RRLaser2Switch_delay   = self.settings['RRLaser2Switch_delay']
        MWmix_duration_short    = self.settings['MWmix_duration_short'];   delay_between_MWmix    = self.settings['delay_between_MWmix']
        nSpinInit               = self.settings['nSpinInit'];              MWI_to_switch_delay    = self.settings['MWI_to_switch_delay']
        spinInit_RR_duration    = self.settings['spinInit_RR_duration'];   spinInit_RR_to_pi_delay= self.settings['spinInit_RR_to_pi_delay']
        pi_time3                = self.settings['pi_time3'];               spinInit_pi_to_RR_delay= self.settings['spinInit_pi_to_RR_delay']
        normalized_style        = self.settings['normalized_style'];       
        
        nMWread   = int(DAQ_duration/(MWmix_duration_short + delay_between_MWmix))
        self.num_loops = num_loops

        # Make pulse sequence
        pulse_sequence = []
        ###################################################################################################################################
        ###################################################################################################################################
        ###################################################################################################################################
        ################## Charge init #######################
        when_init_end           = laser_init_delay + laser_init_duration
        ##################### Spin init ###################### Last MW is to init at ms=0 ############
        spinInit_delay          = when_init_end + laser_to_pi_delay
        spinInit_duration_each  = spinInit_RR_duration + spinInit_RR_to_pi_delay + pi_time3 + spinInit_pi_to_RR_delay
        spinInit_duration_total = nSpinInit*spinInit_duration_each
        spinInit_lastRR_delay   = spinInit_delay + spinInit_duration_total
        spinInit_lastMW_delay   = spinInit_lastRR_delay + spinInit_RR_duration + spinInit_RR_to_pi_delay
        when_spinInit_lastMW_end= spinInit_lastMW_delay + pi_time
        ################## T2E #######################
        firstPiHalf_delay       = when_spinInit_lastMW_end + spinInit_RR_to_pi_delay
        when_firstPiHalf_end    = firstPiHalf_delay + pi_time/2
        pi_delay                = when_firstPiHalf_end + tau_ns/2
        when_pi_end             = pi_delay + pi_time
        secondPiHalf_delay      = when_pi_end + tau_ns/2
        when_secondPiHalf_end   = secondPiHalf_delay + pi_time/2
        ################## Ionization #######################
        ion_RR_delay            = when_secondPiHalf_end + pi_to_ion_delay                 
        when_ion_RR_end         = ion_RR_delay + ion_duration
        ion_strong_delay        = ion_RR_delay + RRLaserSwitch_delay - iznLaserSwitch_delay
        ion_RR2_delay           = ion_RR_delay + RRLaserSwitch_delay - RRLaser2Switch_delay
        ion_strong_duration     = np.max((ion_duration, ion_duration2))
        when_ion_strong_end     = ion_strong_delay + ion_strong_duration
        ################## Read #######################
        laser_read_signal_delay    = when_ion_strong_end + ion_to_read_delay
        DAQ_signal_delay           = laser_read_signal_delay + RRLaserSwitch_delay
        DAQ_signal_duration        = DAQ_duration
        when_DAQ_signal_end        = DAQ_signal_delay + DAQ_signal_duration 
        when_laser_read_signal_end = when_DAQ_signal_end + DAQ_to_laser_off_delay
        laser_read_signal_duration = when_laser_read_signal_end - laser_read_signal_delay

        MWmix_delay                = laser_read_signal_delay + laserRead_to_MWmix
        MWmix_duration             = laser_read_signal_duration - laserRead_to_MWmix + RRLaserSwitch_delay
        when_MWmix_end             = MWmix_delay + MWmix_duration

        when_everything_sig_end = np.max((when_laser_read_signal_end, when_DAQ_signal_end, when_MWmix_end))
        ###################################################################################################################################
        ###################################################################################################################################
        ###################################################################################################################################
        ################## Charge init #######################
        laser_init_ref_delay    = when_everything_sig_end + laser_init_delay
        laser_init_ref_duration = laser_init_duration

        when_init_ref_end       = laser_init_ref_delay + laser_init_ref_duration
        ################## Spin init ####################### Last MW is to init at ms=0 ############
        spinInit_ref_delay           = when_init_ref_end + laser_to_pi_delay
        spinInit_duration_each       = spinInit_RR_duration + spinInit_RR_to_pi_delay + pi_time3 + spinInit_pi_to_RR_delay
        spinInit_duration_total      = nSpinInit*spinInit_duration_each
        spinInit_lastRR_ref_delay    = spinInit_ref_delay + spinInit_duration_total
        spinInit_lastMW_ref_delay    = spinInit_lastRR_ref_delay + spinInit_RR_duration + spinInit_RR_to_pi_delay
        when_spinInit_lastMW_ref_end = spinInit_lastMW_ref_delay + pi_time
        ################## T2E #######################
        firstPiHalf_ref_delay       = when_spinInit_lastMW_ref_end + spinInit_RR_to_pi_delay
        when_firstPiHalf_ref_end     = firstPiHalf_ref_delay + pi_time/2
        pi_ref_delay                 = when_firstPiHalf_ref_end + tau_ns/2
        when_pi_ref_end              = pi_ref_delay + pi_time
        secondPiHalf_ref_delay       = when_pi_ref_end + tau_ns/2
        when_secondPiHalf_ref_end    = secondPiHalf_ref_delay + pi_time/2
        ################## Ionization #######################
        ion_ref_RR_delay        = when_secondPiHalf_ref_end + pi_to_ion_delay
        when_ion_ref_RR_end     = ion_ref_RR_delay + ion_duration
        ion_ref_strong_delay    = ion_ref_RR_delay + RRLaserSwitch_delay - iznLaserSwitch_delay
        ion_ref_RR2_delay       = ion_ref_RR_delay + RRLaserSwitch_delay - RRLaser2Switch_delay
        when_ion_ref_strong_end = ion_ref_strong_delay + ion_strong_duration
        ################## Read #######################
        laser_read_ref_delay    = when_ion_ref_strong_end + ion_to_read_delay
        DAQ_ref_delay           = laser_read_ref_delay + RRLaserSwitch_delay
        DAQ_ref_duration        = DAQ_duration
        when_DAQ_ref_end        = DAQ_ref_delay + DAQ_ref_duration
        when_laser_read_ref_end = when_DAQ_ref_end + DAQ_to_laser_off_delay
        laser_read_ref_duration = when_laser_read_ref_end - laser_read_ref_delay

        MWmix_ref_delay         = laser_read_ref_delay + laserRead_to_MWmix
        MWmix_ref_duration      = laser_read_ref_duration - laserRead_to_MWmix + RRLaserSwitch_delay
        
############################################## Signal ############################################################
##############################################################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_delay,        duration=int(laser_init_duration))] # times are in ns
        
        if self.ifFancySpinInit:
            for i in range(nSpinInit):
                delayRR = spinInit_delay + i*spinInit_duration_each
                delayMW2 = delayRR + spinInit_RR_duration + spinInit_RR_to_pi_delay
                pulse_sequence += [spc.Pulse('LaserRead',  delayRR,            duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch3',  delayMW2,           duration=int(pi_time3))] 
                pulse_sequence += [spc.Pulse('LaserRead2', delayRR,            duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch4',  delayMW2,           duration=int(pi_time3))]
            pulse_sequence += [spc.Pulse('LaserRead',  spinInit_lastRR_delay,  duration=int(spinInit_RR_duration))] 
            pulse_sequence += [spc.Pulse('LaserRead2', spinInit_lastRR_delay,  duration=int(spinInit_RR_duration))] 
         
        if self.ifAntiCorrel == 2:
            pulse_sequence += [spc.Pulse('MWswitch2',  spinInit_lastMW_delay,  duration=int(pi_time))]
        elif self.ifAntiCorrel == 1:
            pulse_sequence += [spc.Pulse('MWswitch',  spinInit_lastMW_delay,  duration=int(pi_time))]

        # T2E: init at ms=-1, ends at ms=-1
        pulse_sequence += [spc.Pulse('MWswitch',      firstPiHalf_delay,       duration=int(pi_time/2))]
        pulse_sequence += [spc.Pulse('MWswitch',      pi_delay,                duration=int(pi_time))]
        pulse_sequence += [spc.Pulse('MWswitch',      secondPiHalf_delay,      duration=int(pi_time/2))]

        pulse_sequence += [spc.Pulse('MWswitch2',     firstPiHalf_delay,       duration=int(pi_time/2))]
        pulse_sequence += [spc.Pulse('MWswitch2',     pi_delay,                duration=int(pi_time))]
        pulse_sequence += [spc.Pulse('MWswitch2',     secondPiHalf_delay,      duration=int(pi_time/2))]
        
        if self.ifSinDetect:
            pulse_sequence += [spc.Pulse('MW_I', secondPiHalf_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]
            pulse_sequence += [spc.Pulse('MW_I2',secondPiHalf_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]

        # Ionization
        pulse_sequence += [spc.Pulse('LaserIon',      ion_strong_delay,        duration=int(ion_strong_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',     ion_RR_delay,            duration=int(ion_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead2',    ion_RR2_delay,           duration=int(ion_duration2))] 
        
        # Read sig NV1
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_signal_delay, duration=int(laser_read_signal_duration))] 
        if self.ifMWDuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch',  delay,          duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch', MWmix_delay,         duration=int(MWmix_duration))]
        if self.ifMW2DuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch3',  delay,         duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch3',MWmix_delay,         duration=int(MWmix_duration))]

        if self.ifMWReadLowDutyCycle == 0:
            pulse_sequence += [spc.Pulse('hiLoMWPwr',    MWmix_delay,         duration=int(MWmix_duration))]

        pulse_sequence += [spc.Pulse('Counter',      DAQ_signal_delay,        duration=int(DAQ_signal_duration))] 
        
        # Read sig NV2
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_signal_delay, duration=int(laser_read_signal_duration))] 
        if self.ifMWDuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch2',  delay,          duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch2', MWmix_delay,         duration=int(MWmix_duration))]
        if self.ifMW2DuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch4',  delay,         duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch4',MWmix_delay,         duration=int(MWmix_duration))]

########################################################### Reference ###################################################################
##############################################################################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',laser_init_ref_delay, duration=int(laser_init_ref_duration))] 
        
        if self.ifFancySpinInit:
            for i in range(nSpinInit):
                delayRR = spinInit_ref_delay + i*spinInit_duration_each
                delayMW2 = delayRR + spinInit_RR_duration + spinInit_RR_to_pi_delay
                pulse_sequence += [spc.Pulse('LaserRead',  delayRR,             duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch3',  delayMW2,            duration=int(pi_time3))]
                pulse_sequence += [spc.Pulse('LaserRead2', delayRR,             duration=int(spinInit_RR_duration))] 
                pulse_sequence += [spc.Pulse('MWswitch4',  delayMW2,            duration=int(pi_time3))] 
            pulse_sequence += [spc.Pulse('LaserRead', spinInit_lastRR_ref_delay,duration=int(spinInit_RR_duration))] 
            pulse_sequence += [spc.Pulse('LaserRead2',spinInit_lastRR_ref_delay,duration=int(spinInit_RR_duration))] 
        
        if self.ifAntiCorrel == 2 or self.ifJustRef_CorrACorr == 2:
            pulse_sequence += [spc.Pulse('MWswitch2', spinInit_lastMW_ref_delay,duration=int(pi_time))]
        if self.ifAntiCorrel == 1 or self.ifJustRef_CorrACorr == 1:
            pulse_sequence += [spc.Pulse('MWswitch', spinInit_lastMW_ref_delay,duration=int(pi_time))]

        ########################## T2E: init at ms=-1, ends at ms=0
        pulse_sequence += [spc.Pulse('MWswitch',      firstPiHalf_ref_delay,       duration=int(pi_time/2))]
        pulse_sequence += [spc.Pulse('MWswitch',      pi_ref_delay,                duration=int(pi_time))]

        pulse_sequence += [spc.Pulse('MWswitch2',     firstPiHalf_ref_delay,       duration=int(pi_time/2))]
        pulse_sequence += [spc.Pulse('MWswitch2',     pi_ref_delay,                duration=int(pi_time))]
        
        if normalized_style == Q_FINAL:
            pulse_sequence += [spc.Pulse('MWswitch',  secondPiHalf_ref_delay,      duration=int(pi_time/2))]
            pulse_sequence += [spc.Pulse('MWswitch2', secondPiHalf_ref_delay,      duration=int(pi_time/2))]
            
            if self.ifJustRef_CorrACorr == 0: # if just measure ref, don't turn IQ on at all
                if self.ifSinDetect:
                    pulse_sequence += [spc.Pulse('MW_Q',  secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]
                    pulse_sequence += [spc.Pulse('MW_Q2', secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]
                else:
                    pulse_sequence += [spc.Pulse('MW_I',  secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]
                    pulse_sequence += [spc.Pulse('MW_Q',  secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))] # times are in ns
                    pulse_sequence += [spc.Pulse('MW_I2', secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]
                    pulse_sequence += [spc.Pulse('MW_Q2', secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))] # times are in ns
            else:
                if self.ifSinDetect:
                    pulse_sequence += [spc.Pulse('MW_I',  secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]
                    pulse_sequence += [spc.Pulse('MW_I2', secondPiHalf_ref_delay-MWI_to_switch_delay, duration=int(pi_time/2 + 2*MWI_to_switch_delay))]
        
        elif normalized_style == THREE_PI_HALF_FINAL:
            pulse_sequence += [spc.Pulse('MWswitch',  secondPiHalf_ref_delay,      duration=int(3*pi_time/2))]
            pulse_sequence += [spc.Pulse('MWswitch2', secondPiHalf_ref_delay,      duration=int(3*pi_time/2))]

        ########################## Ionization
        pulse_sequence += [spc.Pulse('LaserIon',     ion_ref_strong_delay, duration=int(ion_strong_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead',    ion_ref_RR_delay,     duration=int(ion_duration))] 
        pulse_sequence += [spc.Pulse('LaserRead2',   ion_ref_RR2_delay,    duration=int(ion_duration2))]

        ########################## Read ref NV1
        pulse_sequence += [spc.Pulse('LaserRead',    laser_read_ref_delay, duration=int(laser_read_ref_duration))] 
        if self.ifMWDuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_ref_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch',  delay,       duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch', MWmix_ref_delay,  duration=int(MWmix_ref_duration))]
        if self.ifMW2DuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_ref_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch3', delay,       duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch3',MWmix_ref_delay,  duration=int(MWmix_ref_duration))]
        
        if self.ifMWReadLowDutyCycle == 0:
            pulse_sequence += [spc.Pulse('hiLoMWPwr',    MWmix_ref_delay,  duration=int(MWmix_ref_duration))]
        
        pulse_sequence += [spc.Pulse('Counter',      DAQ_ref_delay,        duration=int(DAQ_ref_duration))] 

        ########################## Read ref NV2
        pulse_sequence += [spc.Pulse('LaserRead2',    laser_read_ref_delay, duration=int(laser_read_ref_duration))] 
        if self.ifMWDuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_ref_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch2',  delay,       duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch2', MWmix_ref_delay,  duration=int(MWmix_ref_duration))]
        if self.ifMW2DuringRead:
            if self.ifMWReadLowDutyCycle:
                for i in range(nMWread):
                    delay = MWmix_ref_delay + i*(MWmix_duration_short + delay_between_MWmix)
                    pulse_sequence += [spc.Pulse('MWswitch4', delay,       duration=int(MWmix_duration_short))] 
            else:
                pulse_sequence += [spc.Pulse('MWswitch4',MWmix_ref_delay,  duration=int(MWmix_ref_duration))]
        
       
##############################################################################################################################
##############################################################################################################################

        self.read_duration = DAQ_signal_duration
    
        self.pulse_sequence = pulse_sequence
        self.ifPrintTime = True
        self.pb = spc.B00PulseBlaster("SpinCorePB", settings=self.settings, verbose=False, ifPrintTime=self.ifPrintTime)
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
            samps_per_chan = 5*self.num_reads # x2 to make sure buffer doesn't overflow
            )
        pulseWidthChan.ci_ctr_timebase_src = "/cDAQ1Mod1/PFI0" # counter out PFI str gated/counter PFI channel str
        pulseWidthChan.ci_pulse_width_term = "/cDAQ1Mod1/PFI1" # gate PFI string

        print('Set tau to ' + str(tau_ns) + " ns")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))
        self.ctrtask.stop(); self.ctrtask.close()
        rate = xLineData

        global sig_data;  sig_data  = rate[::self.num_reads_per_iter]
        global ref_data;  ref_data  = rate[1::self.num_reads_per_iter]
        global sig_avg;   sig_avg  = np.average(sig_data)/(self.read_duration/1e6)
        global ref_avg;   ref_avg  = np.average(ref_data)/(self.read_duration/1e6)

        # Line tracking
        if self.settings['laserRead_channel'] == 5 or self.settings['laserRead_channel'] == 14:
            threshold_scanVpz = self.RRtrackingSettings['threshold_scanVpz']
            #######################################################
            # NV 1
            if self.RRtrackingSettings['if_tracking'] == 1 and np.max((np.max(sig_avg), np.max(ref_avg))) < threshold_scanVpz:
                print()
                time_sleep_after_scan = self.RRtrackingSettings['time_sleep_after_scan']
                wvl_correction = self.RRtrackingSettings['wvl_correction']
                scan_lower_margin = self.RRtrackingSettings['scan_lower_margin']
                scan_upper_margin = self.RRtrackingSettings['scan_upper_margin']
                num_point_scan = self.RRtrackingSettings['num_point_scan']
                print('-----------------Start line tracking for NV1---------------------------')
                lockStatusFile = 'C:\\Users\\lukin2dmaterials\\miniconda3\\envs\\NV_control\\B00_codes\\laser1_lockStatus.txt'
                
                # Turn off laser lock
                try:
                    with open(lockStatusFile, 'r+') as file:
                        file.seek(0)
                        file.write(str(0))
                except OSError as e:
                    print(f"Error: {e}")
                    try:
                        with open(lockStatusFile, 'r+') as file:
                            file.seek(0)
                            file.write(str(0))
                    except OSError as e:
                        print(f"Error: {e}")
                
                # Scan vpz
                lockStatusParam = np.loadtxt(lockStatusFile)
                vpz_old = lockStatusParam[2]
                start = vpz_old - scan_lower_margin; stop = vpz_old + scan_upper_margin; num_sweep_points = num_point_scan
                vpzArray = np.linspace(start, stop, num_sweep_points)
                self.RRtrackingSettings['vpzArray'] = vpzArray
                
                ScanRRFreqObject = ScanRRFreq(settings=self.RRtrackingSettings, ifPlotPulse=0)
                vpz_new, wvl = ScanRRFreqObject.runScanInPulseSequenceMonitoringWM()

                # Turn on laser lock again
                with open(lockStatusFile, 'r+') as file:
                    wvl = wvl + wvl_correction
                    file.seek(0)
                    file.write(str(1)+'\n')
                    file.write(str(wvl)+'\n')
                    file.write(str(vpz_new+0.01))
                    file.truncate()
                
                print("Lock laser 1 to " + str(wvl) + ' THz') # set the Velocity's piezo voltage
                print('-----------------End line tracking for NV1---------------------------')
                print()
                
            #######################################################    
            # NV 2
                print()
                time_sleep_after_scan = self.RRtrackingSettings2['time_sleep_after_scan']
                wvl_correction = self.RRtrackingSettings2['wvl_correction']
                scan_lower_margin = self.RRtrackingSettings2['scan_lower_margin']
                scan_upper_margin = self.RRtrackingSettings2['scan_upper_margin']
                num_point_scan = self.RRtrackingSettings2['num_point_scan']
                print('-----------------Start line tracking for NV2---------------------------')
                lockStatusFile = 'C:\\Users\\lukin2dmaterials\\miniconda3\\envs\\NV_control\\B00_codes\\laser2_lockStatus.txt'
                
                # Turn off laser lock
                try:
                    with open(lockStatusFile, 'r+') as file:
                        file.seek(0)
                        file.write(str(0))
                except OSError as e:
                    print(f"Error: {e}")
                    try:
                        with open(lockStatusFile, 'r+') as file:
                            file.seek(0)
                            file.write(str(0))
                    except OSError as e:
                        print(f"Error: {e}")
                
                # Scan vpz
                lockStatusParam = np.loadtxt(lockStatusFile)
                vpz_old = lockStatusParam[2]
                start = vpz_old - scan_lower_margin; stop = vpz_old + scan_upper_margin; num_sweep_points = num_point_scan
                vpzArray = np.linspace(start, stop, num_sweep_points)
                self.RRtrackingSettings2['vpzArray'] = vpzArray
                
                ScanRRFreqObject = ScanRRFreq(settings=self.RRtrackingSettings2, ifPlotPulse=0)
                vpz_new, wvl = ScanRRFreqObject.runScanInPulseSequenceMonitoringWM()

                # Turn on laser lock again
                with open(lockStatusFile, 'r+') as file:
                    wvl = wvl + wvl_correction
                    file.seek(0)
                    file.write(str(1)+'\n')
                    file.write(str(wvl)+'\n')
                    file.write(str(vpz_new + 0.01))
                    file.truncate()
                
                print("Lock laser 2 to " + str(wvl) + ' THz') # set the Velocity's piezo voltage
                print('-----------------End line tracking for NV2, wait for ' + str(time_sleep_after_scan) + ' s---------------------------')
                print()
                time.sleep(time_sleep_after_scan)

        return sig_avg
    
    def plotPulseSequences(self):
        self.readColor = self.settings['LaserRead']['channel']
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False, readColor=self.readColor)
                fig = plotPulseObject.makePulsePlot()
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.T2ESCCRRIrberDualNVSingleReadObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1
    
    def turn_on_at_end(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        pb.turn_on_infinite(channels=channels)


class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg

class SigFullData(Parameter):
    def __init__(self, name='sigFullData',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_data
    
class RefFullData(Parameter):
    def __init__(self, name='refFullData',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_data





