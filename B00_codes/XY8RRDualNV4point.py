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

class XY8RRDualNV4point(Instrument):

    def __init__(self, name='XY8RRDualNV4pointObject', settings=None, ifPlotPulse=True, **kwargs) -> None:

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
        self.AWG2Param =         {'delay_time': 2, 'channel':settings['AWG2_channel']}
        self.hiLoMWPwrParam =    {'delay_time': 2, 'channel':settings['hiLoMWPwr_channel']}
        global laserInitChannel; laserInitChannel = self.LaserInitParam['channel']
    
        settings_extra = {'clock_speed': self.clock_speed, 'Counter': self.CounterParam,
                          'LaserRead': self.LaserReadParam, 'LaserInit': self.LaserInitParam,'hiLoMWPwr': self.hiLoMWPwrParam,
                          'AWG': self.AWGParam, 'AWG2': self.AWG2Param,
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
        self.tausArray = self.settings['tausArray']
        self.SRSnum = self.settings['SRSnum']; MWPower = self.settings['MWPower']; MWFreq = self.settings['MWFreq']
        self.SDGnum = self.settings['SDGnum']; self.srate=self.settings['srate']
        self.ifNeedVel1 = self.settings['ifNeedVel1']; self.velNum = self.settings['velNum']
        vel_current = self.settings['vel_current']; vel_wvl = self.settings['vel_wvl']; 
        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.ifInitVpz = self.settings['ifInitVpz']; self.ifInitWvl = self.settings['ifInitWvl']

        # Microwave and velocity 2
        self.SRSnum2 = self.settings['SRSnum2']; MWPower2 = self.settings['MWPower2']; MWFreq2 = self.settings['MWFreq2']
        self.SDGnum2 = self.settings['SDGnum2']
        self.ifNeedVel2 = self.settings['ifNeedVel2']; self.velNum2 = self.settings['velNum2']
        vel_current2 = self.settings['vel_current2']; vel_wvl2 = self.settings['vel_wvl2']; 
        self.vel_vpz_target2 = self.settings['vel_vpz_target2']

        ifRandomized = self.settings['ifRandomized']
        if ifRandomized: np.random.shuffle(self.tausArray)
        ###########################################################################################

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
            name = "refx",
            parameter_class = ReferenceX,
        )
        self.add_parameter(
            name = "refxm",
            parameter_class = ReferenceXm,
        )
        self.add_parameter(
            name = "sig2",
            parameter_class = Signal2,
        )
        self.add_parameter(
            name = "ref2",
            parameter_class = Reference2,
        )
        self.add_parameter(
            name = "refx2",
            parameter_class = ReferenceX2,
        )
        self.add_parameter(
            name = "refxm2",
            parameter_class = ReferenceXm2,
        )
        
        self.savedPulseSequencePlots = {}

        # SRS object 1
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

        # AWG object 1
        self.AWG = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum, srate=self.srate)
        global AWG; AWG = self.AWG

        # AWG object 2
        self.AWG2 = SDG6022X(name='SDG6022X', SDGnum=self.SDGnum2, srate=self.srate)
        global AWG2; AWG2 = self.AWG2

        # Velocity object 1
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
        global srs; srs = self.srs; global srs2; srs2 = self.srs2
        
    def runScan(self):
        sig = self.sig # this is implemented as a Parameter
        ref = self.ref
        refx = self.refx
        refxm = self.refxm
        sig2 = self.sig2
        ref2 = self.ref2
        refx2 = self.refx2
        refxm2 = self.refxm2

        # For each iteration, sweep frequency (sig.sweep calls set_raw() method of Parameter sig)
        # and measure Parameter sig, ref (each(sig,ref)) by calling get_raw() method of sig, ref
        loop = Loop(
            sig.sweep(keys=self.tausArray),
            delay = 0,
            sleepTimeAfterFinishing=0).each(sig,sig2,ref,ref2,refx,refx2,refxm,refxm2,
                                            qctask(sig.plotPulseSequences),
                                            )

        data = loop.get_data_set(name='XY8RRDualNV4point')
        data.add_metadata(self.settings)
        self.data = data
        
        plot = QtPlot(
            data.XY8RRDualNV4pointObject_sig, # this is implemented as a Parameter
            figsize = (1200, 600),
            interval = 1,
            name = 'sig'
            )
        plot.add(data.XY8RRDualNV4pointObject_sig2, name='sig2')
        plot.add(data.XY8RRDualNV4pointObject_ref, name='ref')
        plot.add(data.XY8RRDualNV4pointObject_ref2, name='ref2')

        loop.with_bg_task(plot.update, bg_final_task=None)
        loop.run()
        print('Data saved to ' + str(data.location) + '/')

        dataPlotFilename = data.location + "/dataPlot.png"
        dataPlotFile = plot.save(filename=dataPlotFilename, type='data')
        # img = Image.open(dataPlotFile)
        # img.show()
        
        if self.settings['ifPlotPulse']: # save the first and last pulse sequence plot
            for index in self.savedPulseSequencePlots:
                fig = self.savedPulseSequencePlots[index]
                pulsePlotFilename = data.location + "/pulsePlot_" + str(index) + ".png"
                fig.savefig(pulsePlotFilename)
        
        self.srs.disable_RFOutput()
        self.srs.disableModulation()
        self.srs2.disable_RFOutput()
        self.srs2.disableModulation()
        self.AWG.turn_off()
        self.AWG2.turn_off()
    
    def getDataFilename(self):
        return 'C:/Users/lukin2dmaterials/' + self.data.location + '/XY8RRDualNV4pointObject_sig_set.dat'

    
class Signal(Parameter):
    def __init__(self, name='sig', measurementObject=None, settings=None, **kwargs):
        super().__init__(name, **kwargs)
        self.loopCounter = 0
        self.settings = settings
        self.tausArray = self.settings['tausArray']
        self.XY8RRDualNV4pointObject = measurementObject
        self.RRtrackingSettings = self.settings['RRtrackingSettings']
        self.RRtrackingSettings2 = self.settings['RRtrackingSettings2']
        self.ifHiloExtra = self.settings['ifHiloExtra']
        self.srate = self.settings['srate']

        self.readColor    = self.settings['LaserRead']['channel']
        self.initColor    = self.settings['LaserInit']['channel']

        self.vel_vpz_target = self.settings['vel_vpz_target']
        self.vel_vpz_target2 = self.settings['vel_vpz_target2']

    def get_raw(self):
        self.ctrtask.start()

        self.pb.start_pulse_seq()
        self.pb.wait()
        self.pb.stop_pulse_seq(); self.pb.close()

        xLineData = np.array(self.ctrtask.read(self.num_reads, timeout=0.5))

        self.ctrtask.stop(); self.ctrtask.close()

        rate = xLineData
        sig = rate[::self.num_reads_per_iter]
        sig2 = rate[1::self.num_reads_per_iter]
        ref = rate[2::self.num_reads_per_iter]
        ref2 = rate[3::self.num_reads_per_iter]
        refx = rate[4::self.num_reads_per_iter]
        refx2 = rate[5::self.num_reads_per_iter]
        refxm = rate[6::self.num_reads_per_iter]
        refxm2 = rate[7::self.num_reads_per_iter]
        global sig_avg;  sig_avg = np.average(sig)/(self.read_duration/1e9)/1e3
        global ref_avg;  ref_avg = np.average(ref)/(self.read_duration/1e9)/1e3
        global refx_avg;  refx_avg = np.average(refx)/(self.read_duration/1e9)/1e3
        global refxm_avg;  refxm_avg = np.average(refxm)/(self.read_duration/1e9)/1e3
        global sig_avg2;  sig_avg2 = np.average(sig2)/(self.read_duration2/1e9)/1e3
        global ref_avg2;  ref_avg2 = np.average(ref2)/(self.read_duration2/1e9)/1e3
        global refx_avg2;  refx_avg2 = np.average(refx2)/(self.read_duration/1e9)/1e3
        global refxm_avg2;  refxm_avg2 = np.average(refxm2)/(self.read_duration/1e9)/1e3
   
        return sig_avg
    
    def set_raw(self, tau_ns):
        # Make pulses, program Pulse Blaster
        print("Loop " + str(self.loopCounter))

        NO_MS_EQUALS_1 = 0
        Q_FINAL = 1
        THREE_PI_HALF_FINAL = 2

        # Pulse parameters
        num_loops               = self.settings['num_loops']
        laser_init_delay        = self.settings['laser_init_delay'];    laser_init_duration = self.settings['laser_init_duration']
        laser_to_MWI_delay      = self.settings['laser_to_MWI_delay'];  laser_to_DAQ_delay  = self.settings['laser_to_DAQ_delay']
        pi_half                 = self.settings['pi_half'];             pi_half2            = self.settings['pi_half2']
        read_duration           = self.settings['read_duration'];       read_laser_duration = self.settings['read_laser_duration']
        read_duration2          = self.settings['read_duration2'];      read_laser_duration2= self.settings['read_laser_duration2']
        MW_to_read_delay        = self.settings['MW_to_read_delay']
        shift_btwn_2NV_MW       = self.settings['shift_btwn_2NV_MW'];   shift_btwn_2NV_read = self.settings['shift_btwn_2NV_read']
        laser_to_DAQ_delay2     = self.settings['laser_to_DAQ_delay2']; normalized_style    = self.settings['normalized_style']
        AWG_buffer              = self.settings['AWG_buffer'];          AWG_output_delay    = self.settings['AWG_output_delay']  
        AWG_buffer2             = self.settings['AWG_buffer2'];         AWG_output_delay2   = self.settings['AWG_output_delay2']         
        phi_IQ                  = self.settings['phi_IQ'];              phi_IQ2             = self.settings['phi_IQ2']
        pi                      = self.settings['pi_pulse'];            pi2                 = self.settings['pi_pulse2']                                                           
        ifStartInY          = self.settings['ifStartInY']
        NXY81                   = self.settings['NXY81'];               NXY82               = self.settings['NXY82']
        last_pi_2time           = self.settings['last_pi_2time'];       
        tauExtra                = self.settings['tauExtra'];            sweepWhich          = self.settings['sweepWhich']
        if sweepWhich=='N':
            NXY81 = NXY82 = int(tau_ns)
            tau_ns              = self.settings['tau']
        elif sweepWhich=='phase_last_pi_half':
            phi_IQ              = tau_ns; phi_IQ2 = phi_IQ
            tau_ns              = self.settings['tau']
        elif sweepWhich=='laser_init_delay':
            laser_init_delay    = tau_ns
            tau_ns              = self.settings['tau']

        MW_duration_for_AWG     = int(2*int((AWG_buffer  + pi_half  + last_pi_2time + NXY81*tau_ns*8 + NXY81*pi*8  + tauExtra - pi_half)/2))
        MW_duration_for_AWG2    = int(2*int((AWG_buffer2 + pi_half2 + last_pi_2time + NXY82*tau_ns*8 + NXY82*pi2*8 + tauExtra - pi_half2)/2))
        # the last -pi_half is bc the first and last delay is (tau-pihalf)/2
        if tau_ns/2 > 30:
            MWI_to_switch_delay  = self.settings['MWI_to_switch_delay']
        else: 
            MWI_to_switch_delay  = 0

        when_init_end              = laser_init_delay + laser_init_duration

        ################# Signal  ################################################################
        ###### NV1 ######################################
        MWI_delay                  = when_init_end + laser_to_MWI_delay # for AWG       
        global MW_del;      MW_del = MWI_delay + AWG_output_delay # true time when pulse appears

        when_pulse_end = MW_del + MW_duration_for_AWG # true time when pulse ends
        
        ###### NV2 ######################################
        MWI_delay2                 = MWI_delay + shift_btwn_2NV_MW # for AWG2
        global MW_del2;    MW_del2 = MWI_delay2 + AWG_output_delay2 # true time when pulse appears

        when_pulse_end2 = MW_del2 + MW_duration_for_AWG2
        
        laser_read_signal_delay    = np.max((when_pulse_end,when_pulse_end2)) + MW_to_read_delay
        laser_read_signal_duration = read_laser_duration
        when_laser_read_signal_end = laser_read_signal_delay + laser_read_signal_duration
        read_signal_delay          = laser_read_signal_delay + laser_to_DAQ_delay
        read_signal_duration       = read_duration
        when_read_signal_end       = read_signal_delay + read_signal_duration
        
        laser_read_signal_delay2    = laser_read_signal_delay + shift_btwn_2NV_read
        laser_read_signal_duration2 = read_laser_duration2
        when_laser_read_signal_end2 = laser_read_signal_delay2 + laser_read_signal_duration2
        read_signal_delay2          = laser_read_signal_delay2 + laser_to_DAQ_delay2
        read_signal_duration2       = read_duration2
        when_read_signal_end2       = read_signal_delay2 + read_signal_duration2

        ################# Reference  ################################################################
        ###### NV1 ######################################
        laser_init_ref_delay        = np.max((when_read_signal_end2,when_read_signal_end)) + laser_init_delay
        when_init_ref_end           = laser_init_ref_delay + laser_init_duration

        MW_ref_delays     = []
        MW_ref_delays.append(when_init_ref_end + laser_to_MWI_delay)         
        MW_ref_delays.append(MW_ref_delays[0] + pi_half + (tau_ns-pi_half)/2)
        for i in range(2,NXY81*8+1):
            MW_ref_delays.append(MW_ref_delays[i-1] + pi + tau_ns)
        MW_ref_delays.append(MW_ref_delays[NXY81*8] + pi + (tau_ns-pi_half)/2 + tauExtra)

        laser_read_ref_delay        = np.max((when_read_signal_end2,when_read_signal_end)) + laser_read_signal_delay
        laser_read_ref_duration     = read_laser_duration
        when_laser_read_ref_end     = laser_read_ref_delay + laser_read_ref_duration
        read_ref_delay              = laser_read_ref_delay + laser_to_DAQ_delay;  
        read_ref_duration           = read_duration
        when_read_ref_end           = read_ref_delay + read_ref_duration

        sig_to_ref_wait             = (MW_ref_delays[0] + AWG_output_delay) - when_pulse_end
        # MW_ref_delays[0] is equivalent to the delay for AWG triggering pulse
        

        ###### NV2 ######################################
        MW_ref_delays2     = []
        MW_ref_delays2.append(when_init_ref_end + laser_to_MWI_delay + shift_btwn_2NV_MW)         
        MW_ref_delays2.append(MW_ref_delays2[0] + pi_half2 + (tau_ns-pi_half2)/2)
        for i in range(2,NXY82*8+1):
            MW_ref_delays2.append(MW_ref_delays2[i-1] + pi2 + tau_ns)
        MW_ref_delays2.append(MW_ref_delays2[NXY82*8] + pi2 + (tau_ns-pi_half2)/2 + tauExtra)

        laser_read_ref_delay2    = np.max((when_read_signal_end2,when_read_signal_end)) + laser_read_signal_delay2
        laser_read_ref_duration2 = read_laser_duration2
        when_laser_read_ref_end2 = laser_read_ref_delay2 + laser_read_ref_duration2
        read_ref_delay2          = laser_read_ref_delay2 + laser_to_DAQ_delay2
        read_ref_duration2       = read_duration2
        when_read_ref_end2       = read_ref_delay2 + read_ref_duration2

        sig_to_ref_wait2         = (MW_ref_delays2[0] + AWG_output_delay) - when_pulse_end2


        ################# ReferenceX  ################################################################
        ###### NV1 ######################################
        laser_init_refx_delay        = np.max((when_read_ref_end2,when_read_ref_end)) + laser_init_delay
        when_init_refx_end           = laser_init_ref_delay + laser_init_duration

        laser_read_refx_delay        = np.max((when_read_ref_end2,when_read_ref_end)) + laser_read_signal_delay
        laser_read_refx_duration     = read_laser_duration
        when_laser_read_refx_end     = laser_read_refx_delay + laser_read_refx_duration
        read_refx_delay              = laser_read_refx_delay + laser_to_DAQ_delay;  
        read_refx_duration           = read_duration
        when_read_refx_end           = read_refx_delay + read_refx_duration
        ###### NV2 ######################################
        laser_read_refx_delay2    = np.max((when_read_ref_end2,when_read_ref_end)) + laser_read_signal_delay2
        laser_read_refx_duration2 = read_laser_duration2
        when_laser_read_refx_end2 = laser_read_refx_delay2 + laser_read_refx_duration2
        read_refx_delay2          = laser_read_refx_delay2 + laser_to_DAQ_delay2
        read_refx_duration2       = read_duration2
        when_read_refx_end2       = read_refx_delay2 + read_refx_duration2


        ################# ReferenceXm  ################################################################
        ###### NV1 ######################################
        laser_init_refxm_delay        = np.max((when_read_refx_end2,when_read_refx_end)) + laser_init_delay
        when_init_refxm_end           = laser_init_refxm_delay + laser_init_duration

        laser_read_refxm_delay        = np.max((when_read_refx_end2,when_read_refx_end)) + laser_read_signal_delay
        laser_read_refxm_duration     = read_laser_duration
        when_laser_read_refxm_end     = laser_read_refxm_delay + laser_read_refxm_duration
        read_refxm_delay              = laser_read_refxm_delay + laser_to_DAQ_delay;  
        read_refxm_duration           = read_duration
        when_read_refxm_end           = read_refxm_delay + read_refxm_duration
        ###### NV2 ######################################
        laser_read_refxm_delay2    = np.max((when_read_refx_end2,when_read_refx_end)) + laser_read_signal_delay2
        laser_read_refxm_duration2 = read_laser_duration2
        when_laser_read_refxm_end2 = laser_read_refxm_delay2 + laser_read_refxm_duration2
        read_refxm_delay2          = laser_read_refxm_delay2 + laser_to_DAQ_delay2
        read_refxm_duration2       = read_duration2
        when_read_refxm_end2       = read_refxm_delay2 + read_refxm_duration2


        self.read_duration = read_signal_duration
        self.read_duration2 = read_signal_duration2
        if read_signal_duration != read_ref_duration:
            raise Exception("Duration of reading signal and reference must be the same")  
        
        #  Set I and Q for special pi/2 pulse
        I = round(32767*np.cos(phi_IQ)) # phi_IQ in radian = 0 means x-axis
        Q = round(32767*np.sin(phi_IQ))
        I2 = round(32767*np.cos(phi_IQ2))
        Q2 = round(32767*np.sin(phi_IQ2))
        
        a0 = np.array((2,4,5,7)); MWI_list1 = a0
        if NXY81 > 1:
            for i in range(1,NXY81):
                MWI_list1 = np.concatenate((MWI_list1,a0+i*8))
        if ifStartInY:
            MWI_list1 = np.insert(MWI_list1,0,values=0)
        a0 = np.array((2,4,5,7)); MWI_list2 = a0
        if NXY82 > 1:
            for i in range(1,NXY82):
                MWI_list2 = np.concatenate((MWI_list2,a0+i*8))
        if ifStartInY:
            MWI_list2 = np.insert(MWI_list2,0,values=0)
        
        if self.srate is not None:
            if self.loopCounter==0: sleepTime = 10
            else: sleepTime = 1
        else:
            sleepTime = 0

        global ch1plot; global ch2plot
        ch1plot, ch2plot = AWG.send_XY84point_seq(numxy8=NXY81, pi_2time=int(pi_half), 
                                            pitime=int(pi), tau = int(tau_ns), buffer=int(AWG_buffer), sig_to_ref_wait=int(sig_to_ref_wait),
                                            special_I=I, special_Q=Q,
                                            last_pi_2time=int(last_pi_2time),
                                            tauExtra=tauExtra, sleepTime=1)
        global ch1plot2; global ch2plot2
        ch1plot2, ch2plot2 = AWG2.send_XY84point_seq(numxy8=NXY82, pi_2time=int(pi_half2), 
                                            pitime=int(pi2), tau = int(tau_ns), buffer=int(AWG_buffer2), sig_to_ref_wait=int(sig_to_ref_wait2),
                                            special_I=I2, special_Q=Q2,
                                            last_pi_2time=int(last_pi_2time),
                                            tauExtra=tauExtra, sleepTime=sleepTime)
        
        ####################################### Signal #########################################################
        # Make pulse sequence (per each freq)
        pulse_sequence = []
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_delay,             duration=int(laser_init_duration))] # times are in ns

        ########################### Make XY8 sequence #############################
        pulse_sequence += [spc.Pulse('AWG',      MWI_delay,      duration=50)]
        pulse_sequence += [spc.Pulse('AWG2',     MWI_delay2,     duration=50)]
        
        ########################### Read #############################
        pulse_sequence += [spc.Pulse('LaserRead',  laser_read_signal_delay,  duration=int(laser_read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',    read_signal_delay,        duration=int(read_signal_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_signal_delay2, duration=int(laser_read_signal_duration2))] # times are in ns
        pulse_sequence += [spc.Pulse('Counter',    read_signal_delay2,       duration=int(read_signal_duration2))] 
        

        ####################################### Reference ################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_ref_delay,         duration=int(laser_init_duration))] # times are in ns

        ########################### Read #############################
        pulse_sequence += [spc.Pulse('LaserRead',  laser_read_ref_delay,     duration=int(laser_read_ref_duration))]
        pulse_sequence += [spc.Pulse('Counter',    read_ref_delay,           duration=int(read_ref_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_ref_delay2,    duration=int(laser_read_ref_duration2))]
        pulse_sequence += [spc.Pulse('Counter',    read_ref_delay2,          duration=int(read_ref_duration2))] 


        ####################################### ReferenceX ################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_refx_delay,         duration=int(laser_init_duration))] # times are in ns

        ########################### Read #############################
        pulse_sequence += [spc.Pulse('LaserRead',  laser_read_refx_delay,     duration=int(laser_read_refx_duration))]
        pulse_sequence += [spc.Pulse('Counter',    read_refx_delay,           duration=int(read_refx_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_refx_delay2,    duration=int(laser_read_refx_duration2))]
        pulse_sequence += [spc.Pulse('Counter',    read_refx_delay2,          duration=int(read_refx_duration2))]


        ####################################### ReferenceXm ################################################################
        if not laser_init_delay == 0:
            pulse_sequence += [spc.Pulse('LaserInit',    laser_init_refxm_delay,         duration=int(laser_init_duration))] # times are in ns

        ########################### Read #############################
        pulse_sequence += [spc.Pulse('LaserRead',  laser_read_refxm_delay,     duration=int(laser_read_refxm_duration))]
        pulse_sequence += [spc.Pulse('Counter',    read_refxm_delay,           duration=int(read_refxm_duration))] # times are in ns
        pulse_sequence += [spc.Pulse('LaserRead2', laser_read_refxm_delay2,    duration=int(laser_read_refxm_duration2))]
        pulse_sequence += [spc.Pulse('Counter',    read_refxm_delay2,          duration=int(read_refxm_duration2))]  

        ########################### Make hilo pulse ###########################
        if self.ifHiloExtra==1:
            plotPulseObject = PlotPulse(pulseSequence=pulse_sequence, ifSave=False)
            _, ch11, ch21 = plotPulseObject.makeTraceAWG(ch1plot, ch2plot, MW_del)
            _, ch12, ch22 = plotPulseObject.makeTraceAWG(ch1plot2, ch2plot2, MW_del2)
            all_ch = ch11 + ch21 + ch12 + ch22
            zeroSegments = plotPulseObject.find_zero_segments(all_ch)
            for start, length in zeroSegments:
                margin_start = self.settings['hilo_margin_start']; margin_end = self.settings['hilo_margin_end']
                hilo_delay = start + margin_start
                hilo_duration = length - (margin_start + margin_end)
                if hilo_duration >= self.settings['hilo_min']:
                    pulse_sequence += [spc.Pulse('hiLoMWPwr', int(hilo_delay), duration=int(hilo_duration))]
                    
        #####################################################################################################
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

        if sweepWhich=='phase_last_pi_half':
            print('Set phi to ' + str(phi_IQ/np.pi))
        elif sweepWhich=='N':
            print('Set N to ' + str((NXY81+NXY82)/2))
        else:
            print('Set tau to ' + str(tau_ns) + " ns")
        if not self.settings['ifPlotPulse']: self.loopCounter += 1
    
    def plotPulseSequences(self):
        if self.settings['ifPlotPulse']:
            if np.mod(self.loopCounter,5) == 0 or self.loopCounter == len(self.tausArray)-1: # plot every 5 sequences and plot the last seq
                plotPulseObject = PlotPulse(pulseSequence=self.pulse_sequence, ifShown=True, ifSave=False,
                                            readColor=self.readColor, 
                                            initColor=self.initColor)
                fig = plotPulseObject.makePulsePlotAWG(ch1plot, ch2plot, MW_del, label1='ch1-AWG1', label2='ch2-AWG1')
                fig = plotPulseObject.makePulsePlotAWG(ch1plot2,ch2plot2, MW_del2, fig=fig,
                                                           offset1=17.75, offset2=17.5, label1='ch1-AWG2', label2='ch2-AWG2')
                
            if self.loopCounter == 0 or self.loopCounter == len(self.tausArray)-1: # only save first and last pulse sequence
                self.XY8RRDualNV4pointObject.savedPulseSequencePlots[self.loopCounter] = deepcopy(fig)
            self.loopCounter += 1

    def close_turnOnAtEnd(self):
        pb = spc.B00PulseBlaster("SpinCorePBFinal", settings=self.settings, verbose=False)
        channels = np.linspace(laserInitChannel,laserInitChannel,1)
        # pb.turn_on_infinite(channels=channels)
    
class Reference(Parameter):
    def __init__(self, name='ref',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg

class ReferenceX(Parameter):
    def __init__(self, name='refx',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return refx_avg

class ReferenceXm(Parameter):
    def __init__(self, name='refxm',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return refxm_avg

class Signal2(Parameter):
    def __init__(self, name='sig2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return sig_avg2

class Reference2(Parameter):
    def __init__(self, name='ref2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return ref_avg2

class ReferenceX2(Parameter):
    def __init__(self, name='refx2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return refx_avg2

class ReferenceXm2(Parameter):
    def __init__(self, name='refxm2',**kwargs):
        super().__init__(name, **kwargs)

    def get_raw(self):
        return refxm_avg2
    
