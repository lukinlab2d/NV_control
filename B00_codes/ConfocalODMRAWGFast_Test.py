"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.ConfocalODMRAWGFast import *

def main(N=8e5):
    ####################################################################################################################
    # ConfocalODMRAWGFast - Pulsed ODMR
    for i in range(1):
        # xArray = np.linspace(-0.52, -0.44, 3); yArray = np.linspace(-0.40,-0.48,3); ifDifferential = 0
        xArray = np.linspace(-0.52, -0.44, 3); yArray = np.linspace(-1.04,-1.12,3); ifDifferential = 0
        # xArray = np.linspace(-0.36,0.44,21); yArray = np.linspace(-1.48, -0.40, 28); ifDifferential = 0
        if True:
            a11 = 2694e6; a1 = a11-16e6; b1 = a1+32e6; num_sweep_points = 9
            f1 = np.linspace(a1, b1, num_sweep_points)
            a22 = 3061e6; a2 = a22-16e6; b2 = a2+32e6; num_sweep_points = 9
            f2 = np.linspace(a2, b2, num_sweep_points)
            f3 = np.linspace(a1-28e6,a1-8e6,  3); f4 = np.linspace(b2+8e6,  b2+28e6,3)
            f5 = np.linspace(b1+8e6 ,b1+28e6, 3); f6 = np.linspace(a2-28e6, a2-8e6, 3)
            freqsArray = np.concatenate((f3,f1,f5,f6,f2,f4))
        if False:
            f0 = 596e6; a1 = f0-20e6; b1 = f0+20e6; num_sweep_points = 11
            f1 = np.linspace(a1, b1, num_sweep_points)
            f3 = np.linspace(a1-45e6,a1-9e6,  5); f5 = np.linspace(b1+9e6 ,b1+45e6, 5)
            freqsArray = np.concatenate((f3,f1,f5))

        SDGnum=1; SRSnum=1; MWPower = -24.2
        laserInit_channel = 3; laserRead_channel = 3; AWG_channel = 18

        num_loops                    = int(N);   settleTime             = 2e9
        laser_init_delay             = 0;        laser_init_duration    = 0
        pitime                       = 24;       laser_to_AWG_delay     = 0; pitime2=0
        laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900}
        laser_to_DAQ_delay           = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)   
        AWG_output_delay             = 1450;     AWG_buffer             = 10
        read_duration                = 300;      wait_btwn_sig_ref      = 2000
        xref_pitime                  = xArray[0];pi_incr_factor         = 6/0.12; pi_incr_factor2=0 #ns/Vx
        padding                      = 500;      MW_to_DAQ_delay        = 100
        padding_green1               = 200;      AWG_delay              = AWG_output_delay#??
        DAQ_error_factor             = 0.99;     freqsArray2 = None

        if ifDifferential:
            xArray2 = -xArray
            interleaved_array = np.empty(xArray.size + xArray2.size, dtype=xArray.dtype)
            interleaved_array[0::2] = xArray; interleaved_array[1::2] = xArray2
            xArray = interleaved_array

        settings = {'num_loops':num_loops, 'freqsArray':freqsArray, 'freqsArray2':freqsArray2, 
                    'SRSnum':SRSnum, 'MWPower':MWPower, 'SDGnum':SDGnum,
                    'laser_init_delay':       laser_init_delay,      'laser_init_duration': laser_init_duration,
                    'laser_to_AWG_delay':     laser_to_AWG_delay,    'pitime':         pitime, 'pitime2':pitime2,
                    'laser_to_DAQ_delay':     laser_to_DAQ_delay,    'read_duration':       read_duration,
                    'AWG_output_delay':       AWG_output_delay,      'AWG_buffer': AWG_buffer,
                    'laserInit_channel':      laserInit_channel,     'laserRead_channel':   laserRead_channel, 
                    'AWG_channel':    AWG_channel, 'AWG_delay': AWG_delay, 'MW_to_DAQ_delay': MW_to_DAQ_delay,
                    'wait_btwn_sig_ref': wait_btwn_sig_ref,
                    'xArray': xArray, 'yArray':yArray,'settleTime': settleTime,
                    'xref_pitime': xref_pitime, 'pi_incr_factor':pi_incr_factor, 'pi_incr_factor2':pi_incr_factor2,
                    'padding':padding,'padding_green1':padding_green1,'ifDifferential':ifDifferential,
                    'DAQ_error_factor':DAQ_error_factor}

        start = time.time()
        ConfocalODMRAWGFastObject = ConfocalODMRAWGFast(settings=settings, ifPlotPulse=0) 
        ConfocalODMRAWGFastObject.runScan()
        print('Total time = ' + str(time.time() - start) + ' s')
        ConfocalODMRAWGFastObject.close()
        

        # for x in xArray:
        #     # smallXArray = xArray[j*10:(j+1)*10] 
        #     # if j==19: smallXArray = xArray[190:201]   
        #     smallXArray = np.linspace(x,x,1)

        # smallXArray = xArray

        # start = 2610e6; stop =2710e6; num_sweep_points = 41
        # freqsArray = np.linspace(start, stop, num_sweep_points)
        # start = 3050e6; stop =3150e6; num_sweep_points = 41
        # freqsArray2 = np.linspace(start, stop, num_sweep_points)
        # freqsArray = np.concatenate((freqsArray,freqsArray2))

if __name__ == "__main__":
    main()