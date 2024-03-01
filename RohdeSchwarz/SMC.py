import pyvisa as visa
import sys
from qcodes.instrument.base import Instrument
from qcodes.instrument.parameter import Parameter
import time

class SMC():
    def __init__(self, name='RohdeSchwarz_SMC100', SRSnum=3, settings=None, **kwargs) -> None:
        rm = visa.ResourceManager()
        if SRSnum == 3:
            self.smc = rm.open_resource('USB0::0x0AAD::0x006E::102580::INSTR')
            self.smc.query('*IDN?')

    def set_freq(self, freq, units='Hz'):
        self.smc.write('FREQ ' + str(freq) + units)

    def set_RFAmplitude(self, RFamplitude, units='dBm'):
        self.smc.write('POW:MODE CW')
        self.smc.write('POW '+ str(RFamplitude))
    
    def enable_RFOutput(self):
        self.smc.write('OUTP 1')

    def disable_RFOutput(self):
        self.smc.write('OUTP 0')
    
    def disableModulation(self):
        self.smc.write('MOD:STAT OFF')
        

# if __name__ == '__main__':
#     smc = SMC()
#     smc.set_RFAmplitude(-15)
#     smc.set_freq(2.8e9)
#     smc.enable_RFOutput()
#     time.sleep(10)
#     smc.disable_RFOutput()
        
    