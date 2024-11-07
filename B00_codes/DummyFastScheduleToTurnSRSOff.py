"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.Rabi import *
from qcodes_contrib_drivers.drivers.Agilent.Agilent_33522A import AG33522A
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS
from qcodes_contrib_drivers.drivers.Siglent.SDG6022X import SDG6022X

####################################################################################################################

for SRSnum in np.array((1,2,3,4)):
    srs = SRS(SRSnum=SRSnum)
    srs.set_freq(2870e6) #Hz
    srs.set_RFAmplitude(-80) #dBm
    srs.disable_RFOutput()
    srs.disableModulation()
    print('SRS %.0f turned off' % SRSnum)
        

# Turn off AWGs
for SDGnum in np.array((1,2)):
    AWG = SDG6022X(name='SDG6022X', SDGnum=SDGnum)
    AWG.turn_off()
    print('SDG %.0f turned off' % SDGnum)

# # Turn off test source
# AG = AG33522A()
# AG.disable_PM()
# AG.disable_RFOutput()
# print('Agilent turned off')

