"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""

import pyvisa as visa
import nidaqmx, time
import sys
from qcodes.instrument.base import Instrument
from qcodes.instrument.parameter import Parameter
from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity

vel = Velocity(velNum=2, ifInitVpz=0, ifInitWvl=1, initWvl=636.7)
vel.set_current(67)
vel.set_vpiezo(3)
