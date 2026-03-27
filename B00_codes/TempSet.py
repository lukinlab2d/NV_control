from qcodes_contrib_drivers.drivers.Lakeshore.Lakeshore335 import Lakeshore335
import pyvisa as visa
import sys
import time
from qcodes.instrument.base import Instrument
from qcodes.instrument.parameter import Parameter
from pyvisa.constants import Parity, StopBits, ControlFlow

def main(temp=2):
    ls = Lakeshore335()
    ls.set_setp(temp)
    ls.set_PID(p=90,i=90,d=0)
    ls.set_heaterRange(range=3)#2
    # ls.allOff()
    ls.close()

if __name__ == "__main__":
    main(temp=9.5)