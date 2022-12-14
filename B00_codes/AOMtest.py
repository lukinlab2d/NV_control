"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
"""

from operator import contains
from typing import Dict, Union
from collections import namedtuple
import numpy as np
import datetime as dt
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS
from B00_codes import TurnOnLaser, TurnOffLaser

import nidaqmx, datetime, ctypes, os, warnings, time, itertools
from nidaqmx.constants import *
from qcodes.instrument.base import Instrument
# from qcodes.instrument.parameter import ParameterWithSetpoints, Parameter
# from qcodes.utils.validators import *
from pylabcontrol.core.read_write_functions import get_config_value
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
import matplotlib as mpl
import matplotlib.pyplot as plt

def ns2cycles(time_in_ns, samp_rate=1e7):
        return int(time_in_ns/1e9*samp_rate)
    
class Trace:
    def __init__(self, arr, name):
        self.name = name
        self.arr = arr

        if "aser" in name: 
            self.vert_offset = 3; self.color = 'C2'
        elif "AFG" in name or "MWswitch" in name: 
            self.vert_offset = 1.5; self.color = 'C0'
        elif "ounter" in name: 
            self.vert_offset = 0; self.color = 'k'

        self.arr = [x + self.vert_offset for x in arr]
        self.length = len(arr)

####################################################################################################################
if __name__ == '__main__':
    print()
    now = dt.datetime.now(); date = now.strftime("%Y-%m-%d")
    print(now)
    
    for i in range(1):
       
        # clock speed is in MHz - is 'status' needed in the dictionary?
        clock_speed = 500 # MHz
        LaserParam = {'delay_time': 2, 'channel':1}
        CounterParam = {'delay_time': 2, 'channel':2}
        # uwFreq = 80e6; uwPower = -7
        # srs = SRS()
        # srs.set_freq(uwFreq) #Hz
        # srs.set_RFAmplitude(uwPower) #dBm
        # srs.enable_RFOutput()


        # num_loops = int(1e6)
        # laser_delay_in_ns  = 0.2e3;  laser_duration_in_ns   = 20
        # read_delay_in_ns   = 0.2e3;  read_duration_in_ns    = 20


        settings = {'clock_speed': clock_speed, 'Laser': LaserParam, 'Counter': CounterParam, 
                    'PB_type': 'USB', 'min_pulse_dur': int(1*1e3/clock_speed)}

        # pb = spc.B00PulseBlaster("SpinCorePB", settings=settings, verbose=False)
        # pb.turn_on_infinite(channel=2)
        # time.sleep(3)
        # pb.turn_off()

        channels = np.linspace(3,0,1)
        pb = TurnOnLaser.turnOnLaser(channels=channels)