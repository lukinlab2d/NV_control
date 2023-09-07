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
from copy import deepcopy
from qcodes.actions import Task as qctask
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot
import numpy as np
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc
from qcodes_contrib_drivers.drivers.StanfordResearchSystems.SG386 import SRS

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

laserTrack_channel     = 3;       if_tracking = 1
xy_scan_read_time      = 7;       xy_scan_settle_time    = 4;  
xy_scan_resolution_hor = 25;      xy_scan_resolution_ver = 25
x_minus_range          = 0.02;   x_plus_range           = 0.02
y_minus_range          = 0.02;   y_plus_range           = 0.02
xy_displacement_limit  = 0.02;   num_of_scans           = 2;    tracking_period = 1e9 # just dummy
time_btwn_trackings    = 1*60

xz_scan_resolution_hor = 20;      xz_scan_resolution_ver = 25
x_minus_range          = 0.02;   x_plus_range           = 0.02
z_minus_range          = 0.1;    z_plus_range           = 0.1
xz_displacement_limit  = 0.1; 

trackingSettings = {'xy_scan_read_time':      xy_scan_read_time,     'xy_scan_settle_time':    xy_scan_settle_time,
                    'xy_scan_resolution_hor': xy_scan_resolution_hor,'xy_scan_resolution_ver': xy_scan_resolution_ver,
                    'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                    'y_minus_range':          y_minus_range ,        'y_plus_range':           y_plus_range,
                    'xy_displacement_limit':  xy_displacement_limit, 'num_of_scans':           num_of_scans,
                    'tracking_period':        tracking_period,       'if_tracking':            if_tracking,
                    'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                    'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                    'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                    'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}


cfcObject = Confocal(settings=trackingSettings, laserChannel=laserTrack_channel)
cfcObject.optimize_xy_fast()
time.sleep(1)
cfcObject.optimize_xz()
time.sleep(1)
cfcObject.optimize_xy()
time.sleep(1)
cfcObject.close()