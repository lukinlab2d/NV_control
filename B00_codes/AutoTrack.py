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

laserTrack_channel     = 7;       if_tracking = 1
xy_scan_read_time      = 50;      xy_scan_settle_time    = 30;  
xy_scan_resolution_hor = 20;      xy_scan_resolution_ver = 20
x_minus_range          = 0.01;    x_plus_range           = 0.01
y_minus_range          = 0.01;    y_plus_range           = 0.01
xy_displacement_limit  = 0.015;    num_of_scans           = 2
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
                    'tracking_period':        1e9,                   'if_tracking':            if_tracking,
                    'xz_scan_resolution_hor': xz_scan_resolution_hor,'xz_scan_resolution_ver': xz_scan_resolution_ver,
                    'x_minus_range':          x_minus_range ,        'x_plus_range':           x_plus_range,
                    'z_minus_range':          z_minus_range ,        'z_plus_range':           z_plus_range,
                    'xz_displacement_limit':  xz_displacement_limit, 'time_btwn_trackings':    time_btwn_trackings}


cfcObject = Confocal(settings=trackingSettings, laserChannel=laserTrack_channel)
# cfcObject.optimize_xy_fast()
# time.sleep(1)
# cfcObject.optimize_xz()
# time.sleep(1)
x1, y1, z = cfcObject.optimize_xy(direction=1)
x2, y2, z = cfcObject.optimize_xy(direction=-1)
cfcObject.set_coordinate_fnc((x1+x2)/2, (y1+y2)/2, z)
time.sleep(1)
cfcObject.close()