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

laserTrack_channel     = 6;      if_tracking = 1
xy_scan_read_time      = 10;     xy_scan_settle_time    = 10;  
xy_scan_resolution_hor = 25;     xy_scan_resolution_ver = 25
x_minus_range          = 0.04;   x_plus_range           = 0.04
y_minus_range          = 0.04;   y_plus_range           = 0.04
xy_displacement_limit  = 0.04;   num_of_scans           = 2;    tracking_period = 1e9 # just dummy
time_btwn_trackings    = 1*60

xz_scan_resolution_hor = 25;     xz_scan_resolution_ver = 25
x_minus_range          = 0.04;   x_plus_range           = 0.04
z_minus_range          = 0.2;    z_plus_range           = 0.2
xz_displacement_limit  = 0.2; 

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