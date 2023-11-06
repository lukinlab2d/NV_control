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
import numpy as np
from nidaqmx.constants import *
from B00_codes.PlotPulse import *
from B00_codes.CalibReadoutPhotonStat import *
import B00_codes.dataReader as dataReader
from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity

####################################################################################################################
vel = Velocity()
for vpiezo in np.linspace(1,100,991):
    # CalibReadoutPhotonStat
    ifRandomized = 0; ifLooped = 0; ifPlotPulse = False
    ifMeaningfulRef=False; ifRefBright=False; ifRefInitAgain = False
    uwPower = -105; uwFreq = 2.870e9

    start = 0.5e6; stop = 5e6; num_sweep_points = 10 # sweep t_read
    tausArray = np.linspace(start, stop, num_sweep_points)

    current = 55;           vel.set_current(current)
    wvl = 637.2;            vel.set_wvl(wvl)
    vel.set_vpiezo(vpiezo); vel.set_ready()
    
    laserInit_channel=7; laserRead_channel=5; laserTrack_channel=7; laserIon_channel=10 # 532 is 3, S589 is 6, W589 is 9

    num_loops               = int(1e3)
    laser_init_delay        = 5e3;    laser_init_duration       = 50e3
    laser_to_MWI_delay      = 100e3;  pi_time                   = 60 # dummy
    laser_to_DAQ_delay_directory = {3: 850, 6: 1150, 9: 1150, 7: 900, 5: 1750}
    laser_to_DAQ_delay = laser_to_DAQ_delay_directory.get(laserRead_channel, 0)
    DAQ_to_laser_off_delay  = 5e2;    MWI_to_switch_delay       = 10 # cannot be between 0 and 10
    ifIon = 0
    laser_ion_duration      = 500;    ion_to_readLaser_delay    = 5e3
    
    if True:
        # For NV tracking
        if_tracking = 1
        xy_scan_read_time      = 5;      xy_scan_settle_time    = 3;  
        xy_scan_resolution_hor = 30;     xy_scan_resolution_ver = 30
        x_minus_range          = 0.05;   x_plus_range           = 0.05
        y_minus_range          = 0.05;   y_plus_range           = 0.05
        xy_displacement_limit  = 0.05;   num_of_scans           = 3;    tracking_period = 1e9 # just dummy
        time_btwn_trackings    = 15000*60

        xz_scan_resolution_hor = 30;     xz_scan_resolution_ver = 30
        x_minus_range          = 0.05;   x_plus_range           = 0.05
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

    settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'tausArray': tausArray,
                'num_loops':num_loops, 'uwPower':uwPower, 'uwFreq': uwFreq, 'wavelength': wvl, 'laser_I': current, 'vpiezo': vpiezo,
                'laser_init_delay':       laser_init_delay,        'laser_init_duration': laser_init_duration,
                'laser_to_MWI_delay':     laser_to_MWI_delay ,     'pi_time':             pi_time,
                'laser_to_DAQ_delay':     laser_to_DAQ_delay ,     
                'DAQ_to_laser_off_delay': DAQ_to_laser_off_delay,  'MWI_to_switch_delay': MWI_to_switch_delay,
                'ifRandomized':           ifRandomized,            'ifRefBright':         ifRefBright,
                'laserRead_channel':      laserRead_channel,       'laserInit_channel': laserInit_channel,
                'trackingSettings':       trackingSettings,        'ifMeaningfulRef':    ifMeaningfulRef,
                'ifRefInitAgain':         ifRefInitAgain,          'laserTrack_channel': laserTrack_channel,
                'laserIon_channel':       laserIon_channel,        'laser_ion_duration': laser_ion_duration,
                'ion_to_readLaser_delay': ion_to_readLaser_delay,  'ifIon': ifIon}

    start = time.time()
    print("wvl ", wvl)
    time.sleep(1)

    CalibReadoutPhotonStatObject = CalibReadoutPhotonStat(settings=settings, ifPlotPulse=ifPlotPulse) # this is implemented as an Instrument
    CalibReadoutPhotonStatObject.runScan()
    print('Total time = ' + str(time.time() - start) + ' s')
    
    dataFilename = CalibReadoutPhotonStatObject.getDataFilename()
    # if ifPlotPulse: dataReader.readDataFullData(dataFilename)
    CalibReadoutPhotonStatObject.close()

    # dataFilename = 'C:/Users/lukin2dmaterials/data/2023-05-23/#059_CalibReadoutPhotonStat_17-02-21/CalibReadoutPhotonStatObject_sig_set.dat'
    # dataReader.readDataFullData(dataFilename)
        




