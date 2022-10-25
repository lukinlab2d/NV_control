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
from PIL import Image
from qcodes.actions import Task as qctask
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot
import numpy as np
import datetime as dt

import nidaqmx, datetime, ctypes, os, warnings, time, itertools
from nidaqmx.constants import *
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
import matplotlib as mpl
import matplotlib.pyplot as plt
from ODMR import *
import sqlite3
from qcodes.utils.dataset.doNd import plot

def ns2cycles(time_in_ns, samp_rate=1e7):
        return int(time_in_ns/1e9*samp_rate)
  

####################################################################################################################

start = 0; stop = 20; num_sweep_points = 21
freqsArray = np.linspace(start, stop, num_sweep_points)

num_loops               = int(3000);            #samp_rate_DAQ              = 1e7 # not important
laser_init_delay_in_ns  = 1e3;                  laser_init_duration_in_ns  = 1e3; when_init_end = laser_init_delay_in_ns+laser_init_duration_in_ns
AFG_delay_in_ns         = when_init_end+40;     AFG_duration_in_ns         = 40; when_pulse_end = AFG_delay_in_ns+AFG_duration_in_ns
laser_read_delay_in_ns  = when_pulse_end;       laser_read_duration_in_ns  = 2e3
read_signal_delay_in_ns = when_pulse_end+500;   read_signal_duration_in_ns = 300
read_ref_delay_in_ns    = when_pulse_end+1500;  read_ref_duration_in_ns    = 300


settings = {'start': start, 'stop': stop, 'num_sweep_points': num_sweep_points, 'num_loops':num_loops, #'samp_rate_DAQ': samp_rate_DAQ,
            'laser_init_delay_in_ns': laser_init_delay_in_ns,'laser_init_duration_in_ns': laser_init_duration_in_ns,
            'laser_read_delay_in_ns': laser_read_delay_in_ns,'laser_read_duration_in_ns': laser_read_duration_in_ns,
            'AFG_delay_in_ns':AFG_delay_in_ns, 'AFG_duration_in_ns':AFG_duration_in_ns,
            'read_signal_delay_in_ns':read_signal_delay_in_ns, 'read_signal_duration_in_ns':read_signal_duration_in_ns,
            'read_ref_delay_in_ns':read_ref_delay_in_ns, 'read_ref_duration_in_ns':read_ref_duration_in_ns}


ODMRObject = ODMR(settings=settings)
freqs = Freqs(name='freqs', instrument=ODMRObject)
doODMRAtOneFreq = ODMRObject.doODMRAtOneFreq # this is implemented as a Parameter

loop = Loop(
    freqs.sweep(start, stop, num=num_sweep_points),
    delay = 0,
    sleepTimeAfterFinishing=0).each(doODMRAtOneFreq).then(qctask(doODMRAtOneFreq.close))

data = loop.get_data_set()
data.add_metadata(settings)
print('data.location = ' + str(data.location))

plot = QtPlot(
    data.ODMRObject_doODMRAtOneFreq, # this is implemented as a Parameter
    figsize = (1200, 600),
    interval = 1,
    # theme = ((255, 255, 255), "black"), # color in quotes is background color
    )

loop.with_bg_task(plot.update)
loop.run()
plotfile = plot.save(type='data')
img = Image.open(plotfile)
img.show()

# dataFolder = os.path.join(os.getcwd(), 'data', ODMRObject.date)
# if not os.path.isdir(dataFolder): os.mkdir(dataFolder)
# dataFilePath = os.path.join(dataFolder, 'test.db')
# print(dataFilePath)
# initialise_or_create_database_at(dataFilePath)
# data.write(write_metadata=True, only_complete=False, filename=dataFilePath)

# print(data.path_to_db)
# data.write()
# print(dir(data))
# data.finalize()

##########################################################################

# Plot pulses
fig2, ax = plt.subplots()
PulseTrace = {}; pulse_sequence = ODMRObject.pulse_sequence
maxLength = np.max([pulse.start_time + pulse.duration for pulse in pulse_sequence])

for pulse in pulse_sequence:
    channel_id = pulse.channel_id
    trace_array = np.concatenate((np.zeros(int(pulse.start_time)),
                            np.ones(int(pulse.duration)),
                            np.zeros(int(maxLength - pulse.start_time - pulse.duration))))
    if not channel_id in PulseTrace: PulseTrace[channel_id] = Trace(trace_array, channel_id)
    else: PulseTrace[channel_id].arr = np.add(PulseTrace[channel_id].arr, trace_array)

for channel_id in PulseTrace:
    tr = PulseTrace[channel_id]
    x_axis = np.array(range(tr.length))/1e9*1e6
    if tr.name != "MWswitch": ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color)
    ax.set_xlabel("Time ($\mu$s)")
    ax.legend(loc='best')
plt.show()