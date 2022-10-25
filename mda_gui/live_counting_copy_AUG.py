# NI_9402 live counter read
# dates initial program written during 07-2022
# updates/changes started: 082222
# last change made on 082922

"""
this program using the QCoDeS experiment control interface to count the rising edges of an
input digital signal. The desired sampling rate can be set. Two channels are used on the NI-9402
module; both can be selected when setting up the instrument. One channel is for the counter
input and the second is for the internal hardware based clock. The counting process is
executed using the legacy QCoDeS Loop() fnc.
"""

################################## importing required modules ####################################

import qcodes as qc
import qcodes.actions
import nidaqmx
import time
from qcodes_contrib_drivers.drivers.NationalInstruments.DAQ import *
import IPython.lib.backgroundjobs as bg
from qcodes.loops import Loop
from qcodes.plots.pyqtgraph import QtPlot

################################## creating ni_4902 instrument ###################################

ni_9402_device_name = "cDAQ1Mod1"
ni_9402_counter_channel = "cDAQ1Mod1/ctr0"
ni_9402_clock_channel = "cDAQ1Mod1/ctr1"
ni_9402_source_channel = "/cDAQ1/Ctr1InternalOutput"
ni_9402_sampling_rate = 1e5 #10
ni_9402_duty_cycle = 0.5                                  # is this for the clock?
ni_9402_integration_time = 0.03 # 0.01
ni_9402_samples_per_channel = int(ni_9402_integration_time * ni_9402_sampling_rate) # 1
ni_9402_timeout = 60

ni_9402 = Counter("ni_9402_module", ni_9402_device_name, ni_9402_counter_channel,
    ni_9402_clock_channel, ni_9402_source_channel, ni_9402_sampling_rate,
    ni_9402_samples_per_channel, ni_9402_duty_cycle,
    ni_9402_integration_time, ni_9402_timeout)

############################## main fnc (counting and plotting) ###################################

p_measure = qc.ManualParameter(name = "counts_num")
p_sweep = qc.Parameter(name = "time_s", set_cmd = p_measure.set)

fnc_2 = ni_9402.integrateavg

loop = Loop(
    # p_sweep.sweep(0, 10000, step = 0.01), # .sweep(start, stop, num values to generate)
    p_sweep.sweep(0, 10000, step = 1), # .sweep(start, stop, num values to generate)
    delay = 0.0
    ).each(fnc_2)

data = loop.get_data_set()
# print(dir(data))

plot = QtPlot(
    data.ni_9402_module_integrateavg, # ni_9402_module_integrateavg
    figsize = (1200, 600),
    interval = 1,
    theme = ((255, 255, 255), "black"), # color in quotes is background color
    )

loop.with_bg_task(plot.update)
loop.run()


# plot = QtPlot(
#     figsize = (1200, 600),
#     interval = 1,
#     theme = ((255, 255, 255), "black"), # color in quotes is background color
#     )

# loop = Loop(
#     # p_sweep.sweep(0, 10000, step = 0.01), # .sweep(start, stop, num values to generate)
#     p_sweep.sweep(0, 10000, step = 1), # .sweep(start, stop, num values to generate)
#     delay = 0.0
#     ).each(
#         fnc_2,
#         print(hasattr(qc.loops.active_data_set, 'data_set')),
#         plot.add(qc.loops.active_data_set()),
#         # plot.add_updater(plot_config=qc.config),
#         plot.update()
#         )

# loop.run()