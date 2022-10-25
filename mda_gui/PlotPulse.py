import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image

class PulseTrace:
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


class PlotPulse():
    def __init__(self, plotFilename=None, measurementObject=None, pulseSequence=None, ifShown=True, ifSave=True):
        self.pulseTrace = {}; 
        if measurementObject is not None:
            self.pulse_sequence = measurementObject.pulse_sequence
        elif pulseSequence is not None:
            self.pulse_sequence = pulseSequence
        else:
            raise Exception("No pulse sequences to plot!!")

        self.maxLength = np.max([pulse.start_time + pulse.duration for pulse in self.pulse_sequence])
        self.ifShown = ifShown; self.ifSave = ifSave
        if ifSave:
            if plotFilename is not None: self.plotFilename = plotFilename
            else:
                raise Exception("Filename for pulse plot missing")

    def makeTraceDict(self):
        for pulse in self.pulse_sequence:
            channel_id = pulse.channel_id
            trace_array = np.concatenate((np.zeros(int(pulse.start_time)),
                                    np.ones(int(pulse.duration)),
                                    np.zeros(int(self.maxLength - pulse.start_time - pulse.duration))))
            if not channel_id in self.pulseTrace: self.pulseTrace[channel_id] = PulseTrace(trace_array, channel_id)
            else: self.pulseTrace[channel_id].arr = np.add(self.pulseTrace[channel_id].arr, trace_array)
    
    def makePulsePlot(self):
        self.makeTraceDict()

        fig = plt.figure(1)
        ax = fig.gca()
        ax.clear()
        for channel_id in self.pulseTrace:
            tr = self.pulseTrace[channel_id]
            x_axis = np.array(range(tr.length))/1e9*1e6
            if tr.name != "MWswitch": ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color)
            ax.set_xlabel("Time ($\mu$s)")
            ax.legend(loc='best')

        if self.ifSave: 
            plt.savefig(self.plotFilename)
            print('Pulse plot saved to ' + self.plotFilename)

            if self.ifShown: 
                plt.show(block=False)
                plt.pause(0.05)
                self.showPulsePlot()

        else:
            if self.ifShown:
                plt.show(block=False)
                plt.pause(0.05)


        return fig
    
    def showPulsePlot(self):
        img = Image.open(self.plotFilename)
        img.show()