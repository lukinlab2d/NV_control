import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image

class PulseTrace:
    def __init__(self, arr, name, readColor=3, initColor=3, ionColor=8, shelveColor=3):
        self.name = name
        self.arr = arr

        if "aserInit" in name or name == "Laser": 
            self.vert_offset = 6
            if initColor == 6 or initColor == 9: self.color = 'C1'
            else: self.color = 'C2'
        elif "aserRead2" in name:
            self.vert_offset = 12
            if readColor == 6 or readColor == 9: self.color = 'C1'
            elif readColor == 5 or readColor == 14: self.color = 'r'
            else: self.color = 'C2'
        elif "aserRead" in name:
            self.vert_offset = 6
            if readColor == 6 or readColor == 9: self.color = 'C1'
            elif readColor == 5 or readColor == 14: self.color = 'r'
            else: self.color = 'C2'
        elif "aserIon2" in name:
            self.vert_offset = 12
            if ionColor == 6 or ionColor == 9: self.color = 'C1'
            elif ionColor == 3: self.color = 'C2'
            else: self.color = 'darkred'
        elif "aserIon" in name:
            self.vert_offset = 6
            if ionColor == 6 or ionColor == 9: self.color = 'C1'
            elif ionColor == 3: self.color = 'C2'
            else: self.color = 'darkred'
        elif "aserShelve2" in name:
            self.vert_offset = 12
            if shelveColor == 6 or shelveColor == 9: self.color = 'C1'
            elif shelveColor == 3: self.color = 'C2'
            else: self.color = 'r'
        elif "aserShelve" in name:
            self.vert_offset = 6
            if shelveColor == 6 or shelveColor == 9: self.color = 'C1'
            elif shelveColor == 3: self.color = 'C2'
            else: self.color = 'r'
        elif "MW_I2" in name or "AFG2" in name: 
            self.vert_offset = 10.5; self.color = 'C0'
        elif "MW_I" in name or "AFG" in name: 
            self.vert_offset = 4.5; self.color = 'C0'
        elif "MW_Q2" in name: 
            self.vert_offset = 9; self.color = 'C1'
        elif "MW_Q" in name: 
            self.vert_offset = 3; self.color = 'C1'
        elif "MWswitch2" in name: 
            self.vert_offset = 7.5; self.color = 'C0'
        elif "MWswitch" in name: 
            self.vert_offset = 1.5; self.color = 'C0'
        elif "ounter" in name: 
            self.vert_offset = 0; self.color = 'k'

        self.arr = [x + self.vert_offset for x in arr]
        self.length = len(arr)


class PlotPulse():
    def __init__(self, plotFilename=None, measurementObject=None, pulseSequence=None, 
                 ifShown=True, ifSave=True, readColor=3, initColor=3, ionColor=8, shelveColor=3):
        self.pulseTrace = {}; 
        if measurementObject is not None:
            self.pulse_sequence = measurementObject.pulse_sequence
        elif pulseSequence is not None:
            self.pulse_sequence = pulseSequence
        else:
            raise Exception("No pulse sequences to plot!!")

        self.maxLength = np.max([pulse.start_time + pulse.duration for pulse in self.pulse_sequence])
        self.ifShown = ifShown; self.ifSave = ifSave; self.readColor = readColor; 
        self.initColor = initColor; self.ionColor = ionColor; self.shelveColor = shelveColor
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
            if not channel_id in self.pulseTrace: 
                self.pulseTrace[channel_id] = PulseTrace(trace_array, channel_id, readColor=self.readColor, 
                                                         initColor=self.initColor, ionColor=self.ionColor, shelveColor=self.shelveColor)
            else: 
                self.pulseTrace[channel_id].arr = np.add(self.pulseTrace[channel_id].arr, trace_array)
    
    def makePulsePlot(self):
        self.makeTraceDict()

        fig = plt.figure(1)
        ax = fig.gca()
        ax.clear()
        for channel_id in self.pulseTrace:
            tr = self.pulseTrace[channel_id]
            x_axis = np.array(range(tr.length))/1e9*1e6
            if '2' in tr.name:
                ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color, linestyle='--')
            else:
                ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color)
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