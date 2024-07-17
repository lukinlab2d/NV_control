"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image

class PulseTrace:
    def __init__(self, arr, name, readColor=3, initColor=3, ionColor=8, shelveColor=3):
        self.name = name
        self.arr = arr        

        if "MWswitch2" in name: 
            self.vert_offset = 16; self.color = 'C3'
        elif "MWswitch4" in name: 
            self.vert_offset = 15.5; self.color = 'C4'

        elif "MW_I2" in name or "AFG2" in name: 
            self.vert_offset = 14; self.color = 'C0'
        elif "MW_Q2" in name: 
            self.vert_offset = 13.5; self.color = 'C1'
        
        elif "aserRead2" in name:
            self.vert_offset = 12
            if readColor == 6 or readColor == 9: self.color = 'C1'
            elif readColor == 5 or readColor == 14: self.color = 'r'
            else: self.color = 'C2'
        elif "aserShelve2" in name:
            self.vert_offset = 12
            if shelveColor == 6 or shelveColor == 9: self.color = 'C1'
            elif shelveColor == 3 or shelveColor == 7: self.color = 'C2'
            else: self.color = 'r'
        elif "aserIon2" in name:
            self.vert_offset = 11.5
            if ionColor == 6 or ionColor == 9: self.color = 'C1'
            elif ionColor == 3: self.color = 'C2'
            else: self.color = 'darkred'
        
        ##########################################################################

        elif "aserInit" in name or name == "Laser": 
            self.vert_offset = 10
            if initColor == 6 or initColor == 9: self.color = 'C1'
            else: self.color = 'C2'
        elif "hiLoMWPwr" in name: 
            self.vert_offset = 8; self.color = 'k'

        ##########################################################################

        elif "aserRead" in name:
            self.vert_offset = 6
            if readColor == 6 or readColor == 9: self.color = 'C1'
            elif readColor == 5 or readColor == 14: self.color = 'r'
            else: self.color = 'C2'
        elif "aserShelve" in name:
            self.vert_offset = 6
            if shelveColor == 6 or shelveColor == 9: self.color = 'C1'
            elif shelveColor == 3 or shelveColor == 7: self.color = 'C2'
            else: self.color = 'r'
        elif "aserIon" in name:
            self.vert_offset = 5.5
            if ionColor == 6 or ionColor == 9: self.color = 'C1'
            elif ionColor == 3: self.color = 'C2'
            else: self.color = 'darkred'

        elif "MW_I" in name or "AFG" in name or "AWG" in name: 
            self.vert_offset = 4; self.color = 'C0'
        elif "MW_Q" in name: 
            self.vert_offset = 3.5; self.color = 'C1'

        elif name == "MWswitch": 
            self.vert_offset = 2; self.color = 'C3'
        elif "MWswitch3" in name: 
            self.vert_offset = 1.5; self.color = 'C4'

        elif "Counter" in name: 
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
            if '2' in tr.name or '4' in tr.name or 'hiLo' in tr.name:
                ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color, linestyle='--')
            else:
                ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color)
            ax.set_xlabel("Time ($\mu$s)")
            ax.legend(loc='right', bbox_to_anchor=(1.35, 0.5))
            plt.tight_layout()

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
    
    def makePulsePlotAWG(self, ch1plot, ch2plot, MW_del):
        self.makeTraceDict()

        fig = plt.figure(1)
        ax = fig.gca()
        ax.clear()

        for channel_id in self.pulseTrace:
            tr = self.pulseTrace[channel_id]
            x_axis = np.array(range(tr.length))/1e9*1e6
            if '2' in tr.name or '4' in tr.name or 'hiLo' in tr.name:
                ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color, linestyle='--')
            else:
                ax.plot(x_axis, tr.arr, label=tr.name, color=tr.color)
            
        if len(x_axis)>=(len(ch1plot)+int(MW_del)):
            diff = len(x_axis) - len(ch1plot) - int(MW_del)
            ch1 = np.concatenate((np.zeros(int(MW_del)), ch1plot, np.zeros(int(diff)))) + 16.75
            ch2 = np.concatenate((np.zeros(int(MW_del)), ch2plot, np.zeros(int(diff)))) + 16.5
            ax.plot(x_axis, ch1, label='MW1')
            ax.plot(x_axis, ch2, label='MW2')
        
        
        ax.set_xlabel("Time ($\mu$s)")
        ax.legend(loc='right', bbox_to_anchor=(1.35, 0.5))
        plt.tight_layout()

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