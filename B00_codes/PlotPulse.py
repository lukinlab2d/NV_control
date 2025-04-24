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

        elif "MW_I2" in name or "AFG2" in name or "AWG2" in name: 
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

        elif name == "MWswitch" or name == "noise": 
            self.vert_offset = 2; self.color = 'C3'
        elif "MWswitch3" in name: 
            self.vert_offset = 1.5; self.color = 'C4'

        elif "Counter" in name: 
            self.vert_offset = 0; self.color = 'k'

        self.arr = arr + self.vert_offset
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
        self.makeTraceDict()

    def makeTraceDict(self):
        for pulse in self.pulse_sequence:
            channel_id = pulse.channel_id
            # print(channel_id)
            trace_array = np.concatenate((np.zeros(int(pulse.start_time)),
                                    np.ones(int(pulse.duration)),
                                    np.zeros(int(self.maxLength - pulse.start_time - pulse.duration))))
            if not channel_id in self.pulseTrace: 
                self.pulseTrace[channel_id] = PulseTrace(trace_array, channel_id, readColor=self.readColor, 
                                                         initColor=self.initColor, ionColor=self.ionColor, shelveColor=self.shelveColor)
            else: 
                self.pulseTrace[channel_id].arr = np.add(self.pulseTrace[channel_id].arr, trace_array)
    
    def makePulsePlot(self):
        # self.makeTraceDict()

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
    
    def makePulsePlotAWG(self, ch1plot, ch2plot, MW_del, fig=None,
                         offset1=16.75, offset2=16.5, label1='MW1', label2='MW2'):
        # self.makeTraceDict()

        if fig is None:
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
                ch1 = np.concatenate((np.zeros(int(MW_del)), ch1plot, np.zeros(int(diff)))) + offset1
                ch2 = np.concatenate((np.zeros(int(MW_del)), ch2plot, np.zeros(int(diff)))) + offset2
            else:
                diff = (len(ch1plot) + int(MW_del)) - len(x_axis)
                ch1 = np.concatenate((np.zeros(int(MW_del)), ch1plot)) + offset1
                ch2 = np.concatenate((np.zeros(int(MW_del)), ch2plot)) + offset2
                step_x = x_axis[1] - x_axis[0]
                next_x = x_axis[-1] + step_x
                extra_xs = np.linspace(next_x, next_x+(diff-1)*step_x ,diff)
                x_axis = np.concatenate((x_axis, extra_xs))
            ax.plot(x_axis, ch1, label=label1)
            ax.plot(x_axis, ch2, label=label2)
            
            ax.set_xlabel("Time ($\mu$s)")
            ax.legend(loc='right', bbox_to_anchor=(1.35, 0.5))
            plt.tight_layout()
        
        else: # to plot more AWG traces
            ax = fig.gca()
            tr = next(iter(self.pulseTrace.values()))
            x_axis = np.array(range(tr.length))/1e9*1e6

            diff = len(x_axis) - (len(ch1plot) + int(MW_del))
            if diff >= 0:
                ch1 = np.concatenate((np.zeros(int(MW_del)), ch1plot, np.zeros(int(diff)))) + offset1
                ch2 = np.concatenate((np.zeros(int(MW_del)), ch2plot, np.zeros(int(diff)))) + offset2
            else:
                ch1 = np.concatenate((np.zeros(int(MW_del)), ch1plot)) + offset1
                ch2 = np.concatenate((np.zeros(int(MW_del)), ch2plot)) + offset2
                x_axis = np.array(range(len(ch1)))/1e9*1e6
                
            ax.plot(x_axis, ch1, label=label1)
            ax.plot(x_axis, ch2, label=label2)

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
    
    def makeTraceAWG(self, ch1plot, ch2plot, MW_del):
        
        tr = next(iter(self.pulseTrace.values()))
        x_axis = np.array(range(tr.length))/1e9*1e6
        
        if len(x_axis)>=(len(ch1plot)+int(MW_del)):
            diff = len(x_axis) - len(ch1plot) - int(MW_del)
            ch1 = np.concatenate((np.zeros(int(MW_del)), ch1plot, np.zeros(int(diff))))
            ch2 = np.concatenate((np.zeros(int(MW_del)), ch2plot, np.zeros(int(diff))))
        else:
            diff = (len(ch1plot) + int(MW_del)) - len(x_axis)
            ch1 = np.concatenate((np.zeros(int(MW_del)), ch1plot))
            ch2 = np.concatenate((np.zeros(int(MW_del)), ch2plot))
            step_x = x_axis[1] - x_axis[0]
            next_x = x_axis[-1] + step_x
            extra_xs = np.linspace(next_x, next_x+(diff-1)*step_x ,diff)
            x_axis = np.concatenate((x_axis, extra_xs))

        return x_axis, ch1, ch2
    
    def findTraceMW34(self):
        MW34 = None
        for channel_id in self.pulseTrace:
            if channel_id == 'MWswitch3' or channel_id == 'MWswitch4':
                tr = self.pulseTrace[channel_id]
                x_axis = np.array(range(tr.length))/1e9*1e6
                if MW34 is None:
                    MW34 = tr.arr - np.min(tr.arr)
                else:
                    MW34 = MW34 + tr.arr - np.min(tr.arr)
        
        return x_axis, MW34

    def find_zero_segments(self,arr):
        # Step 1: Create a boolean array where True represents zero and False represents non-zero
        is_zero = arr == 0
        
        # Step 2: Find the indices where the value changes (from non-zero to zero and vice versa)
        diff = np.diff(is_zero.astype(int))  # Convert boolean to int, then calculate the difference
        
        # Step 3: Find start and end points of zero segments
        starts = np.where(diff == 1)[0] + 1  # Where zero starts (after a non-zero)
        ends = np.where(diff == -1)[0]       # Where zero ends (before a non-zero)
        
        # If the array starts with zeros, include the first segment
        if is_zero[0]:
            starts = np.insert(starts, 0, 0)
        
        # If the array ends with zeros, include the last segment
        if is_zero[-1]:
            ends = np.append(ends, len(arr) - 1)
        
        # Step 4: Calculate lengths of zero segments
        lengths = ends - starts + 1
        
        return list(zip(starts, lengths))

        # # Example usage:
        # arr = np.array([1, 0, 0, 2, 0, 3, 0, 0, 0, 4, 0])
        # segments = find_zero_segments(arr)

        # print("Zero segments (start index, length):")
        # for start, length in segments:
        #     print(f"Start: {start}, Length: {length}")

    def showPulsePlot(self):
        img = Image.open(self.plotFilename)
        img.show()