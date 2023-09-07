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
from nidaqmx.constants import(
    Edge,
    CountDirection,
    AcquisitionType,
    FrequencyUnits
)
import sys
import os
import nidaqmx
import pyqtgraph as pg
import time
from datetime import date
from B00_codes import TurnOnLaser, TurnOffLaser
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc

from qcodes_contrib_drivers.drivers.NationalInstruments.DAQ import *
from qcodes_contrib_drivers.drivers.NationalInstruments.class_file import *

####################################################################################################################
class Confocal():
    def __init__(self, name='CFCScanObject', settings=None, laserChannel=3, **kwargs) -> None:
        self.oldX, self.oldY, self.oldZ = self.read_location()
        self.samp_rate = 1e5
        self.isSnake = 0
        self.settings = settings
        global pb
        channels = np.linspace(laserChannel,laserChannel,1)
        pb = TurnOnLaser.turnOnLaser(channels=channels, instrument_name="PB_NVtracking")
    
    def read_location(self):
        with open('C:/Users/lukin2dmaterials/miniconda3/envs/NV_control/B00_codes/NVlocation.txt', 'r') as f:
            s = f.readlines()
            oldX = float(s[0][0:-1])
            oldY = float(s[1][0:-1])
            oldZ = float(s[2])
            f.close()
        return oldX, oldY, oldZ

    def run_xy_scan_fnc(self, ifSpecifyStartEnd=0):
        self.isStopped = 0
        
        scan_galvo_card_name = "cDAQ1Mod2" # naming the instrument
        scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)} # dictionary of analog output channels
        scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels) # defining the instrument (ni_9263)

        self.scanType = 'xy'

        ############################################################################### def other variables #####################################################################
        ################### setting variales and array ####################

        # acquisition and settling time: 5 0.2 25 20
        xy_scan_read_time      = self.settings['xy_scan_read_time'];      xy_scan_settle_time    = self.settings['xy_scan_settle_time']
        xy_scan_resolution_hor = self.settings['xy_scan_resolution_hor']; xy_scan_resolution_ver = self.settings['xy_scan_resolution_ver']
        
        self.time_acquire = np.round(float(xy_scan_read_time)/1e3, 4)
        self.time_settle = np.round(float(xy_scan_settle_time)/1e3, 5)
        self.time_pts_settle = int(self.time_settle*self.samp_rate)
        self.time_pts_acquire = int(self.time_acquire*self.samp_rate)
        self.time_pts_per_pixel = self.time_pts_acquire  + self.time_pts_settle
        
        # resolution
        self.hor_res = grid_size_x = int(xy_scan_resolution_hor)
        self.ver_res = grid_size_y = int(xy_scan_resolution_ver)

        # initial driving voltages for the x,y-mirror and z piezo stage
        self.hor_init = x_init = self.oldX - self.settings['x_minus_range']
        self.ver_init = y_init = self.oldY - self.settings['y_minus_range']
        self.z = z_init = self.oldZ

        # final driving voltages for the x,y-mirror
        x_final = self.oldX + self.settings['x_plus_range']
        y_final = self.oldY + self.settings['y_plus_range']

        # range of voltage
        self.hor_range = x_final - x_init
        self.ver_range = y_final - y_init

        # array of x,y voltages
        self.x_array = np.linspace(x_init, x_final, grid_size_x); self.x_array_original = self.x_array
        self.y_array = np.linspace(y_init, y_final, grid_size_y)

        if ifSpecifyStartEnd:
            x_init = float(self.settings['xmin']); x_final = float(self.settings['xmax'])
            y_init = float(self.settings['ymin']); y_final = float(self.settings['ymax'])
            self.x_array = np.linspace(x_init, x_final, grid_size_x); self.x_array_original = self.x_array
            self.y_array = np.linspace(y_init, y_final, grid_size_y) 

        X, Y = np.meshgrid(self.x_array, self.y_array)

        # make the x-voltage waveform to pass to aotask
        self.x_array = np.array(np.repeat(self.x_array, self.time_pts_per_pixel))
        # print("len(self.x_array): " + str(len(self.x_array)))

        # create dataset to populate
        global xy_scan_data_array
        xy_scan_data_array = np.zeros((grid_size_y, grid_size_x))
        global most_recent_data_array
        most_recent_data_array = np.zeros((grid_size_y, grid_size_x))

        ################### resetting position of mirrors ####################

        scan_galvo.voltage_cdaq1mod2ao0(x_init) # x-mirror
        scan_galvo.voltage_cdaq1mod2ao1(y_init) # y-mirror
        scan_galvo.voltage_cdaq1mod2ao2(z_init) # z-stage
        
        # the clock to sync AO and counter
        clktask = nidaqmx.Task() 
        clktask.co_channels.add_co_pulse_chan_freq(  # adding dig pulse train chan
            counter = "cDAQ1Mod1/ctr1",
            name_to_assign_to_channel = "",
            units = nidaqmx.constants.FrequencyUnits.HZ,
            idle_state = nidaqmx.constants.Level.LOW,
            initial_delay = 0.0,
            freq = self.samp_rate,
            duty_cycle = 0.5
            )
        clktask.timing.cfg_implicit_timing( # implicit timing by the hardware
            sample_mode = AcquisitionType.CONTINUOUS, # the clock should run continuously in principle
            samps_per_chan = int(len(self.x_array) + 1) # does this matter?
            )
                
        ######################################################################## X and Y scanning #########################################################################
        print('----------------------------------------------------------------')
        print('XY scan started')
        for f in range(grid_size_y): # this loops for each row (y)
            if self.isStopped == 1: break
            start0 = time.time()

            # set initial y location
            scan_galvo.voltage_cdaq1mod2ao1(self.y_array[f])
            time.sleep(self.time_settle)
            
            # snake pattern: scanning right to left
            if self.isSnake == 1 and f % 2 == 1: self.x_array_write = np.ascontiguousarray(self.x_array[::-1])
            else: self.x_array_write = self.x_array

            ##############################################################
            # counter task
            ctrtask = nidaqmx.Task()
            ctrtask.ci_channels.add_ci_count_edges_chan( # define the counter
                counter = "cDAQ1Mod1/ctr0",
                name_to_assign_to_channel = "",
                edge = nidaqmx.constants.Edge.RISING,
                initial_count = 0,
                count_direction = nidaqmx.constants.CountDirection.COUNT_UP
                )
            ctrtask.timing.cfg_samp_clk_timing( # cfg sample clk timing
                rate = self.samp_rate,
                source = "/cDAQ1/Ctr1InternalOutput", # the clock defined above
                active_edge = nidaqmx.constants.Edge.RISING,
                sample_mode = AcquisitionType.FINITE,
                samps_per_chan = int(len(self.x_array) + 1) 
                )

            ##############################################################
            # AO task for the x-galvo
            aotask = nidaqmx.Task() 
            channel_galvo_x = f'{scan_galvo_card_name}/ao{0}'
            aotask.ao_channels.add_ao_voltage_chan(channel_galvo_x)
            aotask.timing.cfg_samp_clk_timing(
                rate = self.samp_rate,
                source = "/cDAQ1/Ctr1InternalOutput", # the clock defined above
                active_edge = nidaqmx.constants.Edge.RISING,
                sample_mode = AcquisitionType.FINITE,
                samps_per_chan = len(self.x_array)
                )
            aotask.write(self.x_array_write, auto_start=False) # assign the x-waveform to aotask, but not start yet

            # start galvo X scan, counter, and clock
            aotask.start()
            ctrtask.start()
            clktask.start() # the clock MUST START AFTER aotask.start() and ctrtask.start()
            
            aotask.wait_until_done()
            xLineData = ctrtask.read(len(self.x_array) + 1)  # +1 is to take the difference later
            # print("len(xLineData): " + str(len(xLineData)))
            
            aotask.stop(); aotask.close()
            ctrtask.stop(); ctrtask.close()
            clktask.stop() 
            diffData = np.diff(xLineData)

            summedData = np.zeros(grid_size_x)
            for i in range(0, grid_size_x):
                summedData[i] = np.sum(
                    diffData[(i*self.time_pts_per_pixel + self.time_pts_settle):(i*self.time_pts_per_pixel + self.time_pts_per_pixel - 1)])
            normedData = summedData * (.001 / self.time_acquire) # normalize to kcounts/sec

            if self.isSnake == 1 and f % 2 == 1: xy_scan_data_array[f] = np.flip(normedData)
            else: xy_scan_data_array[f] = normedData

            end = time.time() - start0
            # if f == 0 or f == grid_size_y-1 or f == int(grid_size_y/2):
            #     # print("Time per row: " + str(end) + " s")
                
        clktask.close(); scan_galvo.close() 

        ############################################################### end scanning script #############################################################################################
        most_recent_data_array = xy_scan_data_array # make temp holding global data_array the same as xy_scan_data_array
        
        indices = np.where(most_recent_data_array == most_recent_data_array.max())
        xMax_idx = indices[1][0]; yMax_idx = indices[0][0]
        xMax = self.x_array_original[xMax_idx]; yMax = self.y_array[yMax_idx]
        # print('xMax: ' + str(xMax))
        # print('yMax: ' + str(yMax))
        print('XY scan finished')
        # print('------------------------------------------------------------------')
        # plt.imshow(most_recent_data_array)
        # plt.show()
        return xMax, yMax, most_recent_data_array

    def optimize_xy(self):
        xMaxArr = []; yMaxArr = []; cMaxArr = []
        for i in range(self.settings['num_of_scans']):
            xMax, yMax, most_recent_data_array = self.run_xy_scan_fnc()
            xMaxArr.append(xMax); yMaxArr.append(yMax); cMaxArr.append(most_recent_data_array.max())
        xMaxArr = np.array(xMaxArr); yMaxArr = np.array(yMaxArr); cMaxArr = np.array(cMaxArr)
        xMaxAvg = np.round(np.average(xMaxArr),3); yMaxAvg = np.round(np.average(yMaxArr),3); cMaxAvg = np.round(np.average(cMaxArr),3)
        print('xMaxAvg: ' + str(xMaxAvg))
        print('yMaxAvg: ' + str(yMaxAvg))
        print('cMaxAvg: ' + str(cMaxAvg))

        if np.abs(xMaxAvg-self.oldX) > self.settings['xy_displacement_limit'] or np.abs(yMaxAvg-self.oldY) > self.settings['xy_displacement_limit']:
            xMaxAvg = self.oldX; yMaxAvg = self.oldY

        self.set_coordinate_fnc(xMaxAvg,yMaxAvg,self.oldZ)
        self.oldX = xMaxAvg
        self.oldY = yMaxAvg
    
    def optimize_xy_fast(self):
        xMaxArr = []; yMaxArr = []; cMaxArr = []
        for i in range(1):
            xMax, yMax, most_recent_data_array = self.run_xy_scan_fnc()
            xMaxArr.append(xMax); yMaxArr.append(yMax); cMaxArr.append(most_recent_data_array.max())
        xMaxArr = np.array(xMaxArr); yMaxArr = np.array(yMaxArr); cMaxArr = np.array(cMaxArr)
        xMaxAvg = np.round(np.average(xMaxArr),3); yMaxAvg = np.round(np.average(yMaxArr),3); cMaxAvg = np.round(np.average(cMaxArr),3)
        print('xMaxAvg: ' + str(xMaxAvg))
        print('yMaxAvg: ' + str(yMaxAvg))
        print('cMaxAvg: ' + str(cMaxAvg))

        if np.abs(xMaxAvg-self.oldX) > self.settings['xy_displacement_limit'] or np.abs(yMaxAvg-self.oldY) > self.settings['xy_displacement_limit']:
            xMaxAvg = self.oldX; yMaxAvg = self.oldY

        self.set_coordinate_fnc(xMaxAvg,yMaxAvg,self.oldZ)
        self.oldX = xMaxAvg
        self.oldY = yMaxAvg




    def run_xz_scan_fnc(self):
        self.isStopped = 0
        
        scan_galvo_card_name = "cDAQ1Mod2" # naming the instrument
        scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)} # dictionary of analog output channels
        scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels) # defining the instrument (ni_9263)
        self.scanType = 'xz'

        ############################################################################### def other variables #####################################################################
        ################### setting variales and array ####################

        # acquisition and settling time:
        xz_scan_read_time      = self.settings['xy_scan_read_time'];      xz_scan_settle_time    = self.settings['xy_scan_settle_time']
        xz_scan_resolution_hor = self.settings['xz_scan_resolution_hor']; xz_scan_resolution_ver = self.settings['xz_scan_resolution_ver']
        
        self.time_acquire = np.round(float(xz_scan_read_time)/1e3, 4)
        self.time_settle = np.round(float(xz_scan_settle_time)/1e3, 5)
        self.time_pts_settle = int(self.time_settle*self.samp_rate)
        self.time_pts_acquire = int(self.time_acquire*self.samp_rate)
        self.time_pts_per_pixel = self.time_pts_acquire  + self.time_pts_settle
        # print('self.time_pts_settle: ' + str(self.time_pts_settle))
        # print('self.time_pts_acquire: ' + str(self.time_pts_acquire))
        
        # resolution
        self.hor_res = grid_size_x = int(xz_scan_resolution_hor)
        self.ver_res = grid_size_z = int(xz_scan_resolution_ver)

        # initial driving voltages for the x,y-mirror and z piezo stage
        self.hor_init = x_init = self.oldX - self.settings['x_minus_range']
        self.ver_init = z_init = self.oldZ - self.settings['z_minus_range']
        self.y = y_init = self.oldY

        # final driving voltages for the x,y-mirror
        x_final = self.oldX + self.settings['x_plus_range']
        z_final = self.oldZ + self.settings['z_plus_range']

        # range of voltage
        self.hor_range = x_final - x_init
        self.ver_range = z_final - z_init

        # array of x,z voltages
        self.x_array = np.linspace(x_init, x_final, grid_size_x); self.x_array_original = self.x_array
        self.z_array = np.linspace(z_init, z_final, grid_size_z)
        X, Z = np.meshgrid(self.x_array, self.z_array)

        # make the x-voltage waveform to pass to aotask
        self.x_array = np.array(np.repeat(self.x_array, self.time_pts_per_pixel))
        # print("len(self.x_array): " + str(len(self.x_array)))

        # create dataset to populate
        global xz_scan_data_array
        xz_scan_data_array = np.zeros((grid_size_z, grid_size_x))
        global most_recent_data_array
        most_recent_data_array = np.zeros((grid_size_z, grid_size_x))

        ################### resetting position of mirrors ####################
        scan_galvo.voltage_cdaq1mod2ao0(x_init) # x-mirror
        scan_galvo.voltage_cdaq1mod2ao1(y_init) # y-mirror
        scan_galvo.voltage_cdaq1mod2ao2(z_init) # z-stage

        # the clock to sync AO and counter
        clktask = nidaqmx.Task() 
        clktask.co_channels.add_co_pulse_chan_freq(  # adding dig pulse train chan
            counter = "cDAQ1Mod1/ctr1",
            name_to_assign_to_channel = "",
            units = nidaqmx.constants.FrequencyUnits.HZ,
            idle_state = nidaqmx.constants.Level.LOW,
            initial_delay = 0.0,
            freq = self.samp_rate,
            duty_cycle = 0.5
            )
        clktask.timing.cfg_implicit_timing( # implicit timing by the hardware
            sample_mode = AcquisitionType.CONTINUOUS, # the clock should run continuously in principle
            samps_per_chan = int(len(self.x_array) + 1) # does this matter?
            ) 
        
                    
        ####################################################################### x and z scanning #########################################################################
        print('----------------------------------------------------------------')
        print('XZ scan started')
        for f in range(grid_size_z): # this loops for each row z
            if self.isStopped == 1: break
            start0 = time.time()

            # set initial z location
            scan_galvo.voltage_cdaq1mod2ao2(self.z_array[f])
            time.sleep(100*self.time_settle)
            
            # snake pattern: scanning right to left
            if self.isSnake == 1 and f % 2 == 1: self.x_array_write = np.ascontiguousarray(self.x_array[::-1])
            else: self.x_array_write = self.x_array

            ##############################################################
            # counter task
            ctrtask = nidaqmx.Task()
            ctrtask.ci_channels.add_ci_count_edges_chan( # define the counter
                counter = "cDAQ1Mod1/ctr0",
                name_to_assign_to_channel = "",
                edge = nidaqmx.constants.Edge.RISING,
                initial_count = 0,
                count_direction = nidaqmx.constants.CountDirection.COUNT_UP
                )
            ctrtask.timing.cfg_samp_clk_timing( # cfg sample clk timing
                rate = self.samp_rate,
                source = "/cDAQ1/Ctr1InternalOutput", # the clock defined above
                active_edge = nidaqmx.constants.Edge.RISING,
                sample_mode = AcquisitionType.FINITE,
                samps_per_chan = int(len(self.x_array) + 1) 
                )

            ##############################################################
            # AO task for the x-galvo
            aotask = nidaqmx.Task() 
            channel_galvo_x = f'{scan_galvo_card_name}/ao{0}'
            aotask.ao_channels.add_ao_voltage_chan(channel_galvo_x)
            aotask.timing.cfg_samp_clk_timing(
                rate = self.samp_rate,
                source = "/cDAQ1/Ctr1InternalOutput", # the clock defined above
                active_edge = nidaqmx.constants.Edge.RISING,
                sample_mode = AcquisitionType.FINITE,
                samps_per_chan = len(self.x_array)
                )
            aotask.write(self.x_array_write, auto_start=False) # assign the x-waveform to aotask, but not start yet

            # start galvo X scan, counter, and clock
            aotask.start()
            ctrtask.start()
            clktask.start() # the clock MUST START AFTER aotask.start() and ctrtask.start()
            
            aotask.wait_until_done()
            xLineData = ctrtask.read(len(self.x_array) + 1)  # +1 is to take the difference later
            # print("len(xLineData): " + str(len(xLineData)))

            aotask.stop(); aotask.close()
            ctrtask.stop(); ctrtask.close()
            clktask.stop() 
            diffData = np.diff(xLineData)

            summedData = np.zeros(grid_size_x)
            for i in range(0, grid_size_x):
                summedData[i] = np.sum(
                    diffData[(i*self.time_pts_per_pixel + self.time_pts_settle):(i*self.time_pts_per_pixel + self.time_pts_per_pixel - 1)])
            normedData = summedData * (.001 / self.time_acquire) # normalize to kcounts/sec

            if self.isSnake == 1 and f % 2 == 1: xz_scan_data_array[f] = np.flip(normedData)
            else: xz_scan_data_array[f] = normedData

            end = time.time() - start0
            # if f == 0 or f == grid_size_z-1 or f == int(grid_size_z/2):
                # print("Time per row: " + str(end) + " s")

        clktask.close(); scan_galvo.close() 

        ############################################ end XZ scanning script #################################################
        ############################################### plotting XZ scan data in plot ####################################################
        most_recent_data_array = xz_scan_data_array

        indices = np.where(most_recent_data_array == most_recent_data_array.max())
        xMax_idx = indices[1][0]; zMax_idx = indices[0][0]
        xMax = self.x_array_original[xMax_idx]; zMax = self.z_array[zMax_idx]
        # print('xMax: ' + str(xMax))
        # print('zMax: ' + str(zMax))
        print('XZ scan finished')
        # print('------------------------------------------------------------------')
        # plt.imshow(most_recent_data_array)
        # plt.show()
        return xMax, zMax, most_recent_data_array

    def optimize_xz(self):
        xMaxArr = []; zMaxArr = []
        for i in range(self.settings['num_of_scans']):
            xMax, zMax, most_recent_data_array = self.run_xz_scan_fnc()
            xMaxArr.append(xMax); zMaxArr.append(zMax)
        xMaxArr = np.array(xMaxArr); zMaxArr = np.array(zMaxArr)
        xMaxAvg = np.round(np.average(xMaxArr),3); zMaxAvg = np.round(np.average(zMaxArr),3)
        print('xMaxAvg: ' + str(xMaxAvg))
        print('zMaxAvg: ' + str(zMaxAvg))

        if np.abs(xMaxAvg-self.oldX) > self.settings['xz_displacement_limit'] or np.abs(zMaxAvg-self.oldZ) > self.settings['xz_displacement_limit']:
            xMaxAvg = self.oldX; zMaxAvg = self.oldZ

        self.set_coordinate_fnc(self.oldX,self.oldY,zMaxAvg)
        self.oldZ = zMaxAvg
    
    def set_coordinate_fnc(self, x,y,z):
        # naming the instrument
        scan_galvo_card_name = "cDAQ1Mod2"

        # dictionary of analog output channels
        scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)}

        # defining the instrument (ni_9263)
        scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels)
        
        scan_galvo.voltage_cdaq1mod2ao0(x) # this is for the x-mirror
        scan_galvo.voltage_cdaq1mod2ao1(y) # this is for the y-mirror
        scan_galvo.voltage_cdaq1mod2ao2(z) # this is for the z-piezo

        # self.spot, = self.sc.axes.plot(self.eventXdata, self.eventYdata, marker="o", markersize=3, markeredgecolor="red", markerfacecolor="red")
        print('Set coordinates to (' + str(x) + ', ' + str(y) + ', ' + str(z) + ')') # EP 10/16

        with open('C:/Users/lukin2dmaterials/miniconda3/envs/NV_control/B00_codes/NVlocation.txt', 'w') as f:
            f.write(str(x) + '\n')
            f.write(str(y) + '\n')
            f.write(str(z))
            f.close()
        scan_galvo.close()  

    def close(self):
        pb.close()

if __name__ == '__main__':
    print()
    cfcObject = Confocal()
    cfcObject.optimize_xy()
    