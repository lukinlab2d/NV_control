"""
B00 gui for room-temp NV center measurement experiment control

Author: Miles D. Ackerman (undergraduate during the summer of 2022). Email: miles.ackerman1@gmail.com. Date: 070622 -> 090222

General info:
This file makes a GUI to control the room-temp NV center-based measurement setup in LISE room B00 in the Lukin Group in the Harvard PHY dept. The setup
is a confocal laser scanning microscope. A galvo steers a laser beam at 518 nm (butteryfly) through a 4F telescope setup. An excitation and a collection
path (both use fiber coupling) are present. The sample stage is manually-controlled by positioning dials. Additionally, the z-axis of the microscope
objective is controlled by a piezo. 

This entire setup is controlled by a NI-cDAQ (model 9147). The NI-DAQmx is the National Insturments API that is used to communicate with the cDAQ device. 
Info about this device can be found at: https://nidaqmx-python.readthedocs.io/en/latest/https://nidaqmx-python.readthedocs.io/en/latest/. Two cards/modules 
are present in the cDAQ chassis: 

1. NI-9402 (four bi-directional digital input/output (I/O) channels -also featuring programmable internal clocks for counting). Info about the NI-9402 
module can be found at: https://www.ni.com/en-us/support/model.ni-9402.htmlhttps://www.ni.com/en-us/support/model.ni-9402.html. This module functions 
as a counter based on an internal clock within the module. The parameters for creating a counter assigned to this module when it is created in each 
scanning script. The counter requires the use of two of the four available channels on this cDAQ card. One channel is for input (from the APD on the
optical bench), and one channel for the hardware-timed clock to "run" the counter. Currently (as of 081522) only one counter (analog input) channel 
can be created at a time. In order to change this functionality, the source code for this contributed drivers package must be edited. 

2. NI-9263 (four analog output channels (-10 V to +10 V). This module's info can be found at: 
https://www.ni.com/en-us/support/model.ni-9263.htmlhttps://www.ni.com/en-us/support/model.ni-9263.html. This module is used for programmable digital 
output. This functionality is largely used for controlling the servo motors on the ThorLab's galvo (by writing a sequence of voltages to them). Here,
all available digital output channels are created in one line: `ni_9263_ao_channels` where a single channel can be accessed at a time. The main 
function of this application is running confocal scanning scripts in different planes (the XY plane at a specified Z height, the XZ plane across a 
specified z-range, and the YZ plane at a specified X range).

The GUI is split into 2 main setions, a left-half and a right-half section. The left-hand section contains the input controls for running a specific scanning 
script. On the bottom of the left-hand section is an option to save the dataset that was just scan and plotted on the right. Additionally, a text output box 
is present at the bottom of the left_window section. This is where the user-entered scanning parameters are printed to. Here, the use has the option to select 
and then copy the scan parameter info for keeping a record. there is also the option (currently commented out) to print this same info w/ its same format to 
the terminal. The right-hand section contains the output plot from a run scan. Include info more: this application implements live plotting for each 
scanning option.

For the user of this application:
1. This application was built within a virtual environment on this computer called: "NV_control" it can be found through the
user/miniconda3/envs/NV_control. Within this enviornment THIS file and additonal reference packages can be found
2. Only limited input-validation/error checking has been implemented so far. Please be advised when running scanning programs. The LIMIT of the voltage 
that can be driven to either mirror (or the objective's z_piezo) is 10 V. DO NOT EXCEED THIS VALUE. The system (as a cause of NI-DAQmx error 
checking) will likely refuse such a request, but it is best practice to avoid.

Helpful info. when reading through this application code:
1. The NI-cDAQ card (chasis) is called (this can be changed using the NI-MAX software when the NI-DAQ card is connected to this computer)
"cDAQ1". Additional info. about the cDAQ device and its inserted modules (as of 081522 only ni_9402 and ni_9263) can be viewed in the NI-MAX software.
This software opens automatically when the cDAQ device is connected to this computer
2. The first module inserted into the cDAQ card (the NI-9402 module) is called "mod1". This naming pattern follows the available slots within the 
cDAQ chasis (there are four available slots)
3. The second module in the cDAQ chasis (the NI-9263 module) is called "mod2" following in suit
"""

# TODO: remove unused imports
# TODO: update naming conventions throughout scripts and entire file
# TODO: can about window be set to the center of the screen regardles of the position of the main window?
# TODO: fix allowing scans to be run one after the other. Will "Try... Except" work?
# 090122 this above is fixed. However, the fix lies in the implementation of validating the resolution input. If the input-validation framework changes 
# (as it should bc it is limited to only checking reoslution now) the allowed-repeated-scanning functionality MUST be re-written
# TODO: cont. below
"""
Input validation (resolution, min and max voltages for axes)
Multile scans w/out closing programs (scan_galvo.close) (use Try... Except?) -fixed; not ideal
Saving scan data/images
Progress (scanning) bar
Live plot updating?
Point size minimum based on wavelength of used laser
input validation method is currently only limited to checking resolutions setting per each scan
"""

#################################################################### imports ###################################################################################

# general packages
import sys
import os
from xml.sax.handler import DTDHandler
import nidaqmx
import numpy
import pyqtgraph as pg
import time
from datetime import date
import threading
from B00_codes import TurnOnLaser, TurnOffLaser
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc

# MatPlotLib plotting packages
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import matplotlib as mpl
import matplotlib.pyplot as plt

# QCoDeS -this only imports needed classes from a document not other QCoDeS functionality
from qcodes_contrib_drivers.drivers.NationalInstruments.DAQ import *
from qcodes_contrib_drivers.drivers.NationalInstruments.class_file import *

# PyQt5, this is the framework that builds the GUI
import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5 import QtWidgets
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QMessageBox

############### global variable ##############
any_script_run_one_Q = False # for multiple scanning

############### prelims #######################
get_todays_date = date.today() # this is used for creating the final plot's plot labels
todays_date = get_todays_date.strftime("%m%d%Y") # this is used for creating the final plot's plot labels

################################################################## "Make_Error_Window_2" Class ######################################################################
class Make_Error_Window_2(QtWidgets.QMainWindow): # create the "Make_Error_Window_2" for displaying a new window with error content

    # ?
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setStyleSheet("background-color: red;")

        self.title = "Error" # define the title of the error window

        self.top = 350 # set the display location of the error window
        self.left = 675 # can this be set to the center of the screen regardles of the position of the main window?

        self.error_window_width  = 205 # define the width of the error window
        self.error_window_height = 50 # define the height of the error window

        self.setMaximumSize(self.error_window_width, self.error_window_height) # set the maximum size of the error window
        self.setMinimumSize(self.error_window_width, self.error_window_height) # set the minimum size of the error window

        # begin content of the error window
        error_window_left_justify_adjust = 5 # optional adjustment parameter for the content of the error window (left justify)

        error_window_top_justify_adjust = 5 # optional adjustment parameter for the content of the error window (top justify)

        error_window_content_line_1 = QLabel("ERROR!", self)
        error_window_content_line_1.move(60 + error_window_left_justify_adjust, 0 + error_window_top_justify_adjust)
        error_window_content_line_1.resize(300, 15)

        error_window_content_line_2 = QLabel("Adjust address to save", self)
        error_window_content_line_2.move(40 + error_window_left_justify_adjust, 15 + error_window_top_justify_adjust)
        error_window_content_line_2.resize(300, 15)

        # end content of the error window

        self.setWindowTitle(self.title) # set the title of the displayed error window
        self.setGeometry(self.left, self.top, self.error_window_width, self.error_window_height) # set the geometry (size) of the displayed error window

################################################################## "Make_Error_Window" Class ######################################################################
class Make_Error_Window(QtWidgets.QMainWindow): # create the "Make_Error_Window" for displaying a new window with error content

    # ?
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setStyleSheet("background-color: red;")

        self.title = "Error" # define the title of the error window

        self.top    = 350 # set the display location of the error window
        self.left   = 675 # can this be set to the center of the screen regardles of the position of the main window?

        self.error_window_width  = 205 # define the width of the error window
        self.error_window_height = 50 # define the height of the error window

        self.setMaximumSize(self.error_window_width, self.error_window_height) # set the maximum size of the error window
        self.setMinimumSize(self.error_window_width, self.error_window_height) # set the minimum size of the error window

        # begin content of the aobut window
        error_window_left_justify_adjust = 5 # optional adjustment parameter for the content of the error window (left justify)

        error_window_top_justify_adjust = 5 # optional adjustment parameter for the content of the error window (top justify)

        error_window_content_line_1 = QLabel("ERROR!", self)
        error_window_content_line_1.move(60 + error_window_left_justify_adjust, 0 + error_window_top_justify_adjust)
        error_window_content_line_1.resize(300, 15)

        error_window_content_line_2 = QLabel("Adjust resolution", self)
        error_window_content_line_2.move(45 + error_window_left_justify_adjust, 15 + error_window_top_justify_adjust)
        error_window_content_line_2.resize(300, 15)

        # end content of the error window

        self.setWindowTitle(self.title) # set the title of the displayed error window
        self.setGeometry(self.left, self.top, self.error_window_width, self.error_window_height) # set the geometry (size) of the displayed error window


################################################################## "Make_About_Window" Class ######################################################################
class Make_About_Window(QtWidgets.QMainWindow): # create the "Make_About_Window" for displaying a new window with about content

    # ?
    def __init__(self, parent=None):
        super().__init__(parent)

        self.title = "About" # define the title of the about window

        self.top    = 350 # set the display location of the about window
        self.left   = 675 # can this be set to the center of the screen regardles of the position of the main window?

        self.about_window_width  = 205 # define the width of the about window
        self.about_window_height = 70 # define the height of the about window

        self.setMaximumSize(self.about_window_width, self.about_window_height) # set the maximum size of the about window
        self.setMinimumSize(self.about_window_width, self.about_window_height) # set the minimum size of the about window

        # begin content of the aobut window
        about_window_left_justify_adjust = 5 # optional adjustment parameter for the content of the about window (left justify)

        about_window_top_justify_adjust = 5 # optional adjustment parameter for the content of the about window (top justify)

        about_window_content_line_1 = QLabel("Application name: mda_b00_gui", self)
        about_window_content_line_1.move(0 + about_window_left_justify_adjust, 0 + about_window_top_justify_adjust)
        about_window_content_line_1.resize(300, 15)

        about_window_content_line_2 = QLabel("Author: Miles D. Ackerman", self)
        about_window_content_line_2.move(0 + about_window_left_justify_adjust, 15 + about_window_top_justify_adjust)
        about_window_content_line_2.resize(300, 15)

        about_window_content_line_3 = QLabel("Last modified: 090222", self)
        about_window_content_line_3.move(0 + about_window_left_justify_adjust, 30 + about_window_top_justify_adjust)
        about_window_content_line_3.resize(300, 15)

        about_window_content_line_3 = QLabel("OS: MS Windows 10 Pro", self)
        about_window_content_line_3.move(0 + about_window_left_justify_adjust, 45 + about_window_top_justify_adjust)
        about_window_content_line_3.resize(300, 15)
        # end content of the about window

        self.setWindowTitle(self.title) # set the title of the displayed about window
        self.setGeometry(self.left, self.top, self.about_window_width, self.about_window_height) # set the geometry (size) of the displayed about window

################################################### MatPlotLib class ######################################################################################
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent = None, width = 10, height = 10, dpi = 1000):

        # fig = Figure(figsize = (width, height), dpi = dpi)
        # self.axes = fig.add_subplot()
        self.fig, self.axes = plt.subplots(figsize=(width, height), dpi=dpi, tight_layout = True)
        super(MplCanvas, self).__init__(self.fig)

########################################################################## "Parent" class #####################################################################
class Parent(QtWidgets.QMainWindow):

    # ?
    def __init__(self, parent = None):
        super().__init__(parent)

        def display_about_window():

            self.Make_About_Window = Make_About_Window()
            self.Make_About_Window.show()

        ######################################################################### GUI prelims ##############################################################################

        # setting up main GUI window
        self.child_widget = Child(parent = self)
        self.setCentralWidget(self.child_widget) # setting the central widget of the main window (Parent class) to the Child class
        gui_window_height = 470 # define the main window height
        gui_window_width = 800 # define the main window width
        self.setGeometry(400, 200, gui_window_width, gui_window_height) # x-coord, y-coord, width, height
        self.setMinimumSize(gui_window_width, gui_window_height) # set the main window min size
        self.setMaximumSize(gui_window_width, gui_window_height) # set the main window max size
        self.setWindowTitle("mda_b00_gui") # set the title of the main window

        ################################################################### menu bar #############################################################################
        main_window_menu_bar = self.menuBar() # this creates the menu abr for the main GUI window
        
        # "File" menu option
        file_menu = main_window_menu_bar.addMenu("File") # this adds the "File" option to the main window's menu bar
        exit_option = QtWidgets.QAction("Exit", self) # adds the "Exit" sub_option to "File" menu option
        exit_option.triggered.connect(qApp.quit) # setting the fnc of clicking "Exit" sub_option to quit application
        file_menu.addAction(exit_option) # adding "Exit" sub_option to "File" option
        
        # "Help" menu option
        help_menu = main_window_menu_bar.addMenu("Help") # this adds the "Help" option to the main window's menu bar
        hep_menu_about = QtWidgets.QAction("About", self) # this adds an "About" sub_option to the "Help" option
        hep_menu_about.triggered.connect(display_about_window) # this connects clicking the "About" to "..."
        help_menu.addAction(hep_menu_about) # adding "About" sub_option to the "Help" menu option

############################################################################# "Child" class #######################################################################
class Child(QtWidgets.QWidget):#, **kwargs): # kwargs needed?

    # ?
    def __init__(self, parent = None):#, **kwargs): # kwargs needed?
        super().__init__(parent)
        global channel_currently_on; channel_currently_on = ()
        hbox = QHBoxLayout(self)
        self.samp_rate = 1e5  # 1e3 (old)
        self.hor_coord = 0; self.ver_coord = 0
        self.isStopped = 0; self.isSnake = 0

        ########################################################## CHILD functions ################################################################################
        
        ########## overall #################

        # display invalid resolution error window fnc
        def display_resolution_error_window_fnc(): # this fnc calls the "Make_Error_Window" class to display an eror message indicating user input is not validated
            self.Make_Error_Window = Make_Error_Window()
            self.Make_Error_Window.show()

        # display invalid saving address error window fnc   
        def display_save_address_length_error_window_fnc(): # this fnc calls the "Make_Error_Window" class to display an eror message indicating user input is not validated
            self.Make_Error_Window_2 = Make_Error_Window_2()
            self.Make_Error_Window_2.show()

        # save most recent scan data fnc
        def save_scan_data_fnc(): # this fnc works for any scanning script

            """
            How this works/applies to each scanning script:
            In any scanning script (XY, XZ, and YZ), a data_array is created according to the user-specified grid_size. It is a Numpy array of zeros that will be
            populated throughout the scanning program as it progresses. At the same time of that data-array created a global variable called "most_recent_data_array"
            is created and is then set to the same size matching the scan-specific data_array. At the end of the scanning scipt this temporary data array is matched to
            the scan's specific data array, value for value. Now "most_recent-data_array" is called (since it is defined to be Global) below for saving. This
            """

            saving_scan_error_bool = False # setting up a bool value for error checking below

            # print("save_scan_data_fnc called")                                               # delete later
            # print("@address :" + save_scan_data_qlineedit.text())                                               # delete later

            # while loop for error checking if address to save data at has length > 0
            while saving_scan_error_bool is False:

                if len(str(save_scan_data_qlineedit.text())) == 0: # checking if length of specified saving address is > 0

                    # print("EXCEPTION!")                                                       # safe to delete
                    display_save_address_length_error_window_fnc()                              # fix to new window. Only a test now
                    break # condition remains False; not saved; exit

                elif len(str(save_scan_data_qlineedit.text())) > 0: # checking the length of the specified address is greater than 0

                    saving_scan_error_bool == True # adjusting the value of the current bool to True
                    address_to_save_scan_data_at = save_scan_data_qlineedit.text() # creating a variable as the specified (now error-checked) address
                    numpy.save(str(address_to_save_scan_data_at), most_recent_data_array) # saving the correct data array
                    break # data has been successfully save; so exit checking loop

        def read_PL(a, b):
            # Calculate the distance between each element and the target (a, b)
            hor_dist = np.abs(self.hor_array - a)
            ver_dist = np.abs(self.ver_array - b)

            # Find the index of the element with the minimum distance
            closest_hor_idx = np.argmin(hor_dist)
            closest_ver_idx = np.argmin(ver_dist)

            # Return the closest element
            PL = self.result[closest_ver_idx, closest_hor_idx]
            return PL
                    
        def mouse_event(event):
            self.hor_coord = np.round(event.xdata,3)
            self.ver_coord = np.round(event.ydata,3)
            self.crosshair.remove()
            self.crosshair, = self.sc.axes.plot(event.xdata, event.ydata, marker="+", markersize=5, markeredgecolor="red", markerfacecolor="red")
            self.sc.draw()
            # print('hor: {} and ver: {}'.format(self.hor_coord, self.ver_coord))

            if self.scanType == 'xy':
                x_qlineedit.setText(str(self.hor_coord))
                y_qlineedit.setText(str(self.ver_coord))
                z_qlineedit.setText(str(self.z))
            elif self.scanType == 'xz':
                x_qlineedit.setText(str(self.hor_coord))
                y_qlineedit.setText(str(self.y))
                z_qlineedit.setText(str(self.ver_coord))
            elif self.scanType == 'yz':
                x_qlineedit.setText(str(self.x))
                y_qlineedit.setText(str(self.hor_coord))
                z_qlineedit.setText(str(self.ver_coord))
            
            PL = read_PL(self.hor_coord, self.ver_coord)
            self.parameters_display_text_box.clear()
            self.parameters_display_text_box.setPlainText(
                                                        'hor: {} and ver: {}'.format(self.hor_coord, self.ver_coord)  + "\n"
                                                        "PL count rate = " + str(PL)
                                                        )

        def set_coordinate_fnc():
            # naming the instrument
            scan_galvo_card_name = "cDAQ1Mod2"

            # dictionary of analog output channels
            scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)}

            # defining the instrument (ni_9263)
            scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels)
            
            x = float(x_qlineedit.text())
            y = float(y_qlineedit.text())
            z = float(z_qlineedit.text())

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

        def stop_fnc():
            self.isStopped = 1
        
        def snake_fnc():
            self.isSnake = int(snake_checkbox.isChecked())
        
        def toggle_laser_fnc_532():
            if toggle_laser_checkbox_532.isChecked():
                global pb_532; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(3,)))
                pb_532 = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_532")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==3)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_532.close()
                pb_532_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_532_off")
                pb_532_off.close() # close the instrument
                print('532 nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")

        def toggle_laser_fnc_589():
            if toggle_laser_checkbox_589.isChecked():
                global pb_589; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(6,)))
                pb_589 = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_589")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==6)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_589.close()
                pb_589_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_589_off")
                pb_589_off.close() # close the instrument
                print('589 nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
        
        def toggle_laser_fnc_589w():
            if toggle_laser_checkbox_589w.isChecked():
                global pb_589w; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(9,)))
                pb_589w = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_589w")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==9)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_589w.close()
                pb_589w_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_589w_off")
                pb_589w_off.close() # close the instrument
                print('589w nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")

        def toggle_laser_fnc_vel():
            if toggle_laser_checkbox_vel.isChecked():
                global pb_vel; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(5,)))
                pb_vel = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_vel")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==5)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_vel.close()
                pb_vel_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_vel_off")
                pb_vel_off.close() # close the instrument
                print('vel nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
        
        def toggle_laser_fnc_pt637():
            if toggle_laser_checkbox_pt637.isChecked():
                global pb_pt637; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(8,)))
                pb_pt637 = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_pt637")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==8)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_pt637.close()
                pb_pt637_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_pt637_off")
                pb_pt637_off.close() # close the instrument
                print('pt637 nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
        
        def toggle_laser_fnc_532b():
            if toggle_laser_checkbox_532b.isChecked():
                global pb_532b; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(7,)))
                pb_532b = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_532b")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==7)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_532b.close()
                pb_532b_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_532b_off")
                pb_532b_off.close() # close the instrument
                print('532b nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
        
        def toggle_laser_fnc_pt637b():
            if toggle_laser_checkbox_pt637b.isChecked():
                global pb_pt637b; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(10,)))
                pb_pt637b = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_pt637b")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==10)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_pt637b.close()
                pb_pt637b_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_pt637b_off")
                pb_pt637b_off.close() # close the instrument
                print('pt637b nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")

        def toggle_laser_fnc_velb():
            if toggle_laser_checkbox_velb.isChecked():
                global pb_velb; global channel_currently_on
                channel_currently_on = np.concatenate((channel_currently_on,(14,)))
                pb_velb = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_velb")
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
            else:
                channel_off_idx = np.where(channel_currently_on==14)[0][0]
                channel_currently_on = np.delete(channel_currently_on, channel_off_idx)
                pb_velb.close()
                pb_velb_off = TurnOnLaser.turnOnLaser(channels=channel_currently_on, instrument_name="PB_velb_off")
                pb_velb_off.close() # close the instrument
                print('velb nm laser turned off')
                print("Channels currently on: " + str(channel_currently_on))
                print("-------------------------------------")
                
        ########### XY scanning #############
        # print/display XY scan parameters fnc
        def print_XY_scan_parameters_fnc(self, parent = Child): # this fnc does...

            # this prints to the QTextBox in the left_window. The output of the user-selected scan parameters is printed below
            self.parameters_display_text_box.clear()
            self.parameters_display_text_box.setPlainText(
                                                        "XY_SCAN PARAMETERS/INFO:\n"
                                                        "XY_scan resolution = " + str(int(xy_scan_resolution_hor_qlineedit.text())) + "\n"
                                                        "XY_scan counter read time = " + str(float(xy_scan_read_time_qlineedit.text())/1e3) + "\n"
                                                        "XY_scan min x driving voltage = " + str(float(xy_scan_x_voltage_min_qlineedit.text())) + "\n"
                                                        "XY_scan max x driving voltage = " + str(float(xy_scan_x_voltage_max_qlineedit.text())) + "\n"
                                                        "XY_scan min y driving voltage = " + str(float(xy_scan_y_voltage_min_qlineedit.text())) + "\n"
                                                        "XY_scan max y driving voltage = " + str(float(xy_scan_y_voltage_max_qlineedit.text())) + "\n"
                                                        "XY_scan z-piezo driving voltage = " + str(float(xy_scan_z_voltage_qlineedit.text()))
                                                        )
        
        # run XY scanning function
        def run_xy_scan_fnc(): # this fnc runs the xy_scan per the user-entered parameters in the xy_scan qlineedits
            self.isStopped = 0
            """
            This runs X and Y only scan. It currently creates and then populates a user defined size numpy array according to a set counter acquisition time and a motor
step voltage setting. Additionally, the initial driving voltage for the X and Y motors can be set according to the desired scanning range. This scanning program runs in a snake pattern, it scan the first row left to right, moves up one row, 
then scans right to left and continues. Alternatives would be scanning left to right, and resetting the position of the laser on the next higher row and scanning again left 
to right OR scanning in a "circular" patter either CW or CCW from outside to inside or inside to outside. The chosen method was picked for simplicity of understanding. The 
scanning loops are present within NI-DAQmx tasks to speed up the program. Starting and stopping a NI-DAQmx task repeatedly slows down the program dramatically. So, the 
counter and hardware clock task are started once, then the scanning program is run, and then the counter and clock tasks are closed -un-reserving the hardware resources. 
This cell uses the "DAQAnalogOutputs" function from a written class file at:
C:/Users/lukin2dmaterials/miniconda3/envs/qcodes/Lib/site-packages/qcodes_contrib_drivers/drivers/NationalInstruments/class_file. Slashes are reversed to run
            """

            ############################################################### begin scanning script #############################################################################################
            ################################################################################## card 2 (AO) ########################################################################
            try:
                scan_galvo_card_name = "cDAQ1Mod2" # naming the instrument
                scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)} # dictionary of analog output channels
                scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels) # defining the instrument (ni_9263)

                get_todays_date = date.today() # this is used for creating the final plot's plot labels
                todays_date = get_todays_date.strftime("%m%d%Y") # this is used for creating the final plot's plot labels
                self.scanType = 'xy'

                ############################################################################### def other variables #####################################################################
                ################### setting variales and array ####################

                # acquisition and settling time:
                self.time_acquire = np.round(float(xy_scan_read_time_qlineedit.text())/1e3, 4)
                self.time_settle = np.round(float(xy_scan_settle_time_qlineedit.text())/1e3, 5)
                self.time_pts_settle = int(self.time_settle*self.samp_rate)
                self.time_pts_acquire = int(self.time_acquire*self.samp_rate)
                self.time_pts_per_pixel = self.time_pts_acquire  + self.time_pts_settle
                print('self.time_pts_settle: ' + str(self.time_pts_settle))
                print('self.time_pts_acquire: ' + str(self.time_pts_acquire))

                # resolution
                self.hor_res = grid_size_x = int(xy_scan_resolution_hor_qlineedit.text())
                self.ver_res = grid_size_y = int(xy_scan_resolution_ver_qlineedit.text())

                # initial driving voltages for the x,y-mirror and z piezo stage
                self.hor_init = x_init = round(float(xy_scan_x_voltage_min_qlineedit.text()), 3)
                self.ver_init = y_init = round(float(xy_scan_y_voltage_min_qlineedit.text()), 3)
                self.z = z_init = round(float(xy_scan_z_voltage_qlineedit.text()), 3)

                # final driving voltages for the x,y-mirror
                x_final = round(float(xy_scan_x_voltage_max_qlineedit.text()), 3)
                y_final = round(float(xy_scan_y_voltage_max_qlineedit.text()), 3)

                # range of voltage
                self.hor_range = x_final - x_init
                self.ver_range = y_final - y_init

                # array of x,y voltages
                self.x_array = np.linspace(x_init, x_final, grid_size_x)
                self.y_array = np.linspace(y_init, y_final, grid_size_y)
                
                # correct for snake shift
                if snake_shift_qlineedit.text() == '': snake_shift = 0
                else: snake_shift = float(snake_shift_qlineedit.text())
                self.x_array_plot = np.linspace(x_init + snake_shift, x_final + snake_shift, grid_size_x)
                self.y_array_plot = np.linspace(y_init, y_final, grid_size_y)
                X_plot, Y_plot = np.meshgrid(self.x_array_plot, self.y_array_plot)

                self.hor_array = self.x_array_plot
                self.ver_array = self.y_array_plot

                # make the x-voltage waveform to pass to aotask
                self.x_array = np.array(np.repeat(self.x_array, self.time_pts_per_pixel))
                print("len(self.x_array): " + str(len(self.x_array)))

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
                        source = "/cDAQ1/Ctr1InternalOutput", # the clktask clock defined above
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
                    if f == 0 or f == grid_size_y-1 or f == int(grid_size_y/2):
                        print("Time per row: " + str(end) + " s")

                    self.sc.axes.pcolormesh(X_plot, Y_plot, xy_scan_data_array, cmap = "inferno")
                    self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
                    self.sc.axes.set_ylabel("y_mirror_driving_voltage_(V)", fontsize = 8, labelpad = -9)
                    self.sc.fig.canvas.draw()
                    self.sc.fig.canvas.flush_events() # this line is very important
                        
                clktask.close(); scan_galvo.close() 

                ############################################################### end scanning script #############################################################################################
                ############################################### plotting XY scan data in plot ###################################################
                self.sc.axes.cla()
                self.result = xy_scan_data_array

                plot = self.sc.axes.pcolormesh(X_plot, Y_plot, xy_scan_data_array, cmap = "inferno")
                self.xy_scan_plot_colorbar = self.sc.fig.colorbar(plot, ax = self.sc.axes, pad = 0.02, aspect = 15)
                self.xy_scan_plot_colorbar.formatter.set_powerlimits((0, 0))
                self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                self.xy_scan_plot_colorbar.ax.tick_params(labelsize = 7)
                self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
                self.sc.axes.set_ylabel("y_mirror_driving_voltage_(V)", fontsize = 8, labelpad = -9)
                self.sc.axes.set_title("XY_scan_%s_z-piezo@%d_microns" % (todays_date, int((z_init * 10))), fontsize = 8)
                self.crosshair, = self.sc.axes.plot(self.x_array[0], self.y_array[0], marker="+", markersize=5, markeredgecolor="red", markerfacecolor="red")
                self.sc.draw()

                cid = self.sc.mpl_connect('button_press_event', mouse_event)
                most_recent_data_array = xy_scan_data_array # make temp holding global data_array the same as xy_scan_data_array

                print('XY scan finished')
                print('----------------------------------------------------------------')
            except nidaqmx.errors.DaqError:
                print('aaaaaa')
                scan_galvo.close()
                aotask.stop(); aotask.close()
                ctrtask.stop(); ctrtask.close()
                clktask.stop(); clktask.close()  

        # xy_scan resolution check then run fnc
        def xy_scan_resolution_validation_fnc():
            self.sc.axes.cla()

            # the try-except frameworks are used to refresh the plotted figs -removing the color bars associated with a previous plotted data
            try:
                self.xy_scan_plot_colorbar.remove()
            except (AttributeError, ValueError):
                pass
        
            try:
                self.yz_scan_plot_colorbar.remove()
            except (AttributeError, ValueError):
                pass

            try:
                self.xz_scan_plot_colorbar.remove()
            except (AttributeError, ValueError):
                pass
        
        #################################### resolution checking ##################################

            res_min_condition = 0 # set the min allowed resolution for scanning
            res_max_condition = 900 # set the max allowed resolution for scanning
            xy_scan_resolution_test_condition = False # define resolution validation bool for xy scan

            while xy_scan_resolution_test_condition is False: # this initiates checking the resolution parameter

                # checking for out of bounds of min and max conditions above
                # TODO: or negative or not a number or too large
                if int(xy_scan_resolution_hor_qlineedit.text()) < res_min_condition or int(xy_scan_resolution_hor_qlineedit.text()) > res_max_condition:
                    display_resolution_error_window_fnc() # call the error message pop-up window
                    break # exit the checking loop: failed

                # if parameter is in bounds; run scan
                elif int(xy_scan_resolution_hor_qlineedit.text()) >= res_min_condition and int(xy_scan_resolution_hor_qlineedit.text()) <= res_max_condition:
                    xy_scan_resolution_test_condition == True
                    print_XY_scan_parameters_fnc(self) # call the print user-entered parameters fnc
                    run_xy_scan_fnc() # call the run xy scan method fnc
                    break # exit the checking loop: passed

        ##########################################################################################################################################################################################################
        ##########################################################################################################################################################################################################
        ########### XZ scanning #############
        # print_XZ_scan_parameters_fnc
        def print_XZ_scan_parameters_fnc(self, parent = Child): # this fnc does...

            # print("XZ_SCAN PARAMETERS/INFO: ", end = "")                                    # this prints to the terminal
            # print("XZ_scan resolution = %d, " % int(xz_scan_resolution_hor_qlineedit.text()), end = "")
            # print("XZ_scan counter read time = %2f, " % round(float(xz_scan_read_time_qlineedit.text())/1e3, 2), end = "")
            # print("XZ_scan min x driving voltage = %2f, " % float(xz_scan_x_voltage_min_qlineedit.text()), end = "")
            # print("XZ_scan max x driving voltage = %2f, " % float(xz_scan_x_voltage_max_qlineedit.text()), end = "")
            # print("XZ_scan y driving voltage = %2f, " % float(xz_scan_y_voltage_qlineedit.text()), end = "")
            # print("XZ_scan z-piezo min driving voltage = %2f, " % float(xz_scan_z_voltage_min_qlineedit.text()), end = "")
            # print("XZ_scan z-piezo max driving voltage = %2f." % float(xz_scan_z_voltage_max_qlineedit.text()))


            # need to clear text box first
            self.parameters_display_text_box.clear()

            # this prints to the QTextBox in the left_window. The output of the user-selected scan parameters is printed below
            self.parameters_display_text_box.setPlainText(
                                                        "XZ_SCAN PARAMETERS/INFO:\n"
                                                        "XZ_scan resolution = " + str(int(xz_scan_resolution_hor_qlineedit.text())) + "\n"
                                                        "XZ_scan counter read time = " + str(float(xz_scan_read_time_qlineedit.text())/1e3) + "\n"
                                                        "XZ_scan min x driving voltage = " + str(float(xz_scan_x_voltage_min_qlineedit.text())) + "\n"
                                                        "XZ_scan max x driving voltage = " + str(float(xz_scan_x_voltage_max_qlineedit.text())) + "\n"
                                                        "XZ_scan y driving voltage = " + str(float(xz_scan_y_voltage_qlineedit.text())) + "\n"
                                                        "XZ_scan z-piezo min driving voltage = " + str(float(xz_scan_z_voltage_min_qlineedit.text())) + "\n"
                                                        "XZ_scan z-piezo max driving voltage = " + str(float(xz_scan_z_voltage_max_qlineedit.text()))
                                                        )
        
        # run XZ scan fnc
        def run_xz_scan_fnc(): # this fnc runs the xz_scan per the user-entered parameters in the xy_scan qlineedits
            self.isStopped = 0
            """
            This cell runs a XZ scan. It currently creates and then populates a user defined size numpy array according to a set counter acquisition time and a motor
step voltage setting. Additionally, the initial driving voltage for the X and Y motors can be set according to the desired scanning range. This scanning program runs in a snake pattern, it scan the first row left to right, moves up one row, 
then scans right to left and continues. Alternatives would be scanning left to right, and resetting the position of the laser on the next higher row and scanning again left 
to right OR scanning in a "circular" patter either CW or CCW from outside to inside or inside to outside. The chosen method was picked for simplicity of understanding. The 
scanning loops are present within NI-DAQmx tasks to speed up the program. Starting and stopping a NI-DAQmx task repeatedly slows down the program dramatically. So, the 
counter and hardware clock task are started once, then the scanning program is run, and then the counter and clock tasks are closed -un-reserving the hardware resources. 
This cell uses the "DAQAnalogOutputs" function from a written class file at:
C:/Users/lukin2dmaterials/miniconda3/envs/qcodes/Lib/site-packages/qcodes_contrib_drivers/drivers/NationalInstruments/class_file. Slashes are reversed to run
            """

            ################################################################################## card 2 (AO) ########################################################################

            try: 
                scan_galvo_card_name = "cDAQ1Mod2" # naming the instrument
                scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)} # dictionary of analog output channels
                scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels) # defining the instrument (ni_9263)

                get_todays_date = date.today() # this is used for creating the final plot's plot labels
                todays_date = get_todays_date.strftime("%m%d%Y") # this is used for creating the final plot's plot labels
                self.scanType = 'xz'

                ############################################################################### def other variables #####################################################################
                ################### setting variales and array ####################

                # acquisition and settling time:
                self.time_acquire = np.round(float(xz_scan_read_time_qlineedit.text())/1e3, 4)
                self.time_settle = np.round(float(xz_scan_settle_time_qlineedit.text())/1e3, 5)
                self.time_pts_settle = int(self.time_settle*self.samp_rate)
                self.time_pts_acquire = int(self.time_acquire*self.samp_rate)
                self.time_pts_per_pixel = self.time_pts_acquire  + self.time_pts_settle
                print('self.time_pts_settle: ' + str(self.time_pts_settle))
                print('self.time_pts_acquire: ' + str(self.time_pts_acquire))
                
                # resolution
                self.hor_res = grid_size_x = int(xz_scan_resolution_hor_qlineedit.text())
                self.ver_res = grid_size_z = int(xz_scan_resolution_ver_qlineedit.text())

                # initial driving voltages for the x,y-mirror and z piezo stage
                self.hor_init = x_init = round(float(xz_scan_x_voltage_min_qlineedit.text()),3)
                self.y = y_init = round(float(xz_scan_y_voltage_qlineedit.text()),3)
                self.ver_init = z_init = round(float(xz_scan_z_voltage_min_qlineedit.text()),3)

                # final driving voltages
                x_final = round(float(xz_scan_x_voltage_max_qlineedit.text()),3)
                z_final = round(float(xz_scan_z_voltage_max_qlineedit.text()),3)

                # range of voltage
                self.hor_range = x_final - x_init
                self.ver_range = z_final - z_init

                # array of x,z voltages
                self.x_array = np.linspace(x_init, x_final, grid_size_x)
                self.z_array = np.linspace(z_init, z_final, grid_size_z)
                X, Z = np.meshgrid(self.x_array, self.z_array)

                self.hor_array = self.x_array
                self.ver_array = self.z_array

                # make the x-voltage waveform to pass to aotask
                self.x_array = np.array(np.repeat(self.x_array, self.time_pts_per_pixel))
                print("len(self.x_array): " + str(len(self.x_array)))

                # create dataset to populate
                global xz_scan_data_array
                xz_scan_data_array = np.zeros((grid_size_z, grid_size_x))
                global most_recent_data_array
                most_recent_data_array = np.zeros((grid_size_z, grid_size_x))

                print('XZ scan finished')
                print('----------------------------------------------------------------')

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
                    if f == 0 or f == grid_size_z-1 or f == int(grid_size_z/2):
                        print(end)

                    ##################### updating plot section ####################
                    self.sc.axes.pcolormesh(X, Z, xz_scan_data_array, cmap = "inferno")
                    self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
                    self.sc.axes.set_ylabel("z_stage_driving_voltage_(V)", fontsize = 8)
                    self.sc.fig.canvas.draw()
                    self.sc.fig.canvas.flush_events()

                clktask.close(); scan_galvo.close() 

                ############################################ end XZ scanning script #################################################
                ############################################### plotting XZ scan data in plot ####################################################
                self.sc.axes.cla()
                self.result = xz_scan_data_array

                plot = self.sc.axes.pcolormesh(X, Z, xz_scan_data_array, cmap = "inferno")
                self.xz_scan_plot_colorbar = self.sc.fig.colorbar(plot, ax = self.sc.axes, pad = 0.02, aspect = 15)
                self.xz_scan_plot_colorbar.formatter.set_powerlimits((0, 0))
                self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                self.xz_scan_plot_colorbar.ax.tick_params(labelsize = 7)
                self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
                self.sc.axes.set_ylabel("z_stage_driving_voltage_(V)", fontsize = 8, labelpad = -10)
                self.sc.axes.set_title("XZ_scan_%s_y_voltage=%f_V" % (todays_date, y_init), fontsize = 8)
                self.crosshair, = self.sc.axes.plot(self.x_array[0], self.z_array[0], marker="+", markersize=5, markeredgecolor="red", markerfacecolor="red")
                self.sc.draw()

                cid = self.sc.mpl_connect('button_press_event', mouse_event)
                most_recent_data_array = xz_scan_data_array

            except nidaqmx.errors.DaqError:
                print('aaaaaa')
                scan_galvo.close()
                aotask.stop(); aotask.close()
                ctrtask.stop(); ctrtask.close()
                clktask.stop(); clktask.close() 

        # xz_scan resolution check then run fnc
        def xz_scan_resolution_validation_fnc():

            self.sc.axes.cla()

            try:
                self.xz_scan_plot_colorbar.remove()

            except (AttributeError, ValueError):
                pass

            try:
                self.yz_scan_plot_colorbar.remove()

            except (AttributeError, ValueError):
                pass
            
            try:
                self.xy_scan_plot_colorbar.remove()
            
            except (AttributeError, ValueError):
                pass

            res_min_condition = 0 # set the min allowed resolution for scanning
            res_max_condition = 900 # set the max allowed resolution for scanning

            xz_scan_resolution_test_condition = False # define resolution validation bool for xz scan

            while xz_scan_resolution_test_condition is False: # this initiates checking the resolution parameter

                # checking for out of bounds of min and max conditions above
                if int(xz_scan_resolution_hor_qlineedit.text()) < res_min_condition or int(xz_scan_resolution_hor_qlineedit.text()) > res_max_condition: # TODO: or negative or not a number or too large

                    display_resolution_error_window_fnc() # call the error message pop-up window
                    break # exit the checking loop: failed

                # if parameter is in bounds; run scan
                elif int(xz_scan_resolution_hor_qlineedit.text()) >= res_min_condition and int(xz_scan_resolution_hor_qlineedit.text()) <= res_max_condition:

                    xz_scan_resolution_test_condition == True
                    print_XZ_scan_parameters_fnc(self) # call the print user-entered parameters fnc
                    run_xz_scan_fnc() # call the run xz scan method fnc
                    break # exit the checking loop: passed

        ########################## YZ scanning ############################
        # print_YZ_scan_parameters_fnc
        def print_YZ_scan_parameters_fnc(self, parent = Child): # this fnc does...

            # need to clear text box first
            self.parameters_display_text_box.clear()

            # this prints to the QTextBox in the left_window. The output of the user-selected scan parameters is printed below
            self.parameters_display_text_box.setPlainText(
                        "YZ_SCAN PARAMETERS/INFO:\n"
                        "YZ_scan resolution = " + str(int(yz_scan_resolution_hor_qlineedit.text())) + "\n"
                        "YZ_scan counter read time = " + str(float(yz_scan_read_time_qlineedit.text())/1e3) + "\n"
                        "YZ_scan min Y driving voltage = " + str(float(yz_scan_y_voltage_min_qlineedit.text())) + "\n"
                        "YZ_scan max Y driving voltage = " + str(float(yz_scan_y_voltage_max_qlineedit.text())) + "\n"
                        "YZ_scan X driving voltage = " + str(float(yz_scan_x_voltage_qlineedit.text())) + "\n'"
                        "YZ_scan z-piezo min driving voltage = " + str(float(yz_scan_z_voltage_min_qlineedit.text())) + "\n"
                        "YZ_scan z-piezo max driving voltage = " + str(float(yz_scan_z_voltage_max_qlineedit.text()))
                                                        )
                
        # run YZ scan function
        def run_yz_scan_fnc(): # this fnc runs the YZ scan script
            self.isStopped = 0
            """
            This cell runs a YZ only scan. It currently creates and then populates a user defined size numpy array according to a set counter acquisition time and a motor
step voltage setting. Additionally, the initial driving voltage for the X and Y motors can be set according to the desired scanning range. This scanning program runs in a snake pattern, it scan the first row left to right, moves up one row, 
then scans right to left and continues. Alternatives would be scanning left to right, and resetting the position of the laser on the next higher row and scanning again left 
to right OR scanning in a "circular" patter either CW or CCW from outside to inside or inside to outside. The chosen method was picked for simplicity of understanding. The 
scanning loops are present within NI-DAQmx tasks to speed up the program. Starting and stopping a NI-DAQmx task repeatedly slows down the program dramatically. So, the 
counter and hardware clock task are started once, then the scanning program is run, and then the counter and clock tasks are closed -un-reserving the hardware resources. 
This cell uses the "DAQAnalogOutputs" function from a written class file at:
C:/Users/lukin2dmaterials/miniconda3/envs/qcodes/Lib/site-packages/qcodes_contrib_drivers/drivers/NationalInstruments/class_file. Slashes are reversed to run
            """

            ################################################################################## card 2 (AO) ########################################################################
            try: 
                scan_galvo_card_name = "cDAQ1Mod2" # naming the instrument
                scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)} # dictionary of analog output channels
                scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels) # defining the instrument (ni_9263)

                get_todays_date = date.today() # this is used for creating the final plot's plot labels
                todays_date = get_todays_date.strftime("%m%d%Y") # this is used for creating the final plot's plot labels
                self.scanType = 'yz'

                ############################################################################### def other variables #####################################################################
                ################### setting variales and array ####################
                
                # acquisition and settling time:
                self.time_acquire = np.round(float(yz_scan_read_time_qlineedit.text())/1e3, 4)
                self.time_settle = np.round(float(yz_scan_settle_time_qlineedit.text())/1e3, 5)
                self.time_pts_settle = int(self.time_settle*self.samp_rate)
                self.time_pts_acquire = int(self.time_acquire*self.samp_rate)
                self.time_pts_per_pixel = self.time_pts_acquire  + self.time_pts_settle
                print('self.time_pts_settle: ' + str(self.time_pts_settle))
                print('self.time_pts_acquire: ' + str(self.time_pts_acquire))

                # resolution
                self.hor_res = grid_size_y = int(yz_scan_resolution_hor_qlineedit.text())
                self.ver_res = grid_size_z = int(yz_scan_resolution_ver_qlineedit.text())

                # initial driving voltages for the x,y-mirror and z piezo stage
                self.x = x_init = round(float(yz_scan_x_voltage_qlineedit.text()), 3)
                self.hor_init = y_init = round(float(yz_scan_y_voltage_min_qlineedit.text()),3)
                self.ver_init = z_init = round(float(yz_scan_z_voltage_min_qlineedit.text()),3)

                # final driving voltages
                y_final = round(float(yz_scan_y_voltage_max_qlineedit.text()),3)
                z_final = round(float(yz_scan_z_voltage_max_qlineedit.text()),3)

                # range of voltage
                self.hor_range = y_final - y_init
                self.ver_range = z_final - z_init

                # array of y, z voltages
                self.y_array = np.linspace(y_init, y_final, grid_size_y)
                self.z_array = np.linspace(z_init, z_final, grid_size_z)
                Y, Z = np.meshgrid(self.y_array, self.z_array)

                self.hor_array = self.y_array
                self.ver_array = self.z_array

                # make the y-voltage waveform to pass to aotask
                self.y_array = np.array(np.repeat(self.y_array, self.time_pts_per_pixel))
                print("len(self.y_array): " + str(len(self.y_array)))
            
                # create dataset to populate
                global yz_scan_data_array
                yz_scan_data_array = np.zeros((grid_size_z, grid_size_y))
                global most_recent_data_array
                most_recent_data_array = np.zeros((grid_size_z, grid_size_y))

                print('YZ scan finished')
                print('----------------------------------------------------------------')

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
                    samps_per_chan = int(len(self.y_array) + 1) # does this matter?
                    )

                ######################################################################## Y and Z scanning #########################################################################
                print('----------------------------------------------------------------')
                print('YZ scan started')
                for f in range(grid_size_z): # this loops for each row z
                    if self.isStopped == 1: break
                    start0 = time.time()

                    # set initial z location
                    scan_galvo.voltage_cdaq1mod2ao2(self.z_array[f])
                    time.sleep(100*self.time_settle)
                    
                    # snake pattern: scanning right to left
                    if self.isSnake == 1 and f % 2 == 1: self.y_array_write = np.ascontiguousarray(self.y_array[::-1])
                    else: self.y_array_write = self.y_array

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
                        samps_per_chan = int(len(self.y_array) + 1) 
                        )

                    ##############################################################
                    # AO task for the y-galvo
                    aotask = nidaqmx.Task() 
                    channel_galvo_y = f'{scan_galvo_card_name}/ao{1}'
                    aotask.ao_channels.add_ao_voltage_chan(channel_galvo_y)
                    aotask.timing.cfg_samp_clk_timing(
                        rate = self.samp_rate,
                        source = "/cDAQ1/Ctr1InternalOutput", # the clock defined above
                        active_edge = nidaqmx.constants.Edge.RISING,
                        sample_mode = AcquisitionType.FINITE,
                        samps_per_chan = len(self.y_array)
                        )
                    aotask.write(self.y_array_write, auto_start=False) # assign the x-waveform to aotask, but not start yet

                    # start galvo Y scan, counter, and clock
                    aotask.start()
                    ctrtask.start()
                    clktask.start() # the clock MUST START AFTER aotask.start() and ctrtask.start()
                    
                    aotask.wait_until_done()
                    yLineData = ctrtask.read(len(self.y_array) + 1)  # +1 is to take the difference later
                    # print("len(yLineData): " + str(len(yLineData)))
                    
                    aotask.stop(); aotask.close()
                    ctrtask.stop(); ctrtask.close()
                    clktask.stop() 
                    diffData = np.diff(yLineData)

                    summedData = np.zeros(grid_size_y)
                    for i in range(0, grid_size_y):
                        summedData[i] = np.sum(
                            diffData[(i*self.time_pts_per_pixel + self.time_pts_settle):(i*self.time_pts_per_pixel + self.time_pts_per_pixel - 1)])
                    normedData = summedData * (.001 / self.time_acquire) # normalize to kcounts/sec

                    if self.isSnake == 1 and f % 2 == 1: yz_scan_data_array[f] = np.flip(normedData)
                    else: yz_scan_data_array[f] = normedData

                    end = time.time() - start0
                    if f == 0 or f == grid_size_z-1 or f == int(grid_size_z/2):
                        print(end)
                    
                    self.sc.axes.pcolormesh(Y, Z, yz_scan_data_array, cmap = "inferno")
                    self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.set_xlabel("y_mirror_driving_voltage_(V)", fontsize = 8)
                    self.sc.axes.set_ylabel("z_stage_driving_voltage_(V)", fontsize = 8, labelpad = -9)
                    self.sc.fig.canvas.draw()
                    self.sc.fig.canvas.flush_events() # this line is very important
                        
                clktask.close(); scan_galvo.close()

                ############################################### plotting YZ scan data in plot ####################################################
                self.sc.axes.cla()
                self.result = yz_scan_data_array

                yz_scan_plot = self.sc.axes.pcolormesh(Y, Z, yz_scan_data_array, cmap = "inferno")
                self.yz_scan_plot_colorbar = self.sc.fig.colorbar(yz_scan_plot, ax = self.sc.axes, pad = 0.02, aspect = 15)
                self.yz_scan_plot_colorbar.formatter.set_powerlimits((0, 0))
                self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                self.yz_scan_plot_colorbar.ax.tick_params(labelsize = 7)
                self.sc.axes.set_xlabel("y_mirror_driving_voltage_(V)", fontsize = 8)
                self.sc.axes.set_ylabel("z_stage_driving_voltage_(V)", fontsize = 8, labelpad=-10)
                self.sc.axes.set_title("YZ_scan_%s_x_voltage=%f_V" % (todays_date, x_init), fontsize = 8)
                self.crosshair, = self.sc.axes.plot(self.y_array[0], self.z_array[0], marker="+", markersize=5, markeredgecolor="red", markerfacecolor="red")
                self.sc.draw()

                cid = self.sc.mpl_connect('button_press_event', mouse_event)
                most_recent_data_array = yz_scan_data_array

            except nidaqmx.errors.DaqError:
                print('aaaaaa')
                scan_galvo.close()
                aotask.stop(); aotask.close()
                ctrtask.stop(); ctrtask.close()
                clktask.stop(); clktask.close() 
    
        # yz_scan resolution check then run fnc
        def yz_scan_resolution_validation_fnc():

            self.sc.axes.cla()

            try:
                self.yz_scan_plot_colorbar.remove()

            except (AttributeError, ValueError):
                pass
            
            try:
                self.xy_scan_plot_colorbar.remove()
            
            except (AttributeError, ValueError):
                pass
                
            try:
                self.xz_scan_plot_colorbar.remove()
            
            except (AttributeError, ValueError):
                pass

            res_min_condition = 0 # set the min allowed resolution for scanning
            res_max_condition = 900 # set the max allowed resolution for scanning

            yz_scan_resolution_test_condition = False # define resolution validation bool for yz scan

            while yz_scan_resolution_test_condition is False: # this initiates checking the resolution parameter

                # checking for out of bounds of min and max conditions above
                if int(yz_scan_resolution_hor_qlineedit.text()) < res_min_condition or int(yz_scan_resolution_hor_qlineedit.text()) > res_max_condition: # TODO: or negative or not a number or too large

                    display_resolution_error_window_fnc() # call the error message pop-up window
                    break # exit the checking loop: failed

                # if parameter is in bounds; run scan
                elif int(yz_scan_resolution_hor_qlineedit.text()) >= res_min_condition and int(yz_scan_resolution_hor_qlineedit.text()) <= res_max_condition:

                    yz_scan_resolution_test_condition == True
                    print_YZ_scan_parameters_fnc(self) # call the print user-entered parameters fnc
                    run_yz_scan_fnc() # call the run yz scan method fnc
                    break # exit the checking loop: passed


        ############################################################# left half window #####################################################################
        left_window = QFrame(self)
        left_window.setFrameShape(QFrame.StyledPanel)
        left_window.setFixedSize(340, 430)

        # QTextEdit display for printing scan parameters to GUI window
        """
        This creates a textbox located in the lower left hand corner of the GUI. It is set (parent) to the left_window. When any scan script is run the 
        parameters enterd by the user are collected from their respective QLineedits and then printed in a list to this QTextBox. This allows the user 
        select and copy the parameters they used quickly in order to save them to an external file as documentation of scans run.

        Room for improvement:
        1. The date (get from "todays_date" that is already implemented) could be added to the first line of any scans parameter printing fnc
        2. The info printed to this QTextBox could also and/or be saved to an external file for record/documentation
        """

        self.parameters_display_text_box = QTextEdit(left_window)
        self.parameters_display_text_box.resize(320, 100)
        self.parameters_display_text_box.move(10, 322)
        self.parameters_display_text_box.setParent(left_window)

        ############################################################### right half window ###################################################################
        right_window = QFrame(self)
        right_window.setFrameShape(QFrame.StyledPanel)
        right_window.setFixedSize(430, 430)

        ##################################################### plot in right window ##################################################################

        plot_res = 3.89
        self.sc = MplCanvas(self, width = plot_res, height = plot_res, dpi = 110)                         # changing dpi does something to scale of figure
        self.sc.move(2, 2)
        self.sc.setParent(right_window)
        self.sc.axes.xaxis.set_tick_params(labelsize = 8)
        self.sc.axes.yaxis.set_tick_params(labelsize = 8)
        # self.sc.fig.colorbar(plot, ax = self.sc.axes)
        # self.sc.axes.set_title("plot_title_here", fontsize = 9)
        # self.sc.axes.set_xlabel("plot_x_label_here", fontsize = 8)
        # self.sc.axes.set_ylabel("plot_y_label_here", fontsize = 8)

        ############################################################# split left and right windows #########################################################
        
        splitter1 = QSplitter(Qt.Horizontal)
        splitter1.addWidget(left_window)
        splitter1.addWidget(right_window)
        hbox.addWidget(splitter1) # set layout and show window
        self.setLayout(hbox)
        self.show()

        ############################################################# scanning (XY, XZ, & YZ) section ############################################################################
        
        ##################################### overall ####################################

        ############ begin save data section ###############
        save_scan_data_button = QPushButton("Save scan data:", self) # create the save scan data button
        save_scan_data_button.setParent(left_window) # set the "parent" bound of the save scan data button
        save_scan_data_button.resize(90, 20) # resize the save scan data button
        save_scan_data_button.move(10, 280) # position the save scan data button in the left winodw
        save_scan_data_button.clicked.connect(save_scan_data_fnc)

        save_scan_data_qlineedit = QLineEdit(self) # qlineedit
        save_scan_data_qlineedit.setParent(left_window)
        save_scan_data_qlineedit.resize(189, 20)
        save_scan_data_qlineedit.move(105, 280)

        save_scan_data_extension_widget = QLabel("\".npy\"", self) # widget
        save_scan_data_extension_widget.setParent(left_window),
        save_scan_data_extension_widget.move(297, 288 - 5)
        ############ ending save data section ###############

        # scan widgets overall "Scanning"
        scan_widget_overall = QLabel("Scanning options:")
        # scan_widget_overall = QLabel()
        # scan_widget_overall = QWebEngineView()
        # scan_widget_overall.setHtml(pageSource)
        scan_widget_overall.setFont(QFont("Times font", 9))
        scan_widget_overall.setParent(left_window)
        scan_widget_overall.move(120, 10)
        # scan_widget_overall.move(30, 10)

        indiv_scan_labels_y_height = 17

        row_y_adjust = 5

        overall_y_adjust = 10

        #########################################################################################
        ########################################### XY scanning #################################
        #########################################################################################

        xy_scan_widgets_left_x_justify = 5

        # XY scan
        xy_scan_label_widget = QLabel("XY scan", self) # widget
        xy_scan_label_widget.setParent(left_window)
        xy_scan_label_widget.move(12 + 10, indiv_scan_labels_y_height + overall_y_adjust)

        # resolution in x
        xy_scan_resolution_hor_widget = QLabel("ResX:", self) # widget
        xy_scan_resolution_hor_widget.setParent(left_window)
        xy_scan_resolution_hor_widget.move(xy_scan_widgets_left_x_justify, 30 + row_y_adjust + overall_y_adjust)

        xy_scan_resolution_hor_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_resolution_hor_qlineedit.setParent(left_window)
        xy_scan_resolution_hor_qlineedit.move(30, 30 + row_y_adjust + overall_y_adjust)
        xy_scan_resolution_hor_qlineedit.resize(55, 15)
        xy_scan_resolution_hor_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # resolution in y
        xy_scan_resolution_ver_widget = QLabel("ResY:", self) # widget
        xy_scan_resolution_ver_widget.setParent(left_window)
        xy_scan_resolution_ver_widget.move(xy_scan_widgets_left_x_justify, 45 + row_y_adjust + overall_y_adjust)

        xy_scan_resolution_ver_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_resolution_ver_qlineedit.setParent(left_window)
        xy_scan_resolution_ver_qlineedit.move(30, 45 + row_y_adjust + overall_y_adjust)
        xy_scan_resolution_ver_qlineedit.resize(55, 15)
        xy_scan_resolution_ver_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # read time
        xy_scan_read_time_widget = QLabel("APD_t (ms)", self) # widget
        xy_scan_read_time_widget.setParent(left_window)
        xy_scan_read_time_widget.move(xy_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)

        xy_scan_read_time_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_read_time_qlineedit.setParent(left_window)
        xy_scan_read_time_qlineedit.move(60, 65 + row_y_adjust + overall_y_adjust)
        xy_scan_read_time_qlineedit.resize(25, 15)
        xy_scan_read_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # settle time
        xy_scan_settle_time_widget = QLabel("settle (ms)", self) # widget
        xy_scan_settle_time_widget.setParent(left_window)
        xy_scan_settle_time_widget.move(xy_scan_widgets_left_x_justify, 80 + row_y_adjust + overall_y_adjust)

        xy_scan_settle_time_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_settle_time_qlineedit.setParent(left_window)
        xy_scan_settle_time_qlineedit.move(60, 80 + row_y_adjust + overall_y_adjust)
        xy_scan_settle_time_qlineedit.resize(25, 15)
        xy_scan_settle_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
        
        # x voltage (min and max)
        xy_scan_x_voltage_min_widget = QLabel("x_V_min:", self) # widget
        xy_scan_x_voltage_min_widget.setParent(left_window)
        xy_scan_x_voltage_min_widget.move(xy_scan_widgets_left_x_justify, 95 + row_y_adjust + overall_y_adjust)

        xy_scan_x_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_x_voltage_min_qlineedit.setParent(left_window)
        xy_scan_x_voltage_min_qlineedit.move(50, 95 + row_y_adjust + overall_y_adjust)
        xy_scan_x_voltage_min_qlineedit.resize(35, 15)
        xy_scan_x_voltage_min_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        xy_scan_x_voltage_max_widget = QLabel("x_V_max:", self) # widget
        xy_scan_x_voltage_max_widget.setParent(left_window)
        xy_scan_x_voltage_max_widget.move(xy_scan_widgets_left_x_justify, 115 + row_y_adjust + overall_y_adjust)

        xy_scan_x_voltage_max_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_x_voltage_max_qlineedit.setParent(left_window)
        xy_scan_x_voltage_max_qlineedit.move(55, 115 + row_y_adjust + overall_y_adjust)
        xy_scan_x_voltage_max_qlineedit.resize(30, 15)
        xy_scan_x_voltage_max_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # y voltage (min and max)
        xy_scan_y_voltage_min_widget = QLabel("y_V_min:", self) # widget
        xy_scan_y_voltage_min_widget.setParent(left_window)
        xy_scan_y_voltage_min_widget.move(xy_scan_widgets_left_x_justify, 140 + row_y_adjust + overall_y_adjust)

        xy_scan_y_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_y_voltage_min_qlineedit.setParent(left_window)
        xy_scan_y_voltage_min_qlineedit.move(50, 140 + row_y_adjust + overall_y_adjust)
        xy_scan_y_voltage_min_qlineedit.resize(35, 15)
        xy_scan_y_voltage_min_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        xy_scan_y_voltage_max_widget = QLabel("y_V_max:", self) # widget
        xy_scan_y_voltage_max_widget.setParent(left_window)
        xy_scan_y_voltage_max_widget.move(xy_scan_widgets_left_x_justify, 165 + row_y_adjust + overall_y_adjust)

        xy_scan_y_voltage_max_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_y_voltage_max_qlineedit.setParent(left_window)
        xy_scan_y_voltage_max_qlineedit.move(55, 165 + row_y_adjust + overall_y_adjust)
        xy_scan_y_voltage_max_qlineedit.resize(30, 15)
        xy_scan_y_voltage_max_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
        
        # z piezo
        xy_scan_z_piezo_voltage_widget = QLabel("z_V:", self) # widget
        xy_scan_z_piezo_voltage_widget.setParent(left_window)
        xy_scan_z_piezo_voltage_widget.move(xy_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)

        xy_scan_z_voltage_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_z_voltage_qlineedit.setParent(left_window)
        xy_scan_z_voltage_qlineedit.move(30, 190 + row_y_adjust + overall_y_adjust)
        xy_scan_z_voltage_qlineedit.resize(55, 15)
        xy_scan_z_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # run xy scan button
        xy_scan_run_button = QPushButton("run\nXY scan", self) # button
        xy_scan_run_button.setParent(left_window)
        xy_scan_run_button.resize(60, 40)
        xy_scan_run_button.move(10, 215 + row_y_adjust + overall_y_adjust)
        xy_scan_run_button.clicked.connect(xy_scan_resolution_validation_fnc) # this framework is limited currently to only validating resolution


        #######################################################################################
        ####################################### XZ scanning ###################################
        #######################################################################################

        xz_scan_widgets_left_x_justify = 130

        # scan widget 2 "XZ"
        scan_widget_2 = QLabel("XZ scan", self)
        scan_widget_2.setParent(left_window)
        scan_widget_2.move(145 + 5, indiv_scan_labels_y_height + overall_y_adjust)

        # resolution in x
        xz_scan_resolution_hor_widget = QLabel("ResX:", self) # widget
        xz_scan_resolution_hor_widget.setParent(left_window)
        xz_scan_resolution_hor_widget.move(xz_scan_widgets_left_x_justify, 30 + row_y_adjust + overall_y_adjust)

        xz_scan_resolution_hor_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_resolution_hor_qlineedit.setParent(left_window)
        xz_scan_resolution_hor_qlineedit.move(25 + xz_scan_widgets_left_x_justify, 30 + row_y_adjust + overall_y_adjust)
        xz_scan_resolution_hor_qlineedit.resize(55, 15)
        xz_scan_resolution_hor_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # resolution in z
        xz_scan_resolution_ver_widget = QLabel("ResZ:", self) # widget
        xz_scan_resolution_ver_widget.setParent(left_window)
        xz_scan_resolution_ver_widget.move(xz_scan_widgets_left_x_justify, 45 + row_y_adjust + overall_y_adjust)

        xz_scan_resolution_ver_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_resolution_ver_qlineedit.setParent(left_window)
        xz_scan_resolution_ver_qlineedit.move(25 + xz_scan_widgets_left_x_justify, 45 + row_y_adjust + overall_y_adjust)
        xz_scan_resolution_ver_qlineedit.resize(55, 15)
        xz_scan_resolution_ver_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # read time
        xz_scan_read_time_widget = QLabel("APD_t (ms):", self) # widget
        xz_scan_read_time_widget.setParent(left_window)
        xz_scan_read_time_widget.move(xz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)

        xz_scan_read_time_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_read_time_qlineedit.setParent(left_window)
        xz_scan_read_time_qlineedit.move(55 + xz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)
        xz_scan_read_time_qlineedit.resize(25, 15)
        xz_scan_read_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # settle time
        xz_scan_settle_time_widget = QLabel("settle (ms)", self) # widget
        xz_scan_settle_time_widget.setParent(left_window)
        xz_scan_settle_time_widget.move(xz_scan_widgets_left_x_justify, 80 + row_y_adjust + overall_y_adjust)

        xz_scan_settle_time_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_settle_time_qlineedit.setParent(left_window)
        xz_scan_settle_time_qlineedit.move(55 + xz_scan_widgets_left_x_justify, 80 + row_y_adjust + overall_y_adjust)
        xz_scan_settle_time_qlineedit.resize(25, 15)
        xz_scan_settle_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
        
        # x voltage (min and max)
        xz_scan_x_voltage_min_widget = QLabel("x_V_min:", self) # widget
        xz_scan_x_voltage_min_widget.setParent(left_window)
        xz_scan_x_voltage_min_widget.move(xz_scan_widgets_left_x_justify, 95 + row_y_adjust + overall_y_adjust)

        xz_scan_x_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_x_voltage_min_qlineedit.setParent(left_window)
        xz_scan_x_voltage_min_qlineedit.move(45 + xz_scan_widgets_left_x_justify, 95 + row_y_adjust + overall_y_adjust)
        xz_scan_x_voltage_min_qlineedit.resize(35, 15)
        xz_scan_x_voltage_min_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        xz_scan_x_voltage_max_widget = QLabel("x_V_max:", self) # widget
        xz_scan_x_voltage_max_widget.setParent(left_window)
        xz_scan_x_voltage_max_widget.move(xz_scan_widgets_left_x_justify, 115 + row_y_adjust + overall_y_adjust)

        xz_scan_x_voltage_max_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_x_voltage_max_qlineedit.setParent(left_window)
        xz_scan_x_voltage_max_qlineedit.move(50 + xz_scan_widgets_left_x_justify, 115 + row_y_adjust + overall_y_adjust)
        xz_scan_x_voltage_max_qlineedit.resize(30, 15)
        xz_scan_x_voltage_max_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # y voltage single setting
        xz_scan_y_voltage_widget = QLabel("y_V:", self) # widget
        xz_scan_y_voltage_widget.setParent(left_window)
        xz_scan_y_voltage_widget.move(xz_scan_widgets_left_x_justify, 140 + row_y_adjust + overall_y_adjust)

        xz_scan_y_voltage_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_y_voltage_qlineedit.setParent(left_window)
        xz_scan_y_voltage_qlineedit.move(25 + xz_scan_widgets_left_x_justify, 140 + row_y_adjust + overall_y_adjust)
        xz_scan_y_voltage_qlineedit.resize(55, 15)
        xz_scan_y_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
                                
        # z piezo (min and max)
        xz_scan_z_piezo_min_voltage_widget = QLabel("z_V_min:", self) # widget
        xz_scan_z_piezo_min_voltage_widget.setParent(left_window)
        xz_scan_z_piezo_min_voltage_widget.move(xz_scan_widgets_left_x_justify, 165 + row_y_adjust + overall_y_adjust)

        xz_scan_z_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_z_voltage_min_qlineedit.setParent(left_window)
        xz_scan_z_voltage_min_qlineedit.move(45 + xz_scan_widgets_left_x_justify, 165 + row_y_adjust + overall_y_adjust)
        xz_scan_z_voltage_min_qlineedit.resize(35, 15)
        xz_scan_z_voltage_min_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        xz_scan_z_piezo_max_voltage_widget = QLabel("z_V_max:", self) # widget
        xz_scan_z_piezo_max_voltage_widget.setParent(left_window)
        xz_scan_z_piezo_max_voltage_widget.move(xz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)

        xz_scan_z_voltage_max_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_z_voltage_max_qlineedit.setParent(left_window)
        xz_scan_z_voltage_max_qlineedit.move(50 + xz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)
        xz_scan_z_voltage_max_qlineedit.resize(30, 15)
        xz_scan_z_voltage_max_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # run xz scan button
        xz_scan_run_button = QPushButton("run\nXZ scan", self) # button
        xz_scan_run_button.setParent(left_window)
        xz_scan_run_button.resize(60, 40)
        xz_scan_run_button.move(10 + xz_scan_widgets_left_x_justify, 215 + row_y_adjust + overall_y_adjust)
        xz_scan_run_button.clicked.connect(xz_scan_resolution_validation_fnc) # this framework is limited currently to only validating resolution

        #####################################################################################
        ####################################### YZ scanning #################################
        #####################################################################################

        yz_scan_widgets_left_x_justify = 255

        # scan widget 3 "YZ"
        scan_widget_3 = QLabel("YZ scan", self)
        scan_widget_3.setParent(left_window)
        scan_widget_3.move(265 + 14, indiv_scan_labels_y_height + overall_y_adjust)

        # resolution in y
        yz_scan_resolution_hor_widget = QLabel("ResY:", self) # widget
        yz_scan_resolution_hor_widget.setParent(left_window)
        yz_scan_resolution_hor_widget.move(yz_scan_widgets_left_x_justify, 30 + row_y_adjust + overall_y_adjust)

        yz_scan_resolution_hor_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_resolution_hor_qlineedit.setParent(left_window)
        yz_scan_resolution_hor_qlineedit.move(25 + yz_scan_widgets_left_x_justify, 30 + row_y_adjust + overall_y_adjust)
        yz_scan_resolution_hor_qlineedit.resize(55, 15)
        yz_scan_resolution_hor_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # resolution in z
        yz_scan_resolution_ver_widget = QLabel("ResZ:", self) # widget
        yz_scan_resolution_ver_widget.setParent(left_window)
        yz_scan_resolution_ver_widget.move(yz_scan_widgets_left_x_justify, 45 + row_y_adjust + overall_y_adjust)

        yz_scan_resolution_ver_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_resolution_ver_qlineedit.setParent(left_window)
        yz_scan_resolution_ver_qlineedit.move(25 + yz_scan_widgets_left_x_justify, 45 + row_y_adjust + overall_y_adjust)
        yz_scan_resolution_ver_qlineedit.resize(55, 15)
        yz_scan_resolution_ver_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # read time
        yz_scan_read_time_widget = QLabel("APD_t (ms):", self) # widget
        yz_scan_read_time_widget.setParent(left_window)
        yz_scan_read_time_widget.move(yz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)

        yz_scan_read_time_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_read_time_qlineedit.setParent(left_window)
        yz_scan_read_time_qlineedit.move(55 + yz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)
        yz_scan_read_time_qlineedit.resize(25, 15)
        yz_scan_read_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # settle time
        yz_scan_settle_time_widget = QLabel("settle (ms)", self) # widget
        yz_scan_settle_time_widget.setParent(left_window)
        yz_scan_settle_time_widget.move(yz_scan_widgets_left_x_justify, 80 + row_y_adjust + overall_y_adjust)

        yz_scan_settle_time_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_settle_time_qlineedit.setParent(left_window)
        yz_scan_settle_time_qlineedit.move(55 + yz_scan_widgets_left_x_justify, 80 + row_y_adjust + overall_y_adjust)
        yz_scan_settle_time_qlineedit.resize(25, 15)
        yz_scan_settle_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
        
        # y voltage (min and max)
        yz_scan_y_voltage_min_widget = QLabel("y_V_min:", self) # widget
        yz_scan_y_voltage_min_widget.setParent(left_window)
        yz_scan_y_voltage_min_widget.move(yz_scan_widgets_left_x_justify, 95 + row_y_adjust + overall_y_adjust)

        yz_scan_y_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_y_voltage_min_qlineedit.setParent(left_window)
        yz_scan_y_voltage_min_qlineedit.move(45 + yz_scan_widgets_left_x_justify, 95 + row_y_adjust + overall_y_adjust)
        yz_scan_y_voltage_min_qlineedit.resize(35, 15)
        yz_scan_y_voltage_min_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        yz_scan_y_voltage_max_widget = QLabel("y_V_max:", self) # widget
        yz_scan_y_voltage_max_widget.setParent(left_window)
        yz_scan_y_voltage_max_widget.move(yz_scan_widgets_left_x_justify, 115 + row_y_adjust + overall_y_adjust)

        yz_scan_y_voltage_max_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_y_voltage_max_qlineedit.setParent(left_window)
        yz_scan_y_voltage_max_qlineedit.move(50 + yz_scan_widgets_left_x_justify, 115 + row_y_adjust + overall_y_adjust)
        yz_scan_y_voltage_max_qlineedit.resize(30, 15)
        yz_scan_y_voltage_max_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # x voltage single setting
        yz_scan_x_voltage_widget = QLabel("x_V:", self) # widget
        yz_scan_x_voltage_widget.setParent(left_window)
        yz_scan_x_voltage_widget.move(yz_scan_widgets_left_x_justify, 140 + row_y_adjust + overall_y_adjust)

        yz_scan_x_voltage_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_x_voltage_qlineedit.setParent(left_window)
        yz_scan_x_voltage_qlineedit.move(25 + yz_scan_widgets_left_x_justify, 140 + row_y_adjust + overall_y_adjust)
        yz_scan_x_voltage_qlineedit.resize(55, 15)
        yz_scan_x_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
                                
        # z piezo (min and max)
        yz_scan_z_piezo_min_voltage_widget = QLabel("z_V_min:", self) # widget
        yz_scan_z_piezo_min_voltage_widget.setParent(left_window)
        yz_scan_z_piezo_min_voltage_widget.move(yz_scan_widgets_left_x_justify, 165 + row_y_adjust + overall_y_adjust)

        yz_scan_z_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_z_voltage_min_qlineedit.setParent(left_window)
        yz_scan_z_voltage_min_qlineedit.move(45 + yz_scan_widgets_left_x_justify, 165 + row_y_adjust + overall_y_adjust)
        yz_scan_z_voltage_min_qlineedit.resize(35, 15)
        yz_scan_z_voltage_min_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        yz_scan_z_piezo_max_voltage_widget = QLabel("z_V_max:", self) # widget
        yz_scan_z_piezo_max_voltage_widget.setParent(left_window)
        yz_scan_z_piezo_max_voltage_widget.move(yz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)

        yz_scan_z_voltage_max_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_z_voltage_max_qlineedit.setParent(left_window)
        yz_scan_z_voltage_max_qlineedit.move(50 + yz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)
        yz_scan_z_voltage_max_qlineedit.resize(30, 15)
        yz_scan_z_voltage_max_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # run yz scan button
        yz_scan_run_button = QPushButton("run\nYZ scan", self) # button
        yz_scan_run_button.setParent(left_window)
        yz_scan_run_button.resize(60, 40)
        yz_scan_run_button.move(15 + yz_scan_widgets_left_x_justify, 215 + row_y_adjust + overall_y_adjust)
        yz_scan_run_button.clicked.connect(yz_scan_resolution_validation_fnc) # this framework is limited currently to only validating resolution

#########################################################################################################
        x_qlineedit = QLineEdit(self) # qclineedit
        x_qlineedit.setParent(left_window)
        x_qlineedit.move(105, 300)
        x_qlineedit.resize(40, 20)
        x_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        y_qlineedit = QLineEdit(self) # qclineedit
        y_qlineedit.setParent(left_window)
        y_qlineedit.move(145, 300)
        y_qlineedit.resize(40, 20)
        y_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        z_qlineedit = QLineEdit(self) # qclineedit
        z_qlineedit.setParent(left_window)
        z_qlineedit.move(185, 300)
        z_qlineedit.resize(40, 20)
        z_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # run set coord button
        set_coord_button = QPushButton("Set coordinate", self) # button
        set_coord_button.setParent(left_window)
        set_coord_button.resize(90, 20)
        set_coord_button.move(10, 300)
        set_coord_button.clicked.connect(set_coordinate_fnc)

        coord_widget = QLabel("x-y-z", self) # widget
        coord_widget.setParent(left_window),
        coord_widget.move(230, 302)

        
        # stop button
        stop_button = QPushButton("Stop scan", self) # button
        stop_button.setParent(left_window)
        stop_button.resize(65,20)
        stop_button.move(265, 300)
        stop_button.clicked.connect(stop_fnc) 

        # snake checkbox
        snake_checkbox = QCheckBox("Snake", self) # button
        snake_checkbox.setParent(left_window)
        snake_checkbox.resize(60,20)
        snake_checkbox.move(210, 25)
        snake_checkbox.clicked.connect(snake_fnc) 

        # snake_shift
        snake_shift_qlineedit = QLineEdit(self) # qclineedit
        snake_shift_qlineedit.setParent(left_window)
        snake_shift_qlineedit.move(210, 45)
        snake_shift_qlineedit.resize(40, 15)
        snake_shift_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # toggle_laser checkbox
        toggle_laser_checkbox_532 = QCheckBox("532", self) # button
        toggle_laser_checkbox_532.setParent(left_window)
        toggle_laser_checkbox_532.resize(60,20)
        toggle_laser_checkbox_532.move(75, 220)
        toggle_laser_checkbox_532.clicked.connect(toggle_laser_fnc_532)

        # toggle_laser checkbox
        toggle_laser_checkbox_589 = QCheckBox("", self) # button
        toggle_laser_checkbox_589.setParent(left_window)
        toggle_laser_checkbox_589.resize(60,20)
        toggle_laser_checkbox_589.move(75-7, 233)
        toggle_laser_checkbox_589.clicked.connect(toggle_laser_fnc_589)
        # toggle_laser checkbox
        toggle_laser_checkbox_589w = QCheckBox("589 S/W", self) # button
        toggle_laser_checkbox_589w.setParent(left_window)
        toggle_laser_checkbox_589w.resize(60,20)
        toggle_laser_checkbox_589w.move(88-7, 233)
        toggle_laser_checkbox_589w.clicked.connect(toggle_laser_fnc_589w)

        # toggle_laser checkbox
        toggle_laser_checkbox_vel = QCheckBox("Velocity", self) # button
        toggle_laser_checkbox_vel.setParent(left_window)
        toggle_laser_checkbox_vel.resize(60,20)
        toggle_laser_checkbox_vel.move(75, 246)
        toggle_laser_checkbox_vel.clicked.connect(toggle_laser_fnc_vel)

        # toggle_laser checkbox
        toggle_laser_checkbox_pt637 = QCheckBox("PT 637", self) # button
        toggle_laser_checkbox_pt637.setParent(left_window)
        toggle_laser_checkbox_pt637.resize(60,20)
        toggle_laser_checkbox_pt637.move(75, 259)
        toggle_laser_checkbox_pt637.clicked.connect(toggle_laser_fnc_pt637)

        # toggle_laser checkbox
        toggle_laser_checkbox_532b = QCheckBox("532b", self) # button
        toggle_laser_checkbox_532b.setParent(left_window)
        toggle_laser_checkbox_532b.resize(60,20)
        toggle_laser_checkbox_532b.move(210, 220)
        toggle_laser_checkbox_532b.clicked.connect(toggle_laser_fnc_532b)

        # toggle_laser checkbox
        toggle_laser_checkbox_pt637b = QCheckBox("PT 637b", self) # button
        toggle_laser_checkbox_pt637b.setParent(left_window)
        toggle_laser_checkbox_pt637b.resize(60,20)
        toggle_laser_checkbox_pt637b.move(210, 233)
        toggle_laser_checkbox_pt637b.clicked.connect(toggle_laser_fnc_pt637b)

        # toggle_laser checkbox
        toggle_laser_checkbox_velb = QCheckBox("Vel 2", self) # button
        toggle_laser_checkbox_velb.setParent(left_window)
        toggle_laser_checkbox_velb.resize(60,20)
        toggle_laser_checkbox_velb.move(210, 246)
        toggle_laser_checkbox_velb.clicked.connect(toggle_laser_fnc_velb)
####################################################################### context menu ######################################################################

    def contextMenuEvent(self, event): # context (right-click) menu

        cmenu = QMenu(self)
        cmenuoneAct = cmenu.addAction("one")
        cmenutwoAct = cmenu.addAction("two")
        cmenuthreeAct = cmenu.addAction("three")
        action = cmenu.exec_(self.mapToGlobal(event.pos()))

############################################################## start gui ################################################################################

def runGUI():
    app = QtWidgets.QApplication(sys.argv)
    mw = Parent()
    mw.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    # t1 = threading.Thread(target=runGUI, args=())
    # t2 = threading.Thread(target=TurnOnLaser.turnOnLaser, args=(2e9,50))
 
    # # starting thread 1
    # t1.start()
    # # starting thread 2
    # t2.start()
    
    # while True:
    #     if t1.is_alive == False: 
    #         t2.join(timeout=1e-3)
    #         break

    # # wait until thread 1 is completely executed
    # t1.join()
    # # wait until thread 2 is completely executed
    # t2.join()

    

    runGUI()

#################################################################### END #######################################################################################
