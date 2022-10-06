"""
B00 gui for room-temp NV center measurement experiment control

Dates: 070622 -> 090222

Author: Miles D. Ackerman (undergraduate during the summer of 2022). Email: miles.ackerman1@gmail.com

General info:
This file makes a GUI to control the room-temp NV center-based measurement setup in LISE room B00 in the Lukin Group in the Harvard PHY dept. The setup
is a confocal laser scanning microscope. A galvo steers a laser beam at 518 nm (butteryfly) through a 4F telescope setup. An excitation and a collection
path (both use fiber coupling) are present. The sample stage is manually-controlled by positioning dials. Additionally, the z-axis of the microscope
objective is controlled by a piezo. This entire setup is controlled by a NI-cDAQ (model 9147). The NI-DAQmx is the National Insturments API that is 
used to communicate with the cDAQ device. Info. about this device can be found at:
https://nidaqmx-python.readthedocs.io/en/latest/https://nidaqmx-python.readthedocs.io/en/latest/. Two cards/modules are present in the cDAQ chassis: 1.
NI-9402 (four bi-directional digital input/output (I/O) channels -also featuring programmable internal clocks for counting). Info about the NI-9402 
module can be found at: https://www.ni.com/en-us/support/model.ni-9402.htmlhttps://www.ni.com/en-us/support/model.ni-9402.html. This module functions 
as a counter based on an internal clock within the module. The parameters for creating a counter assigned to this module when it is created in each 
scanning script. The counter requires the use of two of the four available channels on this cDAQ card. One channel is for input (from the APD on the
optical bench), and one channel for the hardware-timed clock to "run" the counter. Currently (as of 081522) only one counter (analog input) channel 
can be created at a time. In order to change this functionality, the source code for this contributed drivers package must be edited. And 2. NI-9263 
(four analog output channels (-10 V to +10 V). This module's info can be found at: 
https://www.ni.com/en-us/support/model.ni-9263.htmlhttps://www.ni.com/en-us/support/model.ni-9263.html. This module is used for programmable digital 
output. This functionality is largely used for controlling the servo motors on the ThorLab's galvo (by writing a sequence of voltages to them). Here,
all available digital output channels are created in one line: `ni_9263_ao_channels` where a single channel can be accessed at a time. The main 
function of this application is running confocal scanning scripts in different planes (the XY plane at a specified Z height, the XZ plane across a 
specified z-range, and the YZ plane at a specified X range). The GUI is split into 2 main setions, a left-half and a right-half section. The left-hand 
section contains the input controls for running a specific scanning script. On the bottom of the left-hand section is an option to save the dataset 
that was just scan and plotted on the right. Additionally, a text output box is present at the bottom of the left_window section. This is where the 
user-entered scanning parameters are printed to. Here, the use has the option to select and then copy the scan parameter info for keeping a record. 
there is also the option (currently commented out) to print this same info w/ its same format to the terminal. The right-hand section contains the 
output plot from a run scan. Include info more: this application implements live plotting for each scanning option.

For the user of this application:
1. This application was built within a virtual environment on this computer called: "NV_control" it can be found through the
user/miniconda3/envs/NV_control. Within this enviornment THIS file and additonal referenec packages can be found
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
import nidaqmx
import numpy
import pyqtgraph as pg
from datetime import date

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

############################# OLD ####################################
# def clear_plot():
#     print("Hello world!")
#     plot_res = 3.89
#     self.sc = MplCanvas(self, width = plot_res, height = plot_res, dpi = 110)                         # changing dpi does something to scale of figure
#     self.sc.move(2, 2)
#     self.sc.setParent(right_window)
#     self.sc.axes.xaxis.set_tick_params(labelsize = 8)
#     self.sc.axes.yaxis.set_tick_params(labelsize = 8)
###################################################################

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

############################################################################# "Chid" class #######################################################################
class Child(QtWidgets.QWidget):#, **kwargs): # kwargs needed?

    # ?
    def __init__(self, parent = None):#, **kwargs): # kwargs needed?
        super().__init__(parent)
        hbox = QHBoxLayout(self)
        plotwin = pg.GraphicsLayoutWidget(show = True)

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

        ########### XY scanning #############
        # print/display XY scan parameters fnc
        def print_XY_scan_parameters_fnc(self, parent = Child): # this fnc does...

            # print("XY_SCAN PARAMETERS/INFO: ", end = "")                                    # this prints to the terminal
            # print("XY_scan resolution = %d, " % int(xy_scan_resolution_qlineedit.text()), end = "")
            # print("XY_scan counter read time = %2f, " % round(float(xy_scan_read_time_qlineedit.text()), 2), end = "")
            # print("XY_scan min x driving voltage = %2f, " % float(xy_scan_x_voltage_min_qlineedit.text()), end = "")
            # print("XY_scan max x driving voltage = %2f, " % float(xy_scan_x_voltage_max_qlineedit.text()), end = "")
            # print("XY_scan min y driving voltage = %2f, " % float(xy_scan_y_voltage_min_qlineedit.text()), end = "")
            # print("XY_scan max y driving voltage = %2f, " % float(xy_scan_y_voltage_max_qlineedit.text()), end = "")
            # print("XY_scan z-piezo driving voltage = %2f." % float(xy_scan_z_piezo_voltage_qlineedit.text()))

            # need to clear text box first
            self.parameters_dsiplay_text_box.clear()

            # this prints to the QTextBox in the left_window. The output of the user-selected scan parameters is printed below
            self.parameters_dsiplay_text_box.setPlainText(
                                                        "XY_SCAN PARAMETERS/INFO:\n"
                                                        "XY_scan resolution = " + str(int(xy_scan_resolution_qlineedit.text())) + "\n"
                                                        "XY_scan counter read time = " + str(float(xy_scan_read_time_qlineedit.text())) + "\n"
                                                        "XY_scan min x driving voltage = " + str(float(xy_scan_x_voltage_min_qlineedit.text())) + "\n"
                                                        "XY_scan max x driving voltage = " + str(float(xy_scan_x_voltage_max_qlineedit.text())) + "\n"
                                                        "XY_scan min y driving voltage = " + str(float(xy_scan_y_voltage_min_qlineedit.text())) + "\n"
                                                        "XY_scan max y driving voltage = " + str(float(xy_scan_y_voltage_max_qlineedit.text())) + "\n"
                                                        "XY_scan z-piezo driving voltage = " + str(float(xy_scan_z_piezo_voltage_qlineedit.text()))
                                                        )
        
        # run XY scanning function
        def run_xy_scan_fnc(): # this fnc runs the xy_scan per the user-entered parameters in the xy_scan qlineedits

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

            # naming the instrument
            scan_galvo_card_name = "cDAQ1Mod2"

            # dictionary of analog output channels
            scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)}

            # defining the instrument (ni_9263)
            scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels)

            ############################################################################### def other variables #####################################################################

            ################### setting variales and array ####################

            # counter read time:
            scan_counter_acquisition_time = float(xy_scan_read_time_qlineedit.text())                          # note a reading time <0.05 s is likely too short

            # def
            grid_size = int(xy_scan_resolution_qlineedit.text())
            grid_size_x = grid_size_y = grid_size

            # def the initial driving voltage for the x-mirror
            initial_x_driving_voltage = round(float(xy_scan_x_voltage_min_qlineedit.text()), 2)
            # def the initial driving voltage for the y-mirror
            initial_y_driving_voltage = round(float(xy_scan_y_voltage_min_qlineedit.text()), 2)

            z_piezo_set_voltage = round(float(xy_scan_z_piezo_voltage_qlineedit.text()), 2)

            # setup parameter for the program to use based on defined variables
            x_driving_voltage_to_change = initial_x_driving_voltage
            y_driving_voltage_to_change = initial_y_driving_voltage

            # set desired end mirror voltages
            desired_end_x_mirror_voltage = round(float(xy_scan_x_voltage_max_qlineedit.text()), 2)
            desired_end_y_mirror_voltage = round(float(xy_scan_y_voltage_max_qlineedit.text()), 2)

            # def the internal stepping voltages based on user-entered settings above
            x_drive_voltage_step = ((np.absolute(initial_x_driving_voltage)) + (desired_end_x_mirror_voltage)) / grid_size_x
            y_drive_voltage_step = ((np.absolute(initial_y_driving_voltage)) + (desired_end_y_mirror_voltage)) / grid_size_y

            # create dataset to populate
            global xy_scan_data_array
            xy_scan_data_array = np.zeros((grid_size_x, grid_size_y))
            global most_recent_data_array
            most_recent_data_array = np.zeros((grid_size_x, grid_size_y))

            # initializing variables to store read counter info.
            output_value = 0
            counter_value = 0

            ################### resetting position of mirrors ####################

            scan_galvo.voltage_cdaq1mod2ao0(initial_x_driving_voltage) # this is for the x-mirror

            scan_galvo.voltage_cdaq1mod2ao1(initial_y_driving_voltage) # this is for the y-mirror

            scan_galvo.voltage_cdaq1mod2ao2(z_piezo_set_voltage) # setting the z_piezo to z_piezo_set_voltage microns

            ########################################################### setting up NI-DAQmx tasks and channels for counting ########################################################

            with nidaqmx.Task() as task1, nidaqmx.Task() as counter_output_task: # this defines 2 NI-DAQmx tasks (one for the counter and one for the counter's clock)

                # adding dig pulse train chan
                counter_output_task.co_channels.add_co_pulse_chan_freq(
                    counter = "cDAQ1Mod1/ctr1",
                    name_to_assign_to_channel = "",
                    units = nidaqmx.constants.FrequencyUnits.HZ,
                    idle_state = nidaqmx.constants.Level.LOW,
                    initial_delay = 0.0,
                    freq = 1000,
                    duty_cycle = 0.5
                    )

                # cfg implict timing
                counter_output_task.timing.cfg_implicit_timing(
                    sample_mode = AcquisitionType.CONTINUOUS,
                    samps_per_chan = 1
                    )

                # adding count egdes chan
                task1.ci_channels.add_ci_count_edges_chan(
                    counter = "cDAQ1Mod1/ctr0",
                    name_to_assign_to_channel = "",
                    edge = nidaqmx.constants.Edge.RISING,
                    initial_count = 0,
                    count_direction = nidaqmx.constants.CountDirection.COUNT_UP
                    )

                # cfg sample clk timing
                task1.timing.cfg_samp_clk_timing(
                    rate = 1000,
                    source = "/cDAQ1/Ctr1InternalOutput",
                    active_edge = nidaqmx.constants.Edge.RISING,
                    sample_mode = AcquisitionType.CONTINUOUS,
                    samps_per_chan = 1
                    )
                
                counter_output_task.start() # this starts the counter NI-DAQmx task
                task1.start() # this starts the hardware-based internal clock NI-DAQmx task
                                
            ######################################################################## X and Y scanning #########################################################################

                for f in range(grid_size_y): # this loops for rows (y)

                    for k in range(grid_size_x): # this loops for columns (x)

                        ################## important section #################

                        for my_var_not_named_i in range(int(scan_counter_acquisition_time * 1000)): # this reads/lets the counter accumulate for the set time and returns value
                            counter_value = task1.read()
                            output_value += counter_value

                        if f % 2 != 0: # this loop populates the created xy_scan_data_array (the if else strucuture is present bc of the snaking scanning pattern)
                            xy_scan_data_array[f][((-k) + 1)] = (output_value - np.sum(xy_scan_data_array)) # add counter result to data array                           # saving (082422)
                            output_value == 0
                            counter_value == 0
                        else:
                            if f == 0 and k == 0:
                                xy_scan_data_array[0][0] = output_value # add counter result to data array                                                      # saving
                            else:
                                xy_scan_data_array[f][k] = (output_value - np.sum(xy_scan_data_array)) # add counter result to data array                               # saving (082422)
                        output_value = 0
                        counter_value = 0

                        ############# end important section ############

                        if f % 2 == 0: # this loop adjusts for sweeping back and forth along each alternating row

                            if k < (grid_size_x - 1):

                                x_driving_voltage_to_change += x_drive_voltage_step # increment drive voltage forwards
                                scan_galvo.voltage_cdaq1mod2ao0(x_driving_voltage_to_change) # step x motor

                            else:
                                break
                        else:
                            if k < (grid_size_x - 1):

                                x_driving_voltage_to_change -= x_drive_voltage_step # increment drive voltage backwards
                                scan_galvo.voltage_cdaq1mod2ao0(x_driving_voltage_to_change) # step x motor

                            else:
                                break

                    ##################### updating plot section ####################
                    # self.sc.axes.cla() # this does not seem to be needed

                    self.sc.axes.pcolormesh(xy_scan_data_array, cmap = "inferno")
                    self.sc.axes.set_xticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_x_driving_voltage, int((initial_x_driving_voltage + desired_end_x_mirror_voltage) / 2), desired_end_x_mirror_voltage])
                    self.sc.axes.set_yticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_y_driving_voltage, int((initial_y_driving_voltage + desired_end_y_mirror_voltage) / 2), desired_end_y_mirror_voltage])
                    self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.yaxis.set_tick_params(labelsize = 8)

                    self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
                    self.sc.axes.set_ylabel("y_mirror_driving_voltage_(V)", fontsize = 8, labelpad = -9)

                    self.sc.fig.canvas.draw()
                    self.sc.fig.canvas.flush_events() #this line is very important
                    ##################### updating plot section ####################

                    if f < (grid_size_y - 1): # this loop prevents from scanning an upper undesired row

                        y_driving_voltage_to_change += y_drive_voltage_step # increment drive voltage
                        scan_galvo.voltage_cdaq1mod2ao1(y_driving_voltage_to_change) # step y motor

                    else:
                        break

                counter_output_task.stop() # this stops the counter NI-DAQmx task - free-ing the reserved cDAQ card resources
                task1.stop() # this stops the hardware-based internal clock NI-DAQmx task - free-ing the reserved cDAQ card resources

            scan_galvo.close() # this is where reeated scanning hinges on and needs to be re-implemented

            ############################################################### end scanning script #############################################################################################

            ############################################### plotting XY scan data in plot ####################################################
            
            # self.sc.fig.clear() # is this needed after implementing live plot updating?
            self.sc.axes.cla()

            plot = self.sc.axes.pcolormesh(xy_scan_data_array, cmap = "inferno")
            self.sc.axes.set_xticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_x_driving_voltage, int((initial_x_driving_voltage + desired_end_x_mirror_voltage) / 2), desired_end_x_mirror_voltage])
            self.sc.axes.set_yticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_y_driving_voltage, int((initial_y_driving_voltage + desired_end_y_mirror_voltage) / 2), desired_end_y_mirror_voltage])
            self.xy_scan_plot_colorbar = self.sc.fig.colorbar(plot, ax = self.sc.axes, pad = 0.02, aspect = 15)
            self.xy_scan_plot_colorbar.formatter.set_powerlimits((0, 0))
            self.sc.axes.xaxis.set_tick_params(labelsize = 8)
            self.sc.axes.yaxis.set_tick_params(labelsize = 8)
            self.xy_scan_plot_colorbar.ax.tick_params(labelsize = 7)
            self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
            self.sc.axes.set_ylabel("y_mirror_driving_voltage_(V)", fontsize = 8, labelpad = -9)
            self.sc.axes.set_title("XY_scan_%s_z-piezo@%d_microns" % (todays_date, int((z_piezo_set_voltage * 10))), fontsize = 8)
            self.sc.draw()

            most_recent_data_array = xy_scan_data_array # make temp holding global data_array the same as xy_scan_data_array

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

            res_min_condition = 20 # set the min allowed resolution for scanning
            res_max_condition = 900 # set the max allowed resolution for scanning

            xy_scan_resolution_test_condition = False # define resolution validation bool for xy scan

            while xy_scan_resolution_test_condition is False: # this initiates checking the resolution parameter

                # checking for out of bounds of min and max conditions above
                # TODO: or negative or not a number or too large
                if int(xy_scan_resolution_qlineedit.text()) < res_min_condition or int(xy_scan_resolution_qlineedit.text()) > res_max_condition:

                    display_resolution_error_window_fnc() # call the error message pop-up window
                    break # exit the checking loop: failed

                # if parameter is in bounds; run scan
                elif int(xy_scan_resolution_qlineedit.text()) >= res_min_condition and int(xy_scan_resolution_qlineedit.text()) <= res_max_condition:

                    xy_scan_resolution_test_condition == True
                    print_XY_scan_parameters_fnc(self) # call the print user-entered parameters fnc
                    run_xy_scan_fnc() # call the run xy scan method fnc
                    break # exit the checking loop: passed

        ########### XZ scanning #############
        # print_XZ_scan_parameters_fnc
        def print_XZ_scan_parameters_fnc(self, parent = Child): # this fnc does...

            # print("XZ_SCAN PARAMETERS/INFO: ", end = "")                                    # this prints to the terminal
            # print("XZ_scan resolution = %d, " % int(xz_scan_resolution_qlineedit.text()), end = "")
            # print("XZ_scan counter read time = %2f, " % round(float(xz_scan_read_time_qlineedit.text()), 2), end = "")
            # print("XZ_scan min x driving voltage = %2f, " % float(xz_scan_x_voltage_min_qlineedit.text()), end = "")
            # print("XZ_scan max x driving voltage = %2f, " % float(xz_scan_x_voltage_max_qlineedit.text()), end = "")
            # print("XZ_scan y driving voltage = %2f, " % float(xz_scan_y_voltage_qlineedit.text()), end = "")
            # print("XZ_scan z-piezo min driving voltage = %2f, " % float(xz_scan_z_piezo_min_voltage_qlineedit.text()), end = "")
            # print("XZ_scan z-piezo max driving voltage = %2f." % float(xz_scan_z_piezo_max_voltage_qlineedit.text()))


            # need to clear text box first
            self.parameters_dsiplay_text_box.clear()

            # this prints to the QTextBox in the left_window. The output of the user-selected scan parameters is printed below
            self.parameters_dsiplay_text_box.setPlainText(
                                                        "XZ_SCAN PARAMETERS/INFO:\n"
                                                        "XZ_scan resolution = " + str(int(xz_scan_resolution_qlineedit.text())) + "\n"
                                                        "XZ_scan counter read time = " + str(float(xz_scan_read_time_qlineedit.text())) + "\n"
                                                        "XZ_scan min x driving voltage = " + str(float(xz_scan_x_voltage_min_qlineedit.text())) + "\n"
                                                        "XZ_scan max x driving voltage = " + str(float(xz_scan_x_voltage_max_qlineedit.text())) + "\n"
                                                        "XZ_scan y driving voltage = " + str(float(xz_scan_y_voltage_qlineedit.text())) + "\n"
                                                        "XZ_scan z-piezo min driving voltage = " + str(float(xz_scan_z_piezo_min_voltage_qlineedit.text())) + "\n"
                                                        "XZ_scan z-piezo max driving voltage = " + str(float(xz_scan_z_piezo_max_voltage_qlineedit.text()))
                                                        )
        
        # run XZ scan fnc
        def run_xz_scan_fnc(): # this fnc runs the xy_scan per the user-entered parameters in the xy_scan qlineedits

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

            # naming the instrument
            scan_galvo_card_name = "cDAQ1Mod2"

            # dictionary of analog output channels
            scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)}

            # defining the instrument (ni_9263)
            scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels)

            ############################################################################### def other variables #####################################################################

            ################### setting variales and array ####################

            # counter read time:
            scan_counter_acquisition_time = float(xz_scan_read_time_qlineedit.text())                                     # note a reading time <0.05 s is likely too short

            # def grid_size of scan (resolution)
            grid_size = int(xz_scan_resolution_qlineedit.text())

            # def the initial driving voltage for the x-mirror
            initial_x_driving_voltage = float(xz_scan_x_voltage_min_qlineedit.text())
            # def the initial driving voltage for the y-mirror
            initial_y_driving_voltage = float(xz_scan_y_voltage_qlineedit.text())
            #
            initial_z_piezo_driving_voltage = float(xz_scan_z_piezo_min_voltage_qlineedit.text())

            # set desired end mirror voltages
            desired_end_x_mirror_voltage = float(xz_scan_x_voltage_max_qlineedit.text())
            desired_end_z_piezo_voltage = float(xz_scan_z_piezo_max_voltage_qlineedit.text())

            # setup parameter for the program to use based on defined variables
            x_driving_voltage_to_change = initial_x_driving_voltage
            z_piezo_driving_voltage_to_change = initial_z_piezo_driving_voltage

            # def the internal stepping voltages based on user-entered settings above
            x_drive_voltage_step = ((np.absolute(initial_x_driving_voltage)) + (desired_end_x_mirror_voltage)) / grid_size
            z_piezo_drive_voltage_step = (initial_z_piezo_driving_voltage + (desired_end_z_piezo_voltage)) / grid_size

            # create dataset to populate
            global xz_scan_data_array
            xz_scan_data_array = np.zeros((grid_size, grid_size))
            global most_recent_data_array
            most_recent_data_array = np.zeros((grid_size, grid_size))

            # initializing variables to store read counter info.
            output_value = 0
            counter_value = 0

            ################### resetting position of mirrors ####################

            scan_galvo.voltage_cdaq1mod2ao0(initial_x_driving_voltage) # this is for setting the initial x voltage

            scan_galvo.voltage_cdaq1mod2ao1(initial_y_driving_voltage) # this is for fixing the x-mirror at defined setting (voltage)

            scan_galvo.voltage_cdaq1mod2ao2(initial_z_piezo_driving_voltage) # starting z-piezo at 0.0 V (0 microns); this variable will increase during scanning

            ########################################################### setting up NI-DAQmx tasks and channels for counting ########################################################

            with nidaqmx.Task() as task1, nidaqmx.Task() as counter_output_task: # this defines 2 NI-DAQmx tasks (one for the counter and one for the counter's clock)

                # adding dig pulse train chan
                counter_output_task.co_channels.add_co_pulse_chan_freq(
                    counter = "cDAQ1Mod1/ctr1",
                    name_to_assign_to_channel = "",
                    units = FrequencyUnits.HZ,
                    idle_state = nidaqmx.constants.Level.LOW,
                    initial_delay = 0.0,
                    freq = 1000,
                    duty_cycle = 0.5
                    )

                # cfg implict timing
                counter_output_task.timing.cfg_implicit_timing(
                    sample_mode = AcquisitionType.CONTINUOUS,
                    samps_per_chan = 1
                    )

                # adding count egdes chan
                task1.ci_channels.add_ci_count_edges_chan(
                    counter = "cDAQ1Mod1/ctr0",
                    name_to_assign_to_channel = "",
                    edge = Edge.RISING,
                    initial_count = 0,
                    count_direction = CountDirection.COUNT_UP
                    )

                # cfg sample clk timing
                task1.timing.cfg_samp_clk_timing(
                    rate = 1000,
                    source = "/cDAQ1/Ctr1InternalOutput",
                    active_edge = Edge.RISING,
                    sample_mode = AcquisitionType.CONTINUOUS,
                    samps_per_chan = 1
                    )
                
                counter_output_task.start() # this starts the counter NI-DAQmx task
                task1.start() # this starts the hardware-based internal clock NI-DAQmx task
                        
                ####################################################################### x and z scanning #########################################################################

                for f in range(grid_size): # this loops for stacks in z

                    for k in range(grid_size): # this loops for columns (varying y voltage/setting)

                        ################## important section #################

                        for my_var_not_named_i in range(int(scan_counter_acquisition_time * 1000)): # this reads/lets the counter accumulate for the set time and returns value
                            counter_value = task1.read()
                            output_value += counter_value

                        if f % 2 != 0: # this loop populates the created xz_scan_data_array (the if else strucuture is present bc of the snaking scanning pattern)
                            xz_scan_data_array[f][((-k) + 1)] = (output_value - np.sum(xz_scan_data_array)) # add counter result to data array                           # saving (082422)
                            output_value == 0
                            counter_value == 0
                        else:
                            if f == 0 and k == 0:
                                xz_scan_data_array[0][0] = output_value # add counter result to data array                                                      # saving
                            else:
                                xz_scan_data_array[f][k] = (output_value - np.sum(xz_scan_data_array)) # add counter result to data array                               # saving (082422)
                        output_value = 0
                        counter_value = 0

                        ############# end important section ############

                        if f % 2 == 0: # this loop adjusts for sweeping back and forth along each alternating row
                            if k < (grid_size - 1):
                                x_driving_voltage_to_change += x_drive_voltage_step # increment y mirro drive voltage forwards
                                scan_galvo.voltage_cdaq1mod2ao0(x_driving_voltage_to_change) # step y mirror
                            else:
                                break
                        else:
                            if k < (grid_size - 1):
                                x_driving_voltage_to_change -= x_drive_voltage_step # increment y mirror drive voltage backwards
                                scan_galvo.voltage_cdaq1mod2ao0(x_driving_voltage_to_change) # step y mirror
                            else:
                                break
                    
                    ##################### updating plot section ####################
                    # self.sc.axes.cla() # this does not seem to be needed
                    self.sc.axes.pcolormesh(xz_scan_data_array, cmap = "inferno")
                    self.sc.axes.set_xticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_x_driving_voltage, int((initial_x_driving_voltage + desired_end_x_mirror_voltage) / 2), desired_end_x_mirror_voltage])
                    self.sc.axes.set_yticks(np.arange(0, grid_size + 10, grid_size / 2), [int(initial_z_piezo_driving_voltage * 10), (int(initial_z_piezo_driving_voltage * 10) + int(desired_end_z_piezo_voltage * 10)) / 2, int(desired_end_z_piezo_voltage * 10)])
                    self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
                    self.sc.axes.set_ylabel("objective_z-piezo_height_(microns)", fontsize = 8)
                    self.sc.fig.canvas.draw()
                    self.sc.fig.canvas.flush_events()
                    ##################### updating plot section ####################

                    if f < (grid_size - 1): # this loop
                        z_piezo_driving_voltage_to_change += z_piezo_drive_voltage_step # increment drive voltage
                        scan_galvo.voltage_cdaq1mod2ao2(z_piezo_driving_voltage_to_change) # step z
                    else:
                        break

                counter_output_task.stop() # this stops the counter NI-DAQmx task - free-ing the reserved cDAQ card resources
                task1.stop() # this stops the hardware-based internal clock NI-DAQmx task - free-ing the reserved cDAQ card resources

            scan_galvo.close()
            ############################################ end XZ scanning script #################################################

            ############################################### plotting XZ scan data in plot ####################################################
            self.sc.axes.cla()
            xz_scan_plot = self.sc.axes.pcolormesh(xz_scan_data_array, cmap = "inferno")
            xz_scan_plot_colorbar = self.sc.fig.colorbar(xz_scan_plot, ax = self.sc.axes, pad = 0.02, aspect = 15)
            xz_scan_plot_colorbar.formatter.set_powerlimits((0, 0))
            self.sc.axes.xaxis.set_tick_params(labelsize = 8)
            self.sc.axes.yaxis.set_tick_params(labelsize = 8)
            xz_scan_plot_colorbar.ax.tick_params(labelsize = 7)
            self.sc.axes.set_xticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_x_driving_voltage, int((initial_x_driving_voltage + desired_end_x_mirror_voltage) / 2), desired_end_x_mirror_voltage])
            self.sc.axes.set_yticks(np.arange(0, grid_size + 10, grid_size / 2), [int(initial_z_piezo_driving_voltage * 10), (int(initial_z_piezo_driving_voltage * 10) + int(desired_end_z_piezo_voltage * 10)) / 2, int(desired_end_z_piezo_voltage * 10)])
            self.sc.axes.set_xlabel("x_mirror_driving_voltage_(V)", fontsize = 8)
            self.sc.axes.set_ylabel("objective_z-piezo_height_(microns)", fontsize = 8)
            self.sc.axes.set_title("XZ_scan_%s_y_voltage=%f_V" % (todays_date, initial_y_driving_voltage), fontsize = 8)
            self.sc.draw()

            most_recent_data_array = xz_scan_data_array

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

            res_min_condition = 20 # set the min allowed resolution for scanning
            res_max_condition = 900 # set the max allowed resolution for scanning

            xz_scan_resolution_test_condition = False # define resolution validation bool for xz scan

            while xz_scan_resolution_test_condition is False: # this initiates checking the resolution parameter

                # checking for out of bounds of min and max conditions above
                if int(xz_scan_resolution_qlineedit.text()) < res_min_condition or int(xz_scan_resolution_qlineedit.text()) > res_max_condition: # TODO: or negative or not a number or too large

                    display_resolution_error_window_fnc() # call the error message pop-up window
                    break # exit the checking loop: failed

                # if parameter is in bounds; run scan
                elif int(xz_scan_resolution_qlineedit.text()) >= res_min_condition and int(xz_scan_resolution_qlineedit.text()) <= res_max_condition:

                    xz_scan_resolution_test_condition == True
                    print_XZ_scan_parameters_fnc(self) # call the print user-entered parameters fnc
                    run_xz_scan_fnc() # call the run xz scan method fnc
                    break # exit the checking loop: passed

        ########################## YZ scanning ############################

        # print_YZ_scan_parameters_fnc
        def print_YZ_scan_parameters_fnc(self, parent = Child): # this fnc does...

            # print("YZ_SCAN PARAMETERS/INFO: ", end = "")
            # print("YZ_scan resolution = %d, " % int(yz_scan_resolution_qlineedit.text()), end = "")
            # print("YZ_scan counter read time = %2f, " % round(float(yz_scan_read_time_qlineedit.text()), 2), end = "")
            # print("YZ_scan min Y driving voltage = %2f, " % float(yz_scan_y_voltage_min_qlineedit.text()), end = "")
            # print("YZ_scan max Y driving voltage = %2f, " % float(yz_scan_y_voltage_max_qlineedit.text()), end = "")
            # print("YZ_scan X driving voltage = %2f, " % float(yz_scan_x_voltage_qlineedit.text()), end = "")
            # print("YZ_scan z-piezo min driving voltage = %2f, " % float(yz_scan_z_piezo_min_voltage_qlineedit.text()), end = "")
            # print("YZ_scan z-piezo max driving voltage = %2f." % float(yz_scan_z_piezo_max_voltage_qlineedit.text()))

            # need to clear text box first
            self.parameters_dsiplay_text_box.clear()

            # this prints to the QTextBox in the left_window. The output of the user-selected scan parameters is printed below
            self.parameters_dsiplay_text_box.setPlainText(
                        "YZ_SCAN PARAMETERS/INFO:\n"
                        "YZ_scan resolution = " + str(int(yz_scan_resolution_qlineedit.text())) + "\n"
                        "YZ_scan counter read time = " + str(float(yz_scan_read_time_qlineedit.text())) + "\n"
                        "YZ_scan min Y driving voltage = " + str(float(yz_scan_y_voltage_min_qlineedit.text())) + "\n"
                        "YZ_scan max Y driving voltage = " + str(float(yz_scan_y_voltage_max_qlineedit.text())) + "\n"
                        "YZ_scan X driving voltage = " + str(float(yz_scan_x_voltage_qlineedit.text())) + "\n'"
                        "YZ_scan z-piezo min driving voltage = " + str(float(yz_scan_z_piezo_min_voltage_qlineedit.text())) + "\n"
                        "YZ_scan z-piezo max driving voltage = " + str(float(yz_scan_z_piezo_max_voltage_qlineedit.text()))
                                                        )
                
        # run YZ scan function
        def run_yz_scan_fnc(): # this fnc runs the YZ scan script
            # print("YZ scan started")

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

            # naming the instrument
            scan_galvo_card_name = "cDAQ1Mod2"

            # dictionary of analog output channels
            scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)}

            # defining the instrument (ni_9263)
            scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels)

            ############################################################################### def other variables #####################################################################

            ################### setting variales and array ####################

            # counter read time:
            scan_counter_acquisition_time = 0.05                                                                               # note a reading time <0.05 s is likely too short

            # def
            grid_size = int(yz_scan_resolution_qlineedit.text())

            # def the initial driving voltage for the x-mirror
            initial_x_driving_voltage = float(yz_scan_x_voltage_qlineedit.text())

            # def the initial driving voltage for the y-mirror
            initial_y_driving_voltage = float(yz_scan_y_voltage_min_qlineedit.text())
            #
            initial_z_piezo_driving_voltage = float(yz_scan_z_piezo_min_voltage_qlineedit.text())

            # setup parameter for the program to use based on defined variables
            y_driving_voltage_to_change = initial_y_driving_voltage
            z_piezo_driving_voltage_to_change = initial_z_piezo_driving_voltage

            # set desired end mirror voltages
            desired_end_y_mirror_voltage = float(yz_scan_y_voltage_max_qlineedit.text())
            desired_end_z_piezo_voltage = float(yz_scan_z_piezo_max_voltage_qlineedit.text())

            # def the internal stepping voltages based on user-entered settings above
            y_drive_voltage_step = ((np.absolute(initial_y_driving_voltage)) + (desired_end_y_mirror_voltage)) / grid_size
            z_piezo_drive_voltage_step = (initial_z_piezo_driving_voltage + (desired_end_z_piezo_voltage)) / grid_size

            # create dataset to populate
            global yz_scan_data_array
            yz_scan_data_array = np.zeros((grid_size, grid_size))
            global most_recent_data_array
            most_recent_data_array = np.zeros((grid_size, grid_size))

            # initializing variables to store read counter info.
            output_value = 0
            counter_value = 0

            ################### resetting position of mirrors ####################

            scan_galvo.voltage_cdaq1mod2ao0(initial_x_driving_voltage) # this is for fixing the x-mirror at defined setting (voltage)

            scan_galvo.voltage_cdaq1mod2ao1(initial_y_driving_voltage) # this is for setting the initial x voltage

            scan_galvo.voltage_cdaq1mod2ao2(initial_z_piezo_driving_voltage) # starting z-piezo at 0.0 V (0 microns); this variable will increase during scanning

            ########################################################### setting up NI-DAQmx tasks and channels for counting ########################################################

            with nidaqmx.Task() as task1, nidaqmx.Task() as counter_output_task: # this defines 2 NI-DAQmx tasks (one for the counter and one for the counter's clock)

                # adding dig pulse train chan
                counter_output_task.co_channels.add_co_pulse_chan_freq(
                    counter = "cDAQ1Mod1/ctr1",
                    name_to_assign_to_channel = "",
                    units = FrequencyUnits.HZ,
                    idle_state = nidaqmx.constants.Level.LOW,
                    initial_delay = 0.0,
                    freq = 1000,
                    duty_cycle = 0.5
                    )

                # cfg implict timing
                counter_output_task.timing.cfg_implicit_timing(
                    sample_mode = AcquisitionType.CONTINUOUS,
                    samps_per_chan = 1
                    )

                # adding count egdes chan
                task1.ci_channels.add_ci_count_edges_chan(
                    counter = "cDAQ1Mod1/ctr0",
                    name_to_assign_to_channel = "",
                    edge = Edge.RISING,
                    initial_count = 0,
                    count_direction = CountDirection.COUNT_UP
                    )

                # cfg sample clk timing
                task1.timing.cfg_samp_clk_timing(
                    rate = 1000,
                    source = "/cDAQ1/Ctr1InternalOutput",
                    active_edge = Edge.RISING,
                    sample_mode = AcquisitionType.CONTINUOUS,
                    samps_per_chan = 1
                    )
                
                counter_output_task.start() # this starts the counter NI-DAQmx task
                task1.start() # this starts the hardware-based internal clock NI-DAQmx task
                        
            ####################################################################### x and z scanning #########################################################################

                for f in range(grid_size): # this loops for stacks in z

                    for k in range(grid_size): # this loops for columns (varying y voltage/setting)

                        ################## important section #################

                        for my_var_not_named_i in range(int(scan_counter_acquisition_time * 1000)): # this reads/lets the counter accumulate for the set time and returns value
                            counter_value = task1.read()
                            output_value += counter_value

                        if f % 2 != 0: # this loop populates the created yz_scan_data_array (the if else strucuture is present bc of the snaking scanning pattern)
                            yz_scan_data_array[f][((-k) + 1)] = (output_value - np.sum(yz_scan_data_array)) # add counter result to data array                           # saving (082422)
                            output_value == 0
                            counter_value == 0
                        else:
                            if f == 0 and k == 0:
                                yz_scan_data_array[0][0] = output_value # add counter result to data array                                                      # saving
                            else:
                                yz_scan_data_array[f][k] = (output_value - np.sum(yz_scan_data_array)) # add counter result to data array                               # saving (082422)
                        output_value = 0
                        counter_value = 0

                        ############# end important section ############

                        if f % 2 == 0: # this loop adjusts for sweeping back and forth along each alternating row
                            if k < (grid_size - 1):
                                y_driving_voltage_to_change += y_drive_voltage_step # increment y mirro drive voltage forwards
                                scan_galvo.voltage_cdaq1mod2ao1(y_driving_voltage_to_change) # step y mirror
                            else:
                                break
                        else:
                            if k < (grid_size - 1):
                                y_driving_voltage_to_change -= y_drive_voltage_step # increment y mirror drive voltage backwards
                                scan_galvo.voltage_cdaq1mod2ao1(y_driving_voltage_to_change) # step y mirror
                            else:
                                break
                    
                    ##################### updating plot section ####################
                    # self.sc.axes.cla() # this does not seem to be needed
                    self.sc.axes.pcolormesh(yz_scan_data_array, cmap = "inferno")
                    self.sc.axes.set_xticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_y_driving_voltage, int((initial_y_driving_voltage + desired_end_y_mirror_voltage) / 2), desired_end_y_mirror_voltage])
                    self.sc.axes.set_yticks(np.arange(0, grid_size + 10, grid_size / 2), [int(initial_z_piezo_driving_voltage * 10), (int(initial_z_piezo_driving_voltage * 10) + int(desired_end_z_piezo_voltage * 10)) / 2, int(desired_end_z_piezo_voltage * 10)])
                    self.sc.axes.xaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.yaxis.set_tick_params(labelsize = 8)
                    self.sc.axes.set_xlabel("y_mirror_driving_voltage_(V)", fontsize = 8)
                    self.sc.axes.set_ylabel("objective_z-piezo_height_(microns)", fontsize = 8)
                    self.sc.fig.canvas.draw()
                    self.sc.fig.canvas.flush_events()
                    ##################### updating plot section ####################

                    if f < (grid_size - 1): # this loop
                        z_piezo_driving_voltage_to_change += z_piezo_drive_voltage_step # increment drive voltage
                        scan_galvo.voltage_cdaq1mod2ao2(z_piezo_driving_voltage_to_change) # step z
                    else:
                        break

                counter_output_task.stop() # this stops the counter NI-DAQmx task - free-ing the reserved cDAQ card resources
                task1.stop() # this stops the hardware-based internal clock NI-DAQmx task - free-ing the reserved cDAQ card resources

            scan_galvo.close()

            ############################################### plotting YZ scan data in plot ####################################################
            self.sc.axes.cla()
            yz_scan_plot = self.sc.axes.pcolormesh(yz_scan_data_array, cmap = "inferno")
            yz_scan_plot_colorbar = self.sc.fig.colorbar(yz_scan_plot, ax = self.sc.axes, pad = 0.02, aspect = 15)
            yz_scan_plot_colorbar.formatter.set_powerlimits((0, 0))
            self.sc.axes.xaxis.set_tick_params(labelsize = 8)
            self.sc.axes.yaxis.set_tick_params(labelsize = 8)
            yz_scan_plot_colorbar.ax.tick_params(labelsize = 7)
            self.sc.axes.set_xticks(np.arange(0, grid_size + 10, grid_size / 2), [initial_y_driving_voltage, int((initial_y_driving_voltage + desired_end_y_mirror_voltage) / 2), desired_end_y_mirror_voltage])
            self.sc.axes.set_yticks(np.arange(0, grid_size + 10, grid_size / 2), [int(initial_z_piezo_driving_voltage * 10), (int(initial_z_piezo_driving_voltage * 10) + int(desired_end_z_piezo_voltage * 10)) / 2, int(desired_end_z_piezo_voltage * 10)])
            self.sc.axes.set_xlabel("y_mirror_driving_voltage_(V)", fontsize = 8)
            self.sc.axes.set_ylabel("objective_z-piezo_height_(microns)", fontsize = 8)
            self.sc.axes.set_title("YZ_scan_%s_x_voltage=%f_V" % (todays_date, initial_x_driving_voltage), fontsize = 8)
            self.sc.draw()
            most_recent_data_array = yz_scan_data_array
    
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

            res_min_condition = 20 # set the min allowed resolution for scanning
            res_max_condition = 900 # set the max allowed resolution for scanning

            yz_scan_resolution_test_condition = False # define resolution validation bool for yz scan

            while yz_scan_resolution_test_condition is False: # this initiates checking the resolution parameter

                # checking for out of bounds of min and max conditions above
                if int(yz_scan_resolution_qlineedit.text()) < res_min_condition or int(yz_scan_resolution_qlineedit.text()) > res_max_condition: # TODO: or negative or not a number or too large

                    display_resolution_error_window_fnc() # call the error message pop-up window
                    break # exit the checking loop: failed

                # if parameter is in bounds; run scan
                elif int(yz_scan_resolution_qlineedit.text()) >= res_min_condition and int(yz_scan_resolution_qlineedit.text()) <= res_max_condition:

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

        self.parameters_dsiplay_text_box = QTextEdit(left_window)
        self.parameters_dsiplay_text_box.resize(320, 113)
        self.parameters_dsiplay_text_box.move(10, 309)
        self.parameters_dsiplay_text_box.setParent(left_window)

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
        ############ bendegin save data section ###############

        # attempt at displaying Tex text as a QLabel
        pageSource = """
                    <html>

                    <head>
                    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_HTMLorMML">                     
                    </script>
                    </head>

                    <body>
                    <p><mathjax style = "font-size:1.0em">$$normal_{subscript}$$</mathjax></p>
                    </body>
                    
                    </html>
                    """

        # scan widgets overall "Scanning"
        scan_widget_overall = QLabel("Scanning options:")
        # scan_widget_overall = QLabel()
        # scan_widget_overall = QWebEngineView()
        # scan_widget_overall.setHtml(pageSource)
        scan_widget_overall.setFont(QFont("Times font", 9))
        scan_widget_overall.setParent(left_window)
        scan_widget_overall.move(120, 10)
        # scan_widget_overall.move(30, 10)

        indiv_scan_labels_y_height = 25

        row_y_adjust = 5

        overall_y_adjust = 10

        #########################################################################################
        ########################################### XY scanning #################################
        #########################################################################################

        xy_scan_widgets_left_x_justify = 5
        xy_scan_parameters_validated_Q = False

        # XY scan
        xy_scan_label_widget = QLabel("XY scan", self) # widget
        xy_scan_label_widget.setParent(left_window)
        xy_scan_label_widget.move(12 + 10, indiv_scan_labels_y_height + overall_y_adjust)

        # resolution
        xy_scan_resolution_widget = QLabel("Res:", self) # widget
        xy_scan_resolution_widget.setParent(left_window)
        xy_scan_resolution_widget.move(xy_scan_widgets_left_x_justify, 40 + row_y_adjust + overall_y_adjust)

        xy_scan_resolution_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_resolution_qlineedit.setParent(left_window)
        xy_scan_resolution_qlineedit.move(30, 40 + row_y_adjust + overall_y_adjust)
        xy_scan_resolution_qlineedit.resize(55, 15)
        xy_scan_resolution_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # read time
        xy_scan_read_time_widget = QLabel("APD_t:", self) # widget
        xy_scan_read_time_widget.setParent(left_window)
        xy_scan_read_time_widget.move(xy_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)

        xy_scan_read_time_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_read_time_qlineedit.setParent(left_window)
        xy_scan_read_time_qlineedit.move(40, 65 + row_y_adjust + overall_y_adjust)
        xy_scan_read_time_qlineedit.resize(45, 15)
        xy_scan_read_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
        
        # x voltage (min and max)
        xy_scan_x_voltage_min_widget = QLabel("x_V_min:", self) # widget
        xy_scan_x_voltage_min_widget.setParent(left_window)
        xy_scan_x_voltage_min_widget.move(xy_scan_widgets_left_x_justify, 90 + row_y_adjust + overall_y_adjust)

        xy_scan_x_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_x_voltage_min_qlineedit.setParent(left_window)
        xy_scan_x_voltage_min_qlineedit.move(50, 90 + row_y_adjust + overall_y_adjust)
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

        xy_scan_z_piezo_voltage_qlineedit = QLineEdit(self) # qclineedit
        xy_scan_z_piezo_voltage_qlineedit.setParent(left_window)
        xy_scan_z_piezo_voltage_qlineedit.move(30, 190 + row_y_adjust + overall_y_adjust)
        xy_scan_z_piezo_voltage_qlineedit.resize(55, 15)
        xy_scan_z_piezo_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

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

        # resolution
        xz_scan_resolution_widget = QLabel("Res:", self) # widget
        xz_scan_resolution_widget.setParent(left_window)
        xz_scan_resolution_widget.move(xz_scan_widgets_left_x_justify, 40 + row_y_adjust + overall_y_adjust)

        xz_scan_resolution_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_resolution_qlineedit.setParent(left_window)
        xz_scan_resolution_qlineedit.move(25 + xz_scan_widgets_left_x_justify, 40 + row_y_adjust + overall_y_adjust)
        xz_scan_resolution_qlineedit.resize(55, 15)
        xz_scan_resolution_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # read time
        xz_scan_read_time_widget = QLabel("APD_t:", self) # widget
        xz_scan_read_time_widget.setParent(left_window)
        xz_scan_read_time_widget.move(xz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)

        xz_scan_read_time_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_read_time_qlineedit.setParent(left_window)
        xz_scan_read_time_qlineedit.move(35 + xz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)
        xz_scan_read_time_qlineedit.resize(45, 15)
        xz_scan_read_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
        
        # x voltage (min and max)
        xz_scan_x_voltage_min_widget = QLabel("x_V_min:", self) # widget
        xz_scan_x_voltage_min_widget.setParent(left_window)
        xz_scan_x_voltage_min_widget.move(xz_scan_widgets_left_x_justify, 90 + row_y_adjust + overall_y_adjust)

        xz_scan_x_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_x_voltage_min_qlineedit.setParent(left_window)
        xz_scan_x_voltage_min_qlineedit.move(45 + xz_scan_widgets_left_x_justify, 90 + row_y_adjust + overall_y_adjust)
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

        xz_scan_z_piezo_min_voltage_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_z_piezo_min_voltage_qlineedit.setParent(left_window)
        xz_scan_z_piezo_min_voltage_qlineedit.move(45 + xz_scan_widgets_left_x_justify, 165 + row_y_adjust + overall_y_adjust)
        xz_scan_z_piezo_min_voltage_qlineedit.resize(35, 15)
        xz_scan_z_piezo_min_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        xz_scan_z_piezo_max_voltage_widget = QLabel("z_V_max:", self) # widget
        xz_scan_z_piezo_max_voltage_widget.setParent(left_window)
        xz_scan_z_piezo_max_voltage_widget.move(xz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)

        xz_scan_z_piezo_max_voltage_qlineedit = QLineEdit(self) # qclineedit
        xz_scan_z_piezo_max_voltage_qlineedit.setParent(left_window)
        xz_scan_z_piezo_max_voltage_qlineedit.move(50 + xz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)
        xz_scan_z_piezo_max_voltage_qlineedit.resize(30, 15)
        xz_scan_z_piezo_max_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

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

        # resolution
        yz_scan_resolution_widget = QLabel("Res:", self) # widget
        yz_scan_resolution_widget.setParent(left_window)
        yz_scan_resolution_widget.move(yz_scan_widgets_left_x_justify, 40 + row_y_adjust + overall_y_adjust)

        yz_scan_resolution_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_resolution_qlineedit.setParent(left_window)
        yz_scan_resolution_qlineedit.move(25 + yz_scan_widgets_left_x_justify, 40 + row_y_adjust + overall_y_adjust)
        yz_scan_resolution_qlineedit.resize(55, 15)
        yz_scan_resolution_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # read time
        yz_scan_read_time_widget = QLabel("APD_t:", self) # widget
        yz_scan_read_time_widget.setParent(left_window)
        yz_scan_read_time_widget.move(yz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)

        yz_scan_read_time_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_read_time_qlineedit.setParent(left_window)
        yz_scan_read_time_qlineedit.move(35 + yz_scan_widgets_left_x_justify, 65 + row_y_adjust + overall_y_adjust)
        yz_scan_read_time_qlineedit.resize(45, 15)
        yz_scan_read_time_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)
        
        # y voltage (min and max)
        yz_scan_y_voltage_min_widget = QLabel("y_V_min:", self) # widget
        yz_scan_y_voltage_min_widget.setParent(left_window)
        yz_scan_y_voltage_min_widget.move(yz_scan_widgets_left_x_justify, 90 + row_y_adjust + overall_y_adjust)

        yz_scan_y_voltage_min_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_y_voltage_min_qlineedit.setParent(left_window)
        yz_scan_y_voltage_min_qlineedit.move(45 + yz_scan_widgets_left_x_justify, 90 + row_y_adjust + overall_y_adjust)
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

        yz_scan_z_piezo_min_voltage_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_z_piezo_min_voltage_qlineedit.setParent(left_window)
        yz_scan_z_piezo_min_voltage_qlineedit.move(45 + yz_scan_widgets_left_x_justify, 165 + row_y_adjust + overall_y_adjust)
        yz_scan_z_piezo_min_voltage_qlineedit.resize(35, 15)
        yz_scan_z_piezo_min_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        yz_scan_z_piezo_max_voltage_widget = QLabel("z_V_max:", self) # widget
        yz_scan_z_piezo_max_voltage_widget.setParent(left_window)
        yz_scan_z_piezo_max_voltage_widget.move(yz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)

        yz_scan_z_piezo_max_voltage_qlineedit = QLineEdit(self) # qclineedit
        yz_scan_z_piezo_max_voltage_qlineedit.setParent(left_window)
        yz_scan_z_piezo_max_voltage_qlineedit.move(50 + yz_scan_widgets_left_x_justify, 190 + row_y_adjust + overall_y_adjust)
        yz_scan_z_piezo_max_voltage_qlineedit.resize(30, 15)
        yz_scan_z_piezo_max_voltage_qlineedit.setAlignment(PyQt5.QtCore.Qt.AlignRight)

        # run yz scan button
        yz_scan_run_button = QPushButton("run\nYZ scan", self) # button
        yz_scan_run_button.setParent(left_window)
        yz_scan_run_button.resize(60, 40)
        yz_scan_run_button.move(15 + yz_scan_widgets_left_x_justify, 215 + row_y_adjust + overall_y_adjust)
        yz_scan_run_button.clicked.connect(yz_scan_resolution_validation_fnc) # this framework is limited currently to only validating resolution

####################################################################### context menu ######################################################################

    def contextMenuEvent(self, event): # context (right-click) menu

        cmenu = QMenu(self)
        cmenuoneAct = cmenu.addAction("one")
        cmenutwoAct = cmenu.addAction("two")
        cmenuthreeAct = cmenu.addAction("three")
        action = cmenu.exec_(self.mapToGlobal(event.pos()))

############################################################## start gui ################################################################################

# ?
if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)
    mw = Parent()
    mw.show()
    sys.exit(app.exec_())

#################################################################### END #######################################################################################
