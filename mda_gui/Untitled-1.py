
################################################################## old useful code ##########################################################################
# # combo select box 1
# combobox1 = QComboBox(self)
# combobox1.move(12, 30)
# combobox1.resize(98, 20)
# combobox1.addItem("one")
# combobox1.addItem("two")
# combobox1.addItem("three")

# # line edit 1
# qle1 = QLineEdit(self)
# qle1.move(12, 52)
# qle1.resize(395, 20)

# # creating progress bar
# pbar = QProgressBar(self)
# # pbar.resize(80, 25)
# pbar.move(160, 250)

# def scan_update(self, parent = Child):
#     self.x_data = numpy.random.randint(9, size = 5)
#     self.y_data = numpy.random.randint(9, size = 5)
#     self.sc.axes.cla()
#     self.sc.axes.plot(self.x_data, self.y_data)
#     self.sc.draw()

# def button_2_clicked(self, parent = Child):                   # use this for demo
#     # print("button 2 clicked")
#     xy_scan_function()
#     sleep(0.75)
#     scan_update()

# rightwindow.label = QLabel(rightwindow)
# rightwindow.pixmap = QPixmap("082422_img_1.png")
# rightwindow.label.setPixmap(rightwindow.pixmap)
# rightwindow.label.resize(300, 300) # other width, height: rightwindow.pixmap.width(), rightwindow.pixmap.height()

# # button1
# button1 = QPushButton("change plot (basic)", self)
# button1.resize(100, 25)
# button1.move(280, 25)
# button1.clicked.connect(button_1_clicked)

# ################### not yet #####################
# def update_plot():
#     # print("update_call called")
#     self.sc.axes.cla()
#     self.x_data = numpy.random.randint(9, size = 5)
#     self.y_data = numpy.random.randint(9, size = 5)
#     self.sc.axes.plot(self.x_data, self.y_data)
#     self.sc.draw()

# def button_1_clicked(self, parent = Child):
#     # print("button 1 clicked")
#     update_plot()

# def scan_update():                                          # use this for demo
#     self.sc.axes.cla()
#     output_data_set = numpy.load("file_name_by_time.npy")
#     self.sc.axes.pcolormesh(output_data_set, cmap = "RdBu_r")
#     self.sc.draw()
# ################### not yet ######################

# self.x_data = numpy.random.randint(9, size = 5)
# self.y_data = numpy.random.randint(9, size = 5)
# # self.z_data = numpy.random.randint(9, size = 5)
# self.sc.axes.plot(self.x_data, self.y_data)

# data_set = numpy.load("082422_faster_scan_03.npy")
# self.sc.axes.pcolormesh(data_set, cmap = "RdBu_r")
# self.sc.axes.colorbar(plot_1, ax = self.sc.axes)
# self.sc.axes.colorbar(plot_1)

# if resolution < 20:
#     make exception window show
#     hard part: wait for user fix
# else:
#     print params
#     run scan

# xy_scan_run_button.clicked.connect(xy_scan_run_button_clicked_fnc)
# xy_scan_run_button.clicked.connect(check resolution fnc)
# if check resolution fnc():
#     print_XY_scan_parameters_fnc(self)
#     run_xy_scan_fnc()
# else:

# if xy_scan_parameters_validated_Q == True:
#     xy_scan_run_button.clicked.connect(print_XY_scan_parameters_fnc) # call print_XY_scan_parameters_fnc
#     xy_scan_run_button.clicked.connect(run_xy_scan_fnc) # call xy_scan_run_button.clicked.connect(run_xy_scan_fnc)
# else:
#     xy_scan_run_button.clicked.connect(xy_scan_input_validation_fnc)

# global user_scan_parameters_validated_Q                                         # overall try to input validation
# if user_scan_parameters_validated_Q == False:
#     xy_scan_run_button.clicked.connect(scan_input_validation_fnc)                                # can this fnc return a bool checked here?
# else:
#     xy_scan_run_button.clicked.connect(print_XY_scan_parameters_fnc)
#     xy_scan_run_button.clicked.connect(run_xy_scan_fnc)


# try:                                                          # Try... Except for running mulitple scans
#     xy_scan_run_button.clicked.connect(run_xy_scan_fnc)
# except KeyError:
#     scan_galvo.close()
#     xy_scan_run_button.clicked.connect(run_xy_scan_fnc)

#
# def xy_scan_run_button_clicked_fnc(): # this fnc
#     my_bool = False
#     if resolution < 20:
#         show exception window

#     print_XY_scan_parameters_fnc(self)
#     run_xy_scan_fnc()

#
# def xy_scan_input_validation_fnc(): # this fnc...
#     print("EXCEPTION!")
#     var == True

# single option for input validation                                                             # INCOMPLETE- done locally instead
# # scanning input validation fnc
# """
# This fnc is called before any scan (XY, XZ, or YZ) is run. It will confirm valid parameters are entered by the user to prevent erros and non-allowed driving
# voltages sent to the objective z_piezo and/or the mirror motors. A list of checks is present below

# 1. Scanning resolution...
# """
# def scan_input_validation_fnc(): # this fnc validates the input for user-entered scanning (XY, XZ, and YZ) parameters
#     if int(xy_scan_resolution_qlineedit.text()) < 20:
#         print("Exception!")
#         global user_scan_parameters_validated_Q
#         user_scan_parameters_validated_Q == True
#     # print("resolution is: %d" % int(xy_scan_resolution_qlineedit.text()))

from qcodes_contrib_drivers.drivers.NationalInstruments.class_file import *

# naming the instrument
scan_galvo_card_name = "cDAQ1Mod2"

# dictionary of analog output channels
scan_galvo_ao_channels = {f'{scan_galvo_card_name}/ao{i}': i for i in range(4)}

# defining the instrument (ni_9263)
scan_galvo = DAQAnalogOutputs("name_two", scan_galvo_card_name, scan_galvo_ao_channels)

scan_galvo.voltage_cdaq1mod2ao0(-0.23) # this is for the x-mirror
scan_galvo.voltage_cdaq1mod2ao1(-0.35) # this is for the y-mirror
scan_galvo.voltage_cdaq1mod2ao2(3.2) # this is for the z-piezo