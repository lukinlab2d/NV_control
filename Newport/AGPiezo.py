##################################################################################
## This Python sample code demonstrates how to write code in the .NET development 
## environment that interacts with the Agilis Command Library (CmdLibAgilis.dll).
##################################################################################
import sys
import os
import inspect
# Import the .NET Common Language Runtime (CLR) to allow interaction with .NET
import clr
import numpy as np
import time


# Add the DLL folder path to the system search path (before adding references)
agpiezo_directory = 'C:/Users/lukin2dmaterials/miniconda3/envs/NV_control/Newport/Piezo Motion Control/AG-UC2-UC8/Samples/Python/StepAmplitude'
sys.path.append (agpiezo_directory)

# Add a reference to each .NET assembly required
clr.AddReference ("CmdLibAgilis")
clr.AddReference ("VCPIOLib")

# Import a class from a namespace
from Newport.Motion.CmdLibAgilis import *
from Newport.VCPIOLib import *
from System.Text import StringBuilder

class AGPiezo():
    def __init__(self, COMchannel='COM7') -> None:
       
        print ("Waiting for device discovery...\n")

        # Call the Virtual COM Port I/O Library constructor with 
        # true passed in so that logging is turned on for this sample
        self.oDeviceIO = VCPIOLib (True)
        self.oCmdLib = CmdLibAgilis (self.oDeviceIO)

        # Discover the devices that are available for communication
        self.oDeviceIO.DiscoverDevices ()
            
        # Get the list of discovered devices
        strDeviceKeyList = np.array ([])
        strDeviceKeyList = self.oDeviceIO.GetDeviceKeys ()
        for oDeviceKey in strDeviceKeyList :
            strDeviceKey = str (oDeviceKey)
            # print(strDeviceKey)
        self.COMchannel = COMchannel
        self.strDeviceKeyList = strDeviceKeyList

        self.strDeviceKey = self.OpenDeviceWithCOM(self.COMchannel) 
        print ("Device Key = %s" % self.strDeviceKey)

        if (not strDeviceKeyList) :
            print ("No devices discovered.\n")

        self.oCmdLib.SetRemoteMode()


    # <summary>
    # This method opens the first valid device in the list of discovered devices.
    # </summary>
    def OpenFirstValidDevice (self) :
        # For each device key in the list
        for oDeviceKey in self.strDeviceKeyList :
            strDeviceKey = str (oDeviceKey)
            
            # If the device was opened
            if (self.oCmdLib.Open (strDeviceKey) == 0) :
                return strDeviceKey

        # No device was opened
        return ""

    def OpenDeviceWithCOM(self,s) :
        # For each device key in the list
        for oDeviceKey in self.strDeviceKeyList :
            strDeviceKey = str (oDeviceKey)

            if strDeviceKey == s:
                self.oCmdLib.Open (strDeviceKey)
                return strDeviceKey
        # No device was opened
        return ""

    # <summary>
    # This method gets the Step Amplitude in the negative direction from the specified device.
    # </summary>
    # <returns>True for success, false for failure.</returns>
    def GetStepAmplitudeNegative (self, ax) :
        global nStepAmplitudeNeg
        # Get the controller's Step Amplitude in the negative direction
        bStatus, nStepAmplitudeNeg = self.oCmdLib.GetStepAmplitudeNegative (ax, nStepAmplitudeNeg)
        
        if (bStatus) :
            return True

        print ("Could not get the Step Amplitude in the negative direction.\n")
        return False
    
    def GetStepAmplitudePositive (self, ax) :
        global nStepAmplitudePos
        # Get the controller's Step Amplitude in the positive direction
        bStatus, nStepAmplitudePos = self.oCmdLib.GetStepAmplitudePositive (ax, nStepAmplitudePos)
        
        if (bStatus) :
            return True

        print ("Could not get the Step Amplitude in the positive direction.\n")
        return False

    # <summary>
    # This method sets the Step Amplitude in the negative direction from the specified device.
    # </summary>
    # <returns>True for success, false for failure.</returns>
    def SetStepAmplitudeNegative (self, ax, n) :
        # Set the controller's Step Amplitude in the negative direction
        bStatus = self.oCmdLib.SetStepAmplitudeNegative (ax, n)
        
        if (bStatus) :
            return True

        print ("Could not set the Step Amplitude in the negative direction.\n")
        return False

    def SetStepAmplitudePositive (self, ax, n) :
        # Set the controller's Step Amplitude in the positive direction
        bStatus = self.oCmdLib.SetStepAmplitudePositive (ax, n)
        
        if (bStatus) :
            return True

        print ("Could not set the Step Amplitude in the positive direction.\n")
        return False

    def AbsoluteMove (self, ax, tar) :
        # This method performs an absolute move to the specified target position. Tar is from 0 to 1000
        bStatus = self.oCmdLib.AbsoluteMove(ax, tar)
        
        if (bStatus) :
            return True

        print ("Could not do an absolute move.\n")
        return False
    
    def RelativeMove(self, ax, step):
        # This method performs a relative move. Tar is from 0 to 1000
        bStatus = self.oCmdLib.RelativeMove(ax, step)

        if (bStatus) :
            return True

        print ("Could not do a relative move.\n")
        return False
    
    def Close(self):
        self.oCmdLib.Close ()
        print ("Shutting down.")

        # Shut down all communication
        self.oDeviceIO.Shutdown ()

    def Wait(self,t):
        time.sleep(t)











# # <summary>
# # This is the main method.  It shows how to communicate 
# # with multiple devices in the same program.
# # </summary>
# print ("Waiting for device discovery...\n")

# # Call the Virtual COM Port I/O Library constructor with 
# # true passed in so that logging is turned on for this sample
# oDeviceIO = VCPIOLib (True)
# oCmdLib = CmdLibAgilis (oDeviceIO)

# # Discover the devices that are available for communication
# oDeviceIO.DiscoverDevices ()
    
# # Get the list of discovered devices
# strDeviceKeyList = np.array ([])
# strDeviceKeyList = oDeviceIO.GetDeviceKeys ()
# for oDeviceKey in strDeviceKeyList :
#     strDeviceKey = str (oDeviceKey)
#     print(strDeviceKey)

# # If no devices were discovered
# if (not strDeviceKeyList) :
#     print ("No devices discovered.\n")
# else :
# # Open the first valid device in the list of discovered devices
# # strDeviceKey = OpenFirstValidDevice ()
# strDeviceKey = OpenDeviceWithCOM("COM7") 
# print ("Device Key = %s" % strDeviceKey)

# # If the device was opened
# if (strDeviceKey != "") :
#     # Set the controller to Remote Mode
#     if (oCmdLib.SetRemoteMode ()) :
#         # Initialize variables
#         nAxis = 1
#         nChannel = 1
#         knStepAmplitudeMax = 50
#         nStepAmplitudeNeg = 0
#         nStepAmplitudePos = 0
        
#         # Set the current channel
#         if (oCmdLib.SetChannel (nChannel)) :
#             # Get the controller's Step Amplitude in the negative direction
#             if (GetStepAmplitudeNegative ()) :
#                 print ("Negative Step Amplitude = %d" % nStepAmplitudeNeg)

#                 # Get the controller's Step Amplitude in the positive direction
#                 if (GetStepAmplitudePositive ()) :
#                     print ("Positive Step Amplitude = %d" % nStepAmplitudePos)
                    
#                     # Update the Step Amplitude in the positive and negative direction
#                     nStepAmplitudeNeg = knStepAmplitudeMax
#                     nStepAmplitudePos = knStepAmplitudeMax
#                     print ("Updating the Negative Step Amplitude to: %d" % nStepAmplitudeNeg)
#                     print ("Updating the Positive Step Amplitude to: %d\n" % nStepAmplitudePos)
                    
#                     # Set the controller's Step Amplitude in the negative direction
#                     if (SetStepAmplitudeNegative ()) :
#                         # Set the controller's Step Amplitude in the positive direction
#                         if (SetStepAmplitudePositive ()) :
#                             # Get the controller's Step Amplitude in the negative direction
#                             if (GetStepAmplitudeNegative ()) :
#                                 print ("New Negative Step Amplitude = %d" % nStepAmplitudeNeg)

#                                 # Get the controller's Step Amplitude in the positive direction
#                                 if (GetStepAmplitudePositive ()) :
#                                     print ("New Positive Step Amplitude = %d" % nStepAmplitudePos)

#                             # Set the controller's Step Amplitude in the positive direction
#                             nStepAmplitudePos = 35
#                             SetStepAmplitudePositive ()

#                         # Set the controller's Step Amplitude in the negative direction
#                         nStepAmplitudeNeg = 35
#                         SetStepAmplitudeNegative ()
#         else :
#             print ("Could not set the current channel.\n")
#     else :
#         print ("Could not put the controller into Remote Mode.\n")

#     # Close the device
#     oCmdLib.Close ()
# else :
#     print ("Could not open the device.\n")

# print ("Shutting down.")

# # Shut down all communication
# oDeviceIO.Shutdown ()
