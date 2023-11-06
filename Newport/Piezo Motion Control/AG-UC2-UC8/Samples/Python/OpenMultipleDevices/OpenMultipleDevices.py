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

print ("Python %s\n\n" % (sys.version,))

strCurrFile = os.path.abspath (inspect.stack()[0][1])
print ("Executing File = %s\n" % strCurrFile)

# Initialize the DLL folder path to where the DLLs are located
strPathDllFolder = os.path.dirname (strCurrFile)
print ("Executing Dir  = %s\n" % strPathDllFolder)

# Add the DLL folder path to the system search path (before adding references)
sys.path.append (strPathDllFolder)

# Add a reference to each .NET assembly required
clr.AddReference ("CmdLibAgilis")
clr.AddReference ("VCPIOLib")

# Import a class from a namespace
from Newport.Motion.CmdLibAgilis import *
from Newport.VCPIOLib import *
from System.Text import StringBuilder

# <summary>
# This is the main method.  It shows how to communicate 
# with multiple devices in the same program.
# </summary>
print ("Waiting for device discovery...\n")

# Call the Virtual COM Port I/O Library constructor with 
# true passed in so that logging is turned on for this sample
oDeviceIO = VCPIOLib (True)
oCmdLib = CmdLibAgilis (oDeviceIO)

# Discover the devices that are available for communication
oDeviceIO.DiscoverDevices ()
    
# Get the list of discovered devices
strDeviceKeyList = np.array ([])
strDeviceKeyList = oDeviceIO.GetDeviceKeys ()
n = -1

# If no devices were discovered
if (not strDeviceKeyList) :
    print ("No devices discovered.\n")
else :
    # For each device key in the list
    for oDeviceKey in strDeviceKeyList :
        strDeviceKey = str (oDeviceKey)
        n = n + 1
        strOut = "Device Key[{}] = {}"
        print (strOut.format (n, strDeviceKey))
        oCmdLib.WriteLog (strOut.format (n, strDeviceKey))
        
        # If the device was opened
        if (oCmdLib.Open (strDeviceKey) == 0) :
            bStatus = False
            strFirmwareVersion = ""
            bStatus, strFirmwareVersion = oCmdLib.GetFirmwareVersion (strFirmwareVersion)
    
            # If the firmware version was read
            if (bStatus) :
                strOut = "Device ID[{}] = '{}'\n"
                print (strOut.format (n, strFirmwareVersion))
                oCmdLib.WriteLog (strOut.format (n, strFirmwareVersion))
            else :
                print ("Could not get the firmware version.\n")

            # Close the device
            oCmdLib.Close ()

print ("Shutting down.")

# Shut down all communication
oDeviceIO.Shutdown ()
