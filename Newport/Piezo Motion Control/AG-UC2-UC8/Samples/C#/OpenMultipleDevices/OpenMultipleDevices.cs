using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using Newport.VCPIOLib;

namespace Newport.Motion.CmdLibAgilis
{
    class OpenMultipleDevices
    {
        /// <summary>
        /// The entry point for the console application.
        /// </summary>
        /// <param name="args">The command line arguments.</param>
        static void Main (string[] args)
        {
            Console.WriteLine ("Waiting for device discovery...\r\n");
            VCPIOLib.VCPIOLib deviceIO = new VCPIOLib.VCPIOLib (true);
            CmdLibAgilis cmdLib = new CmdLibAgilis (deviceIO);

            // Discover the devices that are available for communication
            deviceIO.DiscoverDevices ();

            // Get the list of discovered devices
            string[] strDeviceKeys = deviceIO.GetDeviceKeys ();

            // If no devices were discovered
            if (strDeviceKeys.Length == 0)
            {
                Console.WriteLine ("No devices discovered.");
            }
            else
            {
                // For each device key in the list
                for (int i = 0; i < strDeviceKeys.Length; i++)
                {
                    string strDeviceKey = strDeviceKeys[i];
                    Console.WriteLine ("Device Key[{0}] = {1}", i, strDeviceKey);
                    cmdLib.WriteLog ("Device Key[{0}] = {1}", i, strDeviceKey);

                    // If the device was opened
                    if (cmdLib.Open (strDeviceKey) == 0)
                    {
                        string strFirmwareVersion = string.Empty;

                        // If the firmware version was read
                        if (cmdLib.GetFirmwareVersion (ref strFirmwareVersion))
                        {
                            Console.WriteLine ("Device ID[{0}] = '{1}'\r\n", i, strFirmwareVersion);
                            cmdLib.WriteLog ("Device ID[{0}] = '{1}'", i, strFirmwareVersion);
                        }
                        else
                        {
                            Console.WriteLine ("Could not get the firmware version.\r\n");
                        }

                        // Close the device
                        cmdLib.Close ();
                    }
                }
            }

            Console.WriteLine ("Shutting down.");

            // Shut down all communication
            deviceIO.Shutdown ();
        }
    }
}
