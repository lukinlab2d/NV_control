using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Newport.VCPIOLib;

namespace Newport.Motion.CmdLibAgilis
{
    class StepAmplitude
    {
        /// <summary>
        /// The maximum step amplitude value.
        /// </summary>
        private const int m_knStepAmplitudeMax = 50;

        /// <summary>
        /// The Agilis Command Library object.
        /// </summary>
        private static CmdLibAgilis m_CmdLib;
        /// <summary>
        /// The current axis.
        /// </summary>
        private static int m_nAxis = 1;
        /// <summary>
        /// The current channel.
        /// </summary>
        private static int m_nChannel = 1;
        /// <summary>
        /// The Step Amplitude in the negative direction.
        /// </summary>
        private static int m_nStepAmplitudeNeg = 0;
        /// <summary>
        /// The Step Amplitude in the positive direction.
        /// </summary>
        private static int m_nStepAmplitudePos = 0;

        /// <summary>
        /// The entry point for the console application.
        /// </summary>
        /// <param name="args">The command line arguments.</param>
        static void Main (string[] args)
        {
            Console.WriteLine ("Waiting for device discovery...\r\n");
            VCPIOLib.VCPIOLib deviceIO = new VCPIOLib.VCPIOLib (true);
            m_CmdLib = new CmdLibAgilis (deviceIO);

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
                // Open the first valid device in the list of discovered devices
                string strDeviceKey = OpenFirstValidDevice (strDeviceKeys);

                // If the device was opened
                if (strDeviceKey != null)
                {
                    // Put the controller into Remote Mode
                    if (m_CmdLib.SetRemoteMode ())
                    {
                        // Set the current channel
                        if (m_CmdLib.SetChannel (m_nChannel))
                        {
                            // Get the controller's Step Amplitude in the negative direction
                            if (GetStepAmplitudeNegative ())
                            {
                                Console.WriteLine ("Negative Step Amplitude = {0}", m_nStepAmplitudeNeg);

                                // Get the controller's Step Amplitude in the positive direction
                                if (GetStepAmplitudePositive ())
                                {
                                    Console.WriteLine ("Positive Step Amplitude = {0}\r\n", m_nStepAmplitudePos);

                                    // Update the Step Amplitude in the positive and negative direction
                                    m_nStepAmplitudeNeg = m_knStepAmplitudeMax;
                                    m_nStepAmplitudePos = m_knStepAmplitudeMax;
                                    Console.WriteLine ("Updating the Negative Step Amplitude to: {0}", m_nStepAmplitudeNeg);
                                    Console.WriteLine ("Updating the Positive Step Amplitude to: {0}\r\n", m_nStepAmplitudePos);

                                    // Set the controller's Step Amplitude in the negative direction
                                    if (SetStepAmplitudeNegative ())
                                    {
                                        // Set the controller's Step Amplitude in the positive direction
                                        if (SetStepAmplitudePositive ())
                                        {
                                            // Get the controller's Step Amplitude in the negative direction
                                            if (GetStepAmplitudeNegative ())
                                            {
                                                Console.WriteLine ("New Negative Step Amplitude = {0}", m_nStepAmplitudeNeg);

                                                // Get the controller's Step Amplitude in the positive direction
                                                if (GetStepAmplitudePositive ())
                                                {
                                                    Console.WriteLine ("New Positive Step Amplitude = {0}\r\n", m_nStepAmplitudePos);
                                                }
                                            }

                                            // Set the controller's Step Amplitude in the positive direction
                                            m_nStepAmplitudePos = 35;
                                            SetStepAmplitudePositive ();
                                        }

                                        // Set the controller's Step Amplitude in the negative direction
                                        m_nStepAmplitudeNeg = 35;
                                        SetStepAmplitudeNegative ();
                                    }
                                }
                            }
                        }
                        else
                        {
                            Console.WriteLine ("Could not set the current channel.\r\n");
                        }
                    }
                    else
                    {
                        Console.WriteLine ("Could not put the controller into Remote Mode.\r\n");
                    }

                    // Close the device
                    m_CmdLib.Close ();
                }
                else
                {
                    Console.WriteLine ("Could not open the device.");
                }
            }

            Console.WriteLine ("Shutting down.");

            // Shut down all communication
            deviceIO.Shutdown ();
        }

        /// <summary>
        /// This method opens the first valid device in the list of discovered devices.
        /// </summary>
        /// <param name="deviceKeys">The device keys of the discovered devices.</param>
        /// <returns>The device key of first device that could be opened, otherwise null.</returns>
        private static string OpenFirstValidDevice (string[] deviceKeys)
        {
            // For each device key in the list
            for (int i = 0; i < deviceKeys.Length; i++)
            {
                string strDeviceKey = deviceKeys[i];

                // If the device was opened
                if (m_CmdLib.Open (strDeviceKey) == 0)
                {
                    return strDeviceKey;
                }
            }

            // No device was opened
            return null;
        }

        /// <summary>
        /// This method gets the Step Amplitude in the negative direction from the specified device.
        /// </summary>
        /// <returns>True for success, false for failure.</returns>
        private static bool GetStepAmplitudeNegative ()
        {
            // Get the controller's Step Amplitude in the negative direction
            if (m_CmdLib.GetStepAmplitudeNegative (m_nAxis, ref m_nStepAmplitudeNeg))
            {
                return true;
            }

            Console.WriteLine ("Could not get the Step Amplitude in the negative direction.\r\n");
            return false;
        }

        /// <summary>
        /// This method sets the Step Amplitude in the negative direction for the specified device.
        /// </summary>
        /// <returns>True for success, false for failure.</returns>
        private static bool SetStepAmplitudeNegative ()
        {
            // Set the controller's Step Amplitude in the negative direction
            if (m_CmdLib.SetStepAmplitudeNegative (m_nAxis, m_nStepAmplitudeNeg))
            {
                return true;
            }

            Console.WriteLine ("Could not set the Step Amplitude in the negative direction.\r\n");
            return false;
        }

        /// <summary>
        /// This method gets the Step Amplitude in the positive direction from the specified device.
        /// </summary>
        /// <returns>True for success, false for failure.</returns>
        private static bool GetStepAmplitudePositive ()
        {
            // Get the controller's Step Amplitude in the positive direction
            if (m_CmdLib.GetStepAmplitudePositive (m_nAxis, ref m_nStepAmplitudePos))
            {
                return true;
            }

            Console.WriteLine ("Could not get the Step Amplitude in the positive direction.\r\n");
            return false;
        }

        /// <summary>
        /// This method sets the Step Amplitude in the positive direction for the specified device.
        /// </summary>
        /// <returns>True for success, false for failure.</returns>
        private static bool SetStepAmplitudePositive ()
        {
            // Set the controller's Step Amplitude in the positive direction
            if (m_CmdLib.SetStepAmplitudePositive (m_nAxis, m_nStepAmplitudePos))
            {
                return true;
            }

            Console.WriteLine ("Could not set the Step Amplitude in the positive direction.\r\n");
            return false;
        }
    }
}
