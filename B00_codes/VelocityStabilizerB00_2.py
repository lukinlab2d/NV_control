"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
from nidaqmx.constants import *
import time
import warnings
warnings.filterwarnings("ignore", category=ResourceWarning)

from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity

####################################################################################################################

def read_file(velNum):
    # Download the file content
    if velNum == 1:
        filepath = 'C:/Users/lukin2dmaterials/pylabnet/b00_wlm/wlm_laser1.txt'
    elif velNum == 2:
        filepath = 'C:/Users/lukin2dmaterials/pylabnet/b00_wlm//wlm_laser2.txt'

    data = []
    with open(filepath, 'r') as file:
        # Loop over each line in the file
        for line in file:
            # Convert the line to a float and append it to the list
            data.append(float(line.strip()))
    data = np.array(data)
    return(data)

       
class PID:
    """Generic class for PID locking"""

    def __init__(self, p=0, i=0, d=0, setpoint=0, memory=20):
        """ Constructor for PID class

        :rtype: object
        :param p: proportional gain
        :param i: integral
        :param d: differential
        :param setpoint: setpoint for process variable
        :param memory: number of samples for integral memory
        """

        self.p = p
        self.i = i
        self.d = d
        self.memory = memory
        self.setpoint = setpoint
        self._pv = np.zeros(self.memory)
        self.cv = 0
        self.error = 0

    def set_parameters(self, p=None, i=None, d=None, setpoint=None, memory=None):
        """ Sets parameters of PID controller

        :param p: proportional gain
        :param i: integral
        :param d: differential
        :param setpoint: setpoint for process variable
        :param memory: number of samples for integral memory
        """

        if p is not None:
            self.p = p
        if i is not None:
            self.i = i
        if d is not None:
            self.d = d
        if memory is not None:
            self.memory = memory
        if setpoint is not None:
            self.setpoint = setpoint

        # Update process variable history to have appropriate length if too short
        pv_length = len(self._pv)
        if pv_length < self.memory:

            # Pad constants onto the beginning of the array
            self._pv = np.hstack((np.ones(self.memory - pv_length) * self._pv[0], self._pv))

    def set_pv(self, pv=np.zeros(10)):
        """ Sets process variable

        :param pv: process variable (measured value of process to be locked).
            This should ideally be a numpy array of recent data points
        """

        # Check the length of the input
        pv_length = len(pv)
        if pv_length!=0:
            # If it's too short, append it onto the current pv
            if pv_length < self.memory:
                self._pv = np.append(self._pv[self.memory - pv_length:], pv)

            # Otherwise just take the last elements
            else:
                self._pv = pv[pv_length - self.memory:]

    def set_cv(self):
        """Calculates the appropriate value of the control variable"""

        # Calculate error
        error = self._pv - self.setpoint
        try:
            self.error = error[-1]

            # Calculate response
            if len(error) > 1:
                self.cv = self.p * error[-1] + self.i * np.sum(error) / self.memory + self.d * (error[-1] - error[-2])
            else:
                #if only have a single error so far, then derivative is undefined so do not calculate it
                self.cv = self.p * error[-1] + self.i * np.sum(error) / self.memory
        except:
            print(error)

if __name__ == '__main__':
    velNum=2
    vel = Velocity(velNum=velNum, ifInitVpz=0, ifInitWvl=0)

    ##########################

    lockStatusFile = 'C:\\Users\\lukin2dmaterials\\miniconda3\\envs\\NV_control\\B00_codes\\laser2_lockStatus.txt'
    lockStatusParam = np.loadtxt(lockStatusFile)
    setpoint = lockStatusParam[1]; vpz_old = lockStatusParam[2]

    pid = PID(p=1,i=1,d=0,setpoint=setpoint,memory=20)
    factor=25

    print('Setpoint = ' + str(setpoint) + '. Vpz guess = ' + str(vpz_old))
    timestamp = time.time()
    while True:
        lockStatusParam = np.loadtxt(lockStatusFile)
        lockStatus = int(lockStatusParam[0]); setpoint = lockStatusParam[1]; vpz_old = lockStatusParam[2]

        if lockStatus == 1:
            data = read_file(velNum=velNum)
            pid.set_parameters(setpoint=setpoint)
            pid.set_pv(data)
            pid.set_cv()

            if np.abs(pid.error*1e6) >= 3:
                vpz_correction = -pid.cv*factor
            else:
                vpz_correction = 0
            
            vpz_new = vpz_old + vpz_correction
            vel.set_vpiezo(vpz_new, ifVerbose=0)
            if time.time()-timestamp>=1:
                print(vpz_new)
                timestamp=time.time()
            with open(lockStatusFile, 'r+') as file:
                # Read all lines from the file
                lines = file.readlines()
            
                # Seek to the beginning of the 3rd row (0-indexed)
                file.seek(len(lines[0]+lines[1])+1)

                # Write the number to the 3rd row
                file.write(str(vpz_new))
                file.truncate()

        else:
            print('Laser unlocked')
            time.sleep(0.5)
