"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
import numpy as np
import requests
from requests.auth import HTTPBasicAuth
from nidaqmx.constants import *
import time
import dropbox
import tracemalloc
import warnings
warnings.filterwarnings("ignore", category=ResourceWarning)

from qcodes_contrib_drivers.drivers.TLB_6700_222.Velocity import Velocity
from dropbox.exceptions import AuthError

####################################################################################################################

def read_file_on_github():
    # GitHub repository details
    owner = 'lukinlab2d'
    repo = 'NV_control'
    path = 'wlm_laser1.txt'
    folder = 'B00_codes'

    # GitHub personal access token
    access_token = 2
    timestamp = int(time.time())
    url = f'https://api.github.com/repos/{owner}/{repo}/contents/{folder}/{path}?timestamp={timestamp}'
    headers = {'Authorization': f'token {access_token}',
               'Content-Type': 'application/json',
               'Accept': 'application/json'}

    # Get existing content from GitHub, if any
    try:
        existing_data = requests.get(url, headers=headers)

        # Fetch ETag from response headers
        etag = existing_data.headers.get('ETag')
        
        # Include ETag in the next request
        headers['If-None-Match'] = etag

        # Make the request
        response = requests.get(existing_data.json()['download_url'], headers=headers)

        # Check if the content has changed
        if response.status_code == 304:
            # Content has not changed, use the previous data
            print("Content not modified.")
            return None
        elif response.status_code == 200:
            # Content has changed, process the new data
            data = response.text
            float_numbers = [float(line) for line in data.split('\n') if line]
            data = np.array(float_numbers)
            return data
        else:
            # Handle other status codes
            print(f"Request failed with status code {response.status_code}")
            return None

    except Exception as e:
        print("Error:", e)
        return None

def read_file_from_dropbox(velNum, dbx):
    # Download the file content from Dropbox
    if velNum == 1:
        dropbox_path = '/wlm_laser1.txt'
    elif velNum == 2:
        dropbox_path = '/wlm_laser2.txt'

    max_retries = 300  # Set the maximum number of retries
    retry_delay = 1 # Set the delay between retries in seconds

    for retry_count in range(max_retries):
        try:
            _, response = dbx.files_download(dropbox_path)
            file_content = response.content

            # Convert bytes to string
            data_str = file_content.decode('utf-8')

            # Split the string into lines
            lines = data_str.split('\r\n')

            # Convert each line to a float
            data_float = [float(line) for line in lines if line]
            data = np.array(data_float)

            return data

        except Exception as e:
            # Handle the exception (e.g., log the error)
            print(f"Attempt {retry_count + 1} failed: {e}")
            
            if retry_count < max_retries - 1:
                # Wait before retrying
                time.sleep(retry_delay)
                print("Retrying...")
            else:
                # Maximum retries reached, raise the exception
                raise

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
        self.error = error[-1]

        # Calculate response
        if len(error) > 1:
            self.cv = self.p * error[-1] + self.i * np.sum(error) / self.memory + self.d * (error[-1] - error[-2])
        else:
            #if only have a single error so far, then derivative is undefined so do not calculate it
            self.cv = self.p * error[-1] + self.i * np.sum(error) / self.memory

if __name__ == '__main__':
    velNum=2
    vel = Velocity(velNum=velNum, ifInitVpz=0, ifInitWvl=0)

    # Create a Dropbox client
    tkFile = 'C:/Users/lukin2dmaterials/data/tk.txt'; lines = []
    with open(tkFile, 'r') as file:
        for line in file:
            if "\n" in line: line = line[0:-1]
            lines.append(line)
    access_token = lines[0]
    app_key = lines[1]
    app_secret = lines[2]
    refresh_token = lines[3]
    
    dbx = dropbox.Dropbox(oauth2_access_token=access_token,
                          app_key=app_key, 
                          app_secret=app_secret, 
                          oauth2_refresh_token=refresh_token)
    timeLastRefreshedToken = time.time()
    ##########################

    lockStatusFile = 'C:\\Users\\lukin2dmaterials\\miniconda3\\envs\\NV_control\\B00_codes\\laser2_lockStatus.txt'
    lockStatusParam = np.loadtxt(lockStatusFile)
    setpoint = lockStatusParam[1]; vpz_old = lockStatusParam[2]

    pid = PID(p=2,i=1,d=0,setpoint=setpoint,memory=20)
    factor=25

    print('Setpoint = ' + str(setpoint) + '. Vpz guess = ' + str(vpz_old))

    while True:
        lockStatusParam = np.loadtxt(lockStatusFile)
        lockStatus = int(lockStatusParam[0]); setpoint = lockStatusParam[1]; vpz_old = lockStatusParam[2]

        if lockStatus == 1:
            data = read_file_from_dropbox(velNum=velNum,dbx=dbx)
            pid.set_parameters(setpoint=setpoint)
            pid.set_pv(data)
            pid.set_cv()

            if np.abs(pid.error*1e6) >= 3:
                vpz_correction = -pid.cv*factor
            else:
                vpz_correction = 0
            
            vpz_new = vpz_old + vpz_correction
            vel.set_vpiezo(vpz_new, ifVerbose=0)
            print(vpz_new)
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

        if time.time() - timeLastRefreshedToken > 3600:
            dbx.refresh_access_token()
            timeLastRefreshedToken = time.time()