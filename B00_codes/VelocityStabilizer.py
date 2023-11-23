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
    access_token = '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'
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
    
def read_file_from_dropbox(dbx):
    # Download the file content from Dropbox
    dropbox_path = '/wlm_laser1.txt'
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
    vel = Velocity(velNum=1, ifInitVpz=0, ifInitWvl=0)
    pid = PID(p=0.2,i=1,d=0,setpoint=470.47650,memory=20)
    vpz_now = 69.84

    # Create a Dropbox client
    access_token = 'sl.BqVkE-T6NM7G3vXCagIvfJURtRLsA13xdkdl5U6K9NaRrvFCN8-uRvdVlKLeHOLOLSl5mZxdAtnY25U83Rn2FLr2NLRBeo1hJkcNsG_EKYaQVfRJ1EkPOlBYq6rplQ3TRHvHHBu5DZ2lpce1FZ4gqEY'
    app_key = 'jr02kdlisnsp67m'
    app_secret = 'ts6lhusxptmz4yk'
    refresh_token = '_YKA5RRymgcAAAAAAAAAAQnPHYunyWVJ6xfos_yVNImpt68T03ae7N_udN8_6FC7'
    
    dbx = dropbox.Dropbox(oauth2_access_token=access_token,
                          app_key=app_key, 
                          app_secret=app_secret, 
                          oauth2_refresh_token=refresh_token)
    timeLastRefreshedToken = time.time()

    while True:
        # print()
        data = read_file_from_dropbox(dbx)
        pid.set_pv(data)
        pid.set_cv()
        # print(pid.cv)

        vpz_correction = -pid.cv*25
        
        vel.set_vpiezo(vpz_now + vpz_correction, ifVerbose=0)
        vpz_now = vpz_now + vpz_correction
        print(vpz_now)

        if time.time() - timeLastRefreshedToken > 3600:
            dbx.refresh_access_token()
            timeLastRefreshedToken = time.time()