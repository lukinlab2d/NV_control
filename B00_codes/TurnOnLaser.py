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

import datetime as dt
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc


def turnOnLaser(channel):
    print()
    print(dt.datetime.now())

    # clock speed is in MHz - is 'status' needed in the dictionary?
    LaserParam = {'delay_time': 2, 'channel':channel}
    clock_speed = 500 # MHz
    settings = {'clock_speed': clock_speed, 'Laser': LaserParam, 'PB_type': 'USB',
                'min_pulse_dur': int(5*1e3/clock_speed)}

    pb = spc.B00PulseBlaster("SpinCorePB_toTurnOn", settings=settings)

    pb.turn_on_infinite(channel=channel, verbose=False)

    print('Laser turned on using PB channel ' + str(LaserParam['channel']))
    return pb


####################################################################################################################
if __name__ == '__main__':
    pb = turnOnLaser(channel=3)

    

    