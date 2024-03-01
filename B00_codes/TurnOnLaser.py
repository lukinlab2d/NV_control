"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""

import datetime as dt
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc


def turnOnLaser(channels, instrument_name):
    print()
    print(dt.datetime.now())

    # clock speed is in MHz - is 'status' needed in the dictionary?
    LaserParam = {'delay_time': 2, 'channel':channels}
    clock_speed = 500 # MHz
    settings = {'clock_speed': clock_speed, 'Laser': LaserParam, 'PB_type': 'USB',
                'min_pulse_dur': int(5*1e3/clock_speed)}

    pb = spc.B00PulseBlaster(instrument_name, settings=settings)

    
    pb.turn_on_infinite(channels=channels, verbose=False)

    print('Laser turned on using PB channel ' + str(LaserParam['channel']))
    return pb


####################################################################################################################
if __name__ == '__main__':
    pb = turnOnLaser(channel=3)

    

    