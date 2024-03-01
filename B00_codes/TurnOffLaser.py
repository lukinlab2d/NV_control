"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""

import datetime as dt
from qcodes_contrib_drivers.drivers.SpinAPI import SpinCore as spc


def turnOffLaser(pb):
    print()
    print(dt.datetime.now())

    pb.turn_off()
    pb.close()

    print('Laser turned off.')


####################################################################################################################
# if __name__ == '__main__':
#     turnOffLaser()

    

    