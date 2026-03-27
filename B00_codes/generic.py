import ConfocalT2E_Test
import ConfocalT1_Test
import ConfocalT2EOff_Test
import ConfocalT1Off_Test
import ConfocalODMRAWGFast_Test
import time
import TempSet
import numpy as np
from qcodes_contrib_drivers.drivers.Lakeshore.Lakeshore335 import Lakeshore335


def main():
    
    # ConfocalT2E_Test.main()
    # ConfocalT1_Test.main()
    
    # ConfocalT2EOff_Test.main()
    # ConfocalT1Off_Test.main()

    # # # ConfocalODMRAWGFast_Test.main(N=3e5)
    # # TempSet.main()
    # # time.sleep(600)
    # ConfocalODMRAWGFast_Test.main(N=5e5)

    for i in np.linspace(70,300,231):
        TempSet.main(temp=i)
        time.sleep(20)
    # TempSet.main(temp=9.5)
    # ls = Lakeshore335()
    # ls.allOff()
    # ls.close()

if __name__ == "__main__":
    main()