'''

Documentation here 

Author(s): ChangHoon Hahn 

'''

import numpy as np
import random

#---- Local ----
import cenque_utility as util
import cenque as cq

if __name__=='__main__': 
    snap13 = CenQue() 
    snap13.ImportSnap(13) 
    snap13.writeout()
    snap13.AssignSFR(13)
    snap13.writeout()
