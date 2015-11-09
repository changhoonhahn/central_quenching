'''
    
Tests for star forming duty cycle

'''
import numpy as np
from scipy import signal

def test_squarewave(t): 
    ''' 
    Test whether Gaussianity is preserved with square wave 
    '''
    dSFR = sfr_A * signal.square(sfr_w * (tcosmic - 6.9048 + sfr_d))    # 
