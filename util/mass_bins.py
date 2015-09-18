"""

Class for mass bins 


"""
import numpy as np

class MassBin:         
    ''' Mass bin class  
    '''
    def __init__(self): 
        self.mass_mid = None 
        self.mass_low = None 
        self.mass_high = None 
        self.nbins = None

def simple_mass_bin(): 
    ''' Simple mass bins 

    output 
    ------
    mass_bin class that contains mass bin information
    '''
    simple_mass_binsize = 0.2
    simple_mass_bin = MassBin()
    simple_mass_bin.mass_low  = np.array([ 9.0 + np.float(i)*simple_mass_binsize for i in xrange(13) ])
    simple_mass_bin.mass_high = simple_mass_bin.mass_low + simple_mass_binsize
    simple_mass_bin.mass_mid  = 0.5 * (simple_mass_bin.mass_low + simple_mass_bin.mass_high) 

    simple_mass_bin.mass_wid = simple_mass_bin.mass_high - simple_mass_bin.mass_low 

    simple_mass_bin.nbins = len(simple_mass_bin.mass_low)

    return simple_mass_bin 
