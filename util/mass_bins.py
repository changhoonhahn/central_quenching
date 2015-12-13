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
    simple_mass_bin = MassBin()
    #mb = np.arange(9.0, 12.0, 0.2)  (mass bins used for the original ABC calculation
    mb = np.arange(6.0, 12.0, 0.2)
    simple_mass_bin.mass_low  = mb[:-1]
    simple_mass_bin.mass_high = mb[1:]
    simple_mass_bin.mass_mid  = 0.5 * (simple_mass_bin.mass_low + simple_mass_bin.mass_high) 

    simple_mass_bin.mass_wid = simple_mass_bin.mass_high - simple_mass_bin.mass_low 

    simple_mass_bin.nbins = len(simple_mass_bin.mass_low)

    return simple_mass_bin 
