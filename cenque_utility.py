'''

Utility functions for CenQue

Author(s): ChangHoon Hahn

'''

import numpy as np

# quiescent fraction ------------------------------------------------------------
class fq: 
    ''' quiescent fraction 
    '''
    def __init__(self): 
        self.mass = None
        self.fq = NoneA

    def import_fq(self, lit='cosmos'): 
        ''' Import quiescent fraction from literature
        '''




# mass bins ------------------------------------------------------------------------ 
class mass_bin:         
    ''' Mass bin class  
    '''
    def __init__(self): 
        self.mass_mid = None 
        self.mass_low = None 
        self.mass_high = None 
        self.nbins = None

    def build_delta_mass(delta_mass, mass_min=9.5, mass_max=11.5): 
        ''' build mass bins given min, max mass and delta 
        '''
        pass

def simple_mass_bin(): 
    ''' Simple mass bins 
    '''
    simple_mass_bin = mass_bin()
    simple_mass_bin.mass_low = [ 9.5, 10.0, 10.5, 11.0]
    simple_mass_bin.mass_high = [ 10.0, 10.5, 11.0, 11.5]
    simple_mass_bin.mass_mid = [
            0.5 * (simple_mass_bin.mass_low[i] + simple_mass_bin.mass_high[i]) 
            for i in range(len(simple_mass_bin.mass_low))]
    simple_mass_bin.nbins = len(simple_mass_bin.mass_low)

    return simple_mass_bin 
