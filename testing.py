'''

Central_Quenchign project

code for testing


Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import h5py
from scipy import signal
from scipy import interpolate

#---- Local ----
import cenque_utility as util
import sf_mainseq as sfms
import bovy_plot as bovy


def sfr_squarewave_testing(mass, tcosmic, **sfr_param): 
    ''' SFR determined by square function TESTING PURPOSE

    Parameters
    ----------
    mass : stellar mass 
    tcosmic : cosmic time 
    sfr_param : SFR parameters  (amplitude, phase, and frequency) 

    Notes
    -----
    * Arbitrary sfr_param will be fed to produce random SFRs

    '''
    sfr_A = sfr_param['amp']
    sfr_d = sfr_param['phase']
    sfr_w = sfr_param['freq']

    if not isinstance(sfr_w, float): 
        if len(sfr_w) != len(mass): 
            return None

    dSFR = (10**sfr_A) * signal.square(sfr_w * (tcosmic - 6.9048 + sfr_d))    # 
     
    return np.log10(dSFR)

def sfr_constant_testing(mass, tcosmic, **sfr_param): 
    ''' Constant SFR to test integator 

    '''
    if not isinstance(mass, float):
        return np.array([0.0 for i in range(len(mass))])
    else: 
        return 0.0

def rk4_integrator_testing(): 
    ''' Test RK4 integator used for integrated stellar mass 

    '''
    ngal = 2    # 10,000 galaxies   

    # mass of galaxies
    mass = np.array([0.0, 0.0]) #np.random.uniform(9.0, 12.0, ngal)
    
    # assign random amp, freq, and phase
    sfr_amp = 0.3 * np.random.randn(ngal) 
    sfr_freq = (2.0 * np.pi)/np.random.uniform(0.01, 0.1, ngal)  
    sfr_phase = np.random.uniform(0.0, 1.0, ngal)
    
    t0 = 6.5
    for tcosmic in np.arange(7.0, 13.0, 0.5): 
        print 'SFR = ', sfr_squarewave_testing(mass, tcosmic, 
                amp = sfr_amp, freq = sfr_freq, phase = sfr_phase) 
        new_mass = util.integrated_mass_rk4(sfr_squarewave_testing, mass, 
                t0, tcosmic, f_retain = 1.0, 
                amp = sfr_amp, freq = sfr_freq, phase = sfr_phase) 
        #new_mass = util.integrated_mass_rk4(sfr_constant_testing, mass, 
        #        t0, tcosmic, f_retain = 1.0) 
        #print new_mass - mass 
        print new_mass 
        print np.log10(((tcosmic - 6.5) * 10**9))

        mass = new_mass 
        t0 = tcosmic 
        #bovy.bovy_hist(mass) 

    #plt.show() 

if __name__=='__main__': 
    rk4_integrator_testing()
