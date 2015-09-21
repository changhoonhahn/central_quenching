'''

Quenching E-fold tau


'''
import numpy as np

def get_quenching_efold(mstar, tau_param = {'name': 'instant'}): 
    ''' Return quenching efold based on stellar mass of galaxy 
    '''
    type = tau_param['name']

    if type == 'constant':      # constant tau 

        n_arr = len(mstar) 
        tau = np.array([0.5 for i in xrange(n_arr)]) 

    elif type == 'linear':      # lienar tau(mass) 

        tau = -(0.8 / 1.67) * ( mstar - 9.5) + 1.0
        #if np.min(tau) < 0.1: #    tau[ tau < 0.1 ] = 0.1
         
    elif type == 'instant':     # instant quenching 

        n_arr = len(mstar) 
        tau = np.array([0.001 for i in range(n_arr)]) 

    elif type == 'discrete': 
        # param will give 4 discrete tau at the center of mass bins 
        masses = np.array([9.75, 10.25, 10.75, 11.25]) 

        if param is None: 
            raise ValueError('asdfasdfa') 

        tau = np.interp(mstar, masses, param) 
        tau[ tau < 0.05 ] = 0.05

    elif type == 'linefit': 
        # param will give slope and yint of pivoted tau line 
        
        tau = param[0] * (mstar - 11.0) + param[1]
        #tau = param[0] * (mstar - 10.5) + param[1]  # this was the previous tau linefit (changed on 6/1/2015)
        tau[ tau < 0.1 ] = 0.1

    elif type == 'satellite':   # quenching e-fold of satellite

        tau = -0.57 * ( mstar - 9.78) + 0.8
        if np.min(tau) < 0.0:     
            tau[np.where( tau < 0.0 )] = 0.0
    else: 
        raise NotImplementedError('asdf')

    return tau 

