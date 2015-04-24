'''

Fit Central_Quenching to SDSS Group Catalog 


Author(s): ChangHoon Hahn 

'''

import numpy as np
import random
import h5py

#---- Local ----
import cenque_utility as util
import cenque as cq 

def cq_evolution_param_grid(): 
    ''' Run EvolveCenQue for different parameters
    '''
    cq.build_cenque_importsnap(fq='wetzel') 
    
    for alpha in np.arange(0.6, 1, 0.1): 
        for beta in np.arange(0.4, 0.9, 0.1): 
            for gamma in np.arange(0.1, 0.5, 0.1): 
                for delta in np.arange(0.0, 0.4, 0.1): 
                    print alpha, beta, gamma, delta
                    cq.EvolveCenQue(13, 1, fq='wetzel', tau=[alpha, beta, gamma, delta]) 

if __name__=='__main__': 
    cq_evolution_param_grid()
