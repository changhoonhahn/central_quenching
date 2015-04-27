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
import cenque_groupcat as cq_group

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

def cq_evolution_groupcat_fit(i_nsnap, Mrcut=18, alpha=[0.6, 1.0], beta=[0.4, 0.9], gamma=[0.1, 0.5], delta=[0.0, 0.4]): 
    ''' Fit Evolve CenQue sSFR to SDSS group catalog sSFR. 
    Uses chi^2. Very naive calculation 

    input
    ------
    i_nsnap: final evolution snapshot #
    alpha, beta, gamma, delta: Tau parameter range to calculate from 

    '''
    # mass bins for comparison 
    panel_mass_bins = [
            [10.0, 10.5], [10.5, 11.0], [11.0, 11.5]
            ]

    # import central SSFR 
    centrals = cq_group.central_catalog(Mrcut=Mrcut) 
    central_ssfrs = [] 
    for panel_mass in panel_mass_bins:      # loop through panel mass bins 

        mass_limit = (centrals.mass >= panel_mass[0]) & (centrals.mass < panel_mass[1]) 

        central_ssfr_hist, central_ssfr_bin_edges = np.histogram(centrals.ssfr[mass_limit], 
                range=[-13.0, -7], bins=40, normed=True)
        
        # SSFR bins 
        ssfr_bin_low = central_ssfr_bin_edges[:-1]
        ssfr_bin_high = central_ssfr_bin_edges[1:]
        ssfr_bin_mid = [ 0.5*(ssfr_bin_low[i] + ssfr_bin_high[i]) 
                for i in range(len(ssfr_bin_low)) ] 

        central_ssfrs.append( [ central_ssfr_hist, ssfr_bin_mid ] )     # append to list 

    
    param_dict = {'alpha': [], 'beta': [], 'gamma': [], 'delta': [], 'least_square':[]} 
    for a in np.arange(alpha[0], alpha[1], 0.1): 
        for b in np.arange(beta[0], beta[1], 0.1): 
            for c in np.arange(gamma[0], gamma[1], 0.1): 
                for d in np.arange(delta[0], delta[1], 0.1): 
                    # import snapshot with parameters (a,b,c,d) 
                    cenq = cq.CenQue() 
                    cenq.readin(nsnap=i_nsnap, file_type='evol from 13', 
                            fq='wetzel', tau=[a, b, c, d]) 

                    least_sqr = 0.0
                    for i_panel, panel_mass in enumerate(panel_mass_bins):  
                        # loop through panel mass bins 

                        mass_limit = (cenq.mass >= panel_mass[0]) & \
                                (cenq.mass < panel_mass[1]) 

                        cenq_ssfr_hist, cenq_ssfr_bin_edges = \
                                np.histogram(cenq.ssfr[mass_limit], 
                                        range=[-13.0, -7], bins=40, normed=True)

                        # calculate least square
                        least_sqr += np.sum(((central_ssfrs[i_panel])[0] - cenq_ssfr_hist)**2)
                    
                    param_dict['alpha'].append(a) 
                    param_dict['beta'].append(b) 
                    param_dict['gamma'].append(c) 
                    param_dict['delta'].append(d) 

                    param_dict['least_square'].append(least_sqr) 

    # find minimum least square
    min_index = param_dict['least_square'].index( min(param_dict['least_square']) ) 

    print 'alpha =', (param_dict['alpha'])[min_index]
    print 'beta =', (param_dict['beta'] )[min_index]
    print 'gamma =', (param_dict['gamma'])[min_index]
    print 'delta =', (param_dict['delta'])[min_index]

if __name__=='__main__': 
    #cq_evolution_param_grid()
    cq_evolution_groupcat_fit(1)
