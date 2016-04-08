'''

Fit Central_Quenching to SDSS Group Catalog 


Author(s): ChangHoon Hahn 

'''

import numpy as np
import random
import h5py
import scipy.optimize as spopt

#---- Local ----
import cenque_utility as util
import cenque as cq 
import cenque_groupcat as cq_group

def cq_evolution_param_grid(type='linefit'): 
    ''' Run EvolveCenQue for different parameters
    '''
    
    for slope in np.arange(0.5, 0.65, 0.05):    # go through slope, yint 
        for yint in np.arange(-0.05, 0.1, 0.05): 

            print 'slope = ', slope, ', y-int = ', yint 
            cq.build_cenque_importsnap(fq='wetzel', sfms_slope=slope, sfms_yint=yint) 

            print type
            if type == 'discrete': 
                # loop throuhg tau parameters
                for alpha in np.arange(0.6, 1.0, 0.1): 
                    for beta in np.arange(0.4, 0.9, 0.1): 
                        for gamma in np.arange(0.1, 0.5, 0.1): 
                            for delta in np.arange(0.1, 0.4, 0.1): 
                                print 'a = ', alpha, ' b = ', beta, ' c = ', gamma, ' d = ', delta
                                
                                cq.EvolveCenQue(13, 1, fq='wetzel', 
                                        tau='discrete', tau_param=[alpha, beta, gamma, delta], 
                                        sfms_slope=slope, sfms_yint=yint) 

            elif type == 'linefit': 
                for alpha in np.arange(0.1, 0.4, 0.1): 
                    for beta in np.arange(-0.8, -0.4, 0.1): 
                        print 'tau slope', beta, ' yint', alpha 

                        cq.EvolveCenQue(13, 1, 
                                fq='wetzel', tau='linefit', tau_param=[beta, alpha], 
                                sfms_slope=slope, sfms_yint=yint) 

            else: 
                raise NameError("asdfasdlfkjasldfkjasdjfkaldjflk") 

def cq_evolution_groupcat_fit(i_nsnap, Mrcut=18, sfms={'param':[]}, tau='linefit', tau_param={'param':[]}): 
    ''' Fit Evolve CenQue sSFR to SDSS group catalog sSFR. 
    Uses chi^2. Very naive calculation 

    Paramters
    ---------
    i_nsnap : final evolution snapshot # (default is 1 for z ~ 0.05) 
    Mrcut : absolute magntiude of the SDSS central catalog (18, 19, 20) 
    tau : 
    
    Returns 
    -------
    '''

    # mass bins for comparison 
    panel_mass_bins = [
            [10.0, 10.5], [10.5, 11.0], [11.0, 11.5]
            ]

    centrals = cq_group.central_catalog(Mrcut=Mrcut)    # read central catlaog 
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

    bestfit = {'params': [], 'least_square': []} 
    
    for sfms_par in sfms['params']: 

        slope, yint = sfms_par 

        if tau == 'discrete': 
        
            for param in tau_param['params']: 
                # import snapshot with parameters (a,b,c,d) 
                cenq = cq.CenQue() 
                cenq.readin(nsnap=i_nsnap, file_type='evol from 13', 
                        fq='wetzel', tau='discrete', tau_param=param, 
                        sfms_slope=slope, sfms_yint=yint) 

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
    
                bestfit['params'].append([slope, yint, param]) 
                bestfit['least_square'].append(least_sqr) 

        elif tau == 'linefit': 

            for param in tau_param['params']: 
                cenq = cq.CenQue() 
                cenq.readin(nsnap=i_nsnap, file_type='evol from 13', 
                        fq='wetzel', tau='linefit', tau_param=param, 
                        sfms_slope=slope, sfms_yint=yint) 
                
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
  
                bestfit['params'].append([slope, yint, param]) 
                bestfit['least_square'].append(least_sqr) 
                print [slope, yint, param]
                print least_sqr 


    # find minimum least square
    min_index = bestfit['least_square'].index( min(bestfit['least_square']) ) 

    print min(bestfit['least_square'])
    return (bestfit['params'])[min_index]

def cenq_evol_wrapper(alpha, beta, gamma, delta, sfms_slope, sfms_yint): 
    ''' Function 
    '''

    slope = sfms_slope
    yint = sfms_yint 

    cenq_file = util.cenque_file(nsnap=1, file_type='evol from 13', 
            fq='wetzel', tau=[alpha, beta, gamma, delta], 
            sfms_slope=slope, sfms_yint=yint) 

    if os.path.isfile(cenq_file) == False:      # if not yet computed then re-compute 
        sfassign_file = util.cenque_file(nsnap=13, file_type='sf assign', fq='wetzel', 
                sfms_slope=slope, sfms_yint=yint)

        if os.path.isfile(sfassign_file) == False: 
            print 'Building ',sfassign_file  
            cq.build_cenque_importsnap(fq='wetzel', sfms_slope=slope, sfms_yint=yint) 

        print 'Building ', cenq_file

        cq.EvolveCenQue(13, 1, 
                fq='wetzel', tau=[alpha, beta, gamma, delta], 
                sfms_slope=slope, sfms_yint=yint) 
    
    # readin cenque file 
    cenq = cq.CenQue() 
    cenq.readin(nsnap=1, file_type='evol from 13', 
            fq='wetzel', tau=[alpha, beta, gamma, delta], 
            sfms_slope=slope, sfms_yint=yint) 
    return cenq


def cq_evolution_groupcat_obj_func(params): 
    ''' Least Square Objective Function 
    '''
    alpha, beta, gamma, delta, slope, yint = params
    
    cenq = cenq_evol_wrapper(alpha, beta, gamma, delta, slope, yint)

    centrals = cq_group.central_catalog(Mrcut=18)    # read in SDSS group cat

    panel_mass_bins = [
            [10.0, 10.5], [10.5, 11.0], [11.0, 11.5]
            ]

    least_sqr = 0.0
    for i_panel, panel_mass in enumerate(panel_mass_bins):  
        # loop through panel mass bins 

        ct_mass_limit = (centrals.mass >= panel_mass[0]) & (centrals.mass < panel_mass[1]) 

        central_ssfr_hist, central_ssfr_bin_edges = np.histogram(centrals.ssfr[ct_mass_limit], 
                range=[-13.0, -7], bins=40, normed=True)
        
        cq_mass_limit = (cenq.mass >= panel_mass[0]) & (cenq.mass < panel_mass[1]) 

        cenq_ssfr_hist, cenq_ssfr_bin_edges = \
                np.histogram(cenq.ssfr[cq_mass_limit], 
                        range=[-13.0, -7], bins=40, normed=True)

        # calculate least square
        least_sqr += np.sum((central_ssfr_hist - cenq_ssfr_hist)**2)

    return least_sqr

def cq_evol_grpcat_optimize():#maxfun): 
    ''' Use SciPy Optimize fmin_l_bfgs_b (modified BFGS) to find best-fit parameters
    '''
    #param_guess = [ 0.8, 0.4, 0.4, 0.1, 0.56, 0.07 ] 
    param_guess = np.array([ 1.0, 0.7, 0.4, 0.1, 0.56, 0.07 ])
    #param_bounds = [ (0.0, None), (0.0, None), (0.0, None), (0.0, None), 
    #        (None, None), (None,None) ]
    return spopt.fmin_l_bfgs_b(cq_evolution_groupcat_obj_func, x0=param_guess, approx_grad=True)#, maxfun=maxfun)

if __name__=='__main__': 
    #cq_evolution_param_grid(type='linefit')
    sfms_param = {} 
    sfms_param['params'] = [] 
    for slope in np.arange(0.5, 0.65, 0.05):    # go through slope, yint 
        for yint in np.arange(-0.05, 0.1, 0.05): 
            sfms_param['params'].append([slope, yint]) 

    tau_param = {} 
    tau_param['params'] = [] 
    for alpha in np.arange(-0.8, -0.4, 0.1): 
        for beta in np.arange(0.1, 0.4, 0.1): 
            tau_param['params'].append([alpha, beta]) 

    print cq_evolution_groupcat_fit(1, sfms=sfms_param, tau='linefit', tau_param=tau_param)
    #print cq_evol_grpcat_optimize()#maxfun=300) 
