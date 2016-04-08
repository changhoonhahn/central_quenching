'''

Mass evolution for CenQue and Lineage class objects

Author(s): ChangHoon Hahn 


'''
import time
import numpy as np

def integrated(logsfr, logmass0, t0, tf, massevol_prop=None, kwargs_sfr=None):
    ''' Integrated star formation rate stellar mass using Euler or RK4 integration

    M* = M*0 + f_retain * Int(SFR(t) dt, t0, tf)
    
    Parameters
    ----------
    type : 
        'euler' or 'rk4' to specify Euler or RK4 integration 
    sfr :   
        SFR function that accepts mass and t_cosmic as inputs 
    mass : 
        initial stellar mass  
    t0 : 
        initial cosmic time 
    tf : 
        final cosmic time
    f_retain : 
        fraction of stellar mass not lost from SNe and winds from Wetzel Paper
    '''
    type = massevol_prop['type']
    f_retain = massevol_prop['f_retain']
    delt = massevol_prop['t_step']

    if isinstance(t0, float): 
        niter = int(np.round( (tf-t0)/delt )) 
        delt = (tf - t0)/np.float(niter) 
    elif isinstance(t0, np.ndarray): 
        niter = int(np.round( (tf-np.min(t0))/delt )) 
        delt = (tf - t0)/np.float(niter) 
    else: 
        raise TypeError
    
    t_n_1 = t0 
    logSFR_n_1, M_q = logsfr(logmass0, t0, **kwargs_sfr)
    kwargs_sfr['M_q'] = M_q
    logM_n_1 = logmass0
    
    if niter > 0: 
        print niter, ' ', type, ' iterations'
        print 'f_reatin = ', f_retain, 'delta_t = ', delt
        for i in xrange(niter): 
            iter_time = time.time()
            t_n = t_n_1 + delt
                
            if type == 'rk4': 
                k1 = (10.0 ** logSFR_n_1)
        
                k2_sfr, tmp_Mq = logsfr(np.log10(10.0**logM_n_1 + (10**9 * delt)/2.0 * k1), t_n_1 + delt/2.0, **kwargs_sfr)
                k2 = (10.0 ** k2_sfr)

                k3_sfr, tmp_Mq = logsfr(np.log10(10.0**logM_n_1 + (10**9 * delt)/2.0 * k2), t_n_1 + delt/2.0, **kwargs_sfr)
                k3 = (10.0 ** k3_sfr)

                k4_sfr, tmp_Mq = logsfr(np.log10(10.0**logM_n_1 + (10**9 * delt) * k3), t_n_1 + delt, **kwargs_sfr)
                k4 = (10.0 ** k4_sfr)

                logM_n_1 = np.log10(10.0 ** logM_n_1 + f_retain/6.0 * (delt * 10**9) * (k1 + 2.0*k2 + 2.0*k3 + k4)) 
            elif type == 'euler': 
                #M_n = (10. ** logM_n_1) + delt * 10.**9. * f_retain * (10.** logSFR_n_1)
                #logM_n = np.log10(M_n)
                logM_n_1 = np.log10((10. ** logM_n_1) + delt * 10.**9. * f_retain * (10.** logSFR_n_1))
            else: 
                raise NotImplementedError
            
            if np.sum(np.isnan(logM_n_1)) > 0: 
                raise ValueError('There are NaNs') 
    
            # update log(SFR), and t from step n-1
            logSFR_n_1, M_q_new = logsfr(logM_n_1, t_n, **kwargs_sfr)
            t_n_1 = t_n

            kwargs_sfr['M_q'] = M_q_new      # update M_q
    
    # sanity check
    if np.min(logM_n_1 - logmass0) < 0.0: 
        if np.min(logM_n_1 - logmass0) > -0.001: 
            pass
        else: 
            raise ValueError("integrated mass cannot decrease over cosmic time")

    return logM_n_1, logSFR_n_1
