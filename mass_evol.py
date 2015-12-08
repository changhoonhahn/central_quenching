'''

Mass evolution for CenQue and Lineage class objects

Author(s): ChangHoon Hahn 


'''
import numpy as np

def integrated_rk4(logsfr, logmass0, t0, tf, f_retain=0.6, delt=0.025, **sfrkwargs):
    ''' 
    Integrated star formation rate stellar mass using RK4 integration

    M* = M*0 + f_retain * Int(SFR(t) dt, t0, tf)
    
    Parameters
    ----------
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

    niter = int(np.round( (tf-t0)/delt )) 
    delt = (tf - t0)/np.float(niter) 
    
    t_n_1 = t0 
    logSFR_n_1 = logsfr(logmass0, t0, **sfrkwargs)
    logM_n_1 = logmass0
    
    if niter > 0: 
        #print niter, ' RK4 iterations', f_retain, delt
        for i in xrange(niter): 
            t_n = t_n_1 + delt
                
            k1 = (10.0 ** logSFR_n_1)
            k2 = (10.0 ** logsfr(np.log10(10.0**logM_n_1 + (10**9 * delt)/2.0 * k1), t_n_1 + delt/2.0, **sfrkwargs))
            k3 = (10.0 ** logsfr(np.log10(10.0**logM_n_1 + (10**9 * delt)/2.0 * k2), t_n_1 + delt/2.0, **sfrkwargs) )
            k4 = (10.0 ** logsfr(np.log10(10.0**logM_n_1 + (10**9 * delt) * k3), t_n_1 + delt, **sfrkwargs))

            logM_n = np.log10(10.0 ** logM_n_1 + f_retain/6.0 * (delt * 10**9) * (k1 + 2.0*k2 + 2.0*k3 + k4)) 
            
            if np.sum(np.isnan(logM_n)) > 0: 
                raise ValueError('There are NaNs') 

            logSFR_n_1 = logsfr(logM_n, t_n, **sfrkwargs)
            logM_n_1 = logM_n
            t_n_1 = t_n
    #print tf, t_n_1
    #print 'logsfr', logSFR_n_1
    if np.min(logM_n_1 - logmass0) < 0.0: 
        raise ValueError("integrated mass cannot decrease over cosmic time")

    return logM_n_1, logSFR_n_1

def integrated_euler(logsfr, logmass0, t0, tf, f_retain=0.6, delt=0.025, **sfrkwargs):
    ''' 
    Integrated star formation rate stellar mass using Forward euler integration
    
    M* = M*0 + f_retain * Int(SFR(t) dt, t0, tf)
    d(M*(t)) = SFR(M*(t), t)
    
    Parameters
    ----------
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

    niter = int(np.round( (tf-t0)/delt )) 
    delt = (tf - t0)/np.float(niter) 

    t_n_1 = t0 
    logSFR_n_1 = logsfr(logmass0, t0, **sfrkwargs)
    logM_n_1 = logmass0
    
    if niter > 0: 
        #print niter, ' RK4 iterations', f_retain, delt
        for i in xrange(niter): 
            t_n = t_n_1 + delt

            logM_n = np.log10(
                    (10.0 ** logM_n_1) + delt * f_retain * (10.** logsfr(np.log10(10.0 **logM_n_1), t_n_1, **sfrkwargs))
                    )
                
            if np.sum(np.isnan(logM_n)) > 0: 
                raise ValueError('There are NaNs') 

            logSFR_n_1 = logsfr(logM_n, t_n, **sfrkwargs)
            logM_n_1 = logM_n
            t_n_1 = t_n
    #print tf, t_n_1
    #print 'logsfr', logSFR_n_1
    if np.min(logM_n_1 - logmass0) < 0.0: 
        raise ValueError("integrated mass cannot decrease over cosmic time")

    return logM_n_1, logSFR_n_1
