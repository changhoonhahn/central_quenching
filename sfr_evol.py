'''

Star formation rate evolution for CenQue and Lineage class objects

'''
import time
import numpy as np
from scipy import signal

from util.cenque_utility import get_zsnap
from sfms.fitting import get_param_sfr_mstar_z

def get_sfrevol_param(ngal, indices, **sfrevol_prop): 
    ''' 
    Generate SFR evolution parameters used to evolve SFR based on sfrevol_prop 

    '''
    
    if sfrevol_prop['name'] == 'notperiodic':
        
        return []  

    elif sfrevol_prop['name'] == 'squarewave': 
        period = np.array([-999. for i in xrange(ngal)])
        phase = np.array([-999. for i in xrange(ngal)])

        period[indices] = np.random.uniform(sfrevol_prop['freq_range'][0], sfrevol_prop['freq_range'][1], size=len(indices))
        phase[indices] = np.random.uniform(sfrevol_prop['phase_range'][0], sfrevol_prop['phase_range'][1], size=len(indices))

        return [period, phase]

    else: 
        raise NotImplementedError("Only squarewave implemented") 

def sfr_evol(t_cosmic = None, indices = None, sfrevol_param = None, ancestor_sfr = None, ancestor_delta_sfr = None, **sfrevol_prop):
    '''
    Calculate evolved SFR excluding the SFMS evolution portion. Designed to work with 
    sf_inherit module.

    Parameters
    ----------
    t_cosmic : 
        cosmic time for how long the evolution happened
    indices : 
        indices of the ancestor_sfr, ancestor_delta_sfr, sfrevol_param parameters to evolve
    sfrevol_param : 
        parameters that specify the SFR evolution. For the case of square wave, it's the frequency and phase. 
    ancestor_sfr : 
        Initial SFR that sf_inherit works with. SFRs of the ancestor 
    ancestor_delta_sfr : 
        Initial delta SFR of the SF assignment in the ancestor
    sfrevol_prop : 
        Dictionary that specifies the star forming evolution properties

    '''
    if not isinstance(t_cosmic, float): 
        if len(t_cosmic) != len(indices):
            raise ValueError

    if not isinstance(indices, np.ndarray):
        raise ValueError('indices have to be a numpy array')
    
    if sfrevol_param is None: 
        raise ValueError

    if ancestor_sfr is None: 
        raise ValeuError

    if ancestor_delta_sfr is None: 
        raise ValeuError

    if sfrevol_prop['name'] == 'notperiodic': 

        return ancestor_sfr[indices]
    
    elif sfrevol_prop['name'] == 'squarewave': 

        freq, phase = sfrevol_param  

        #amp = ancestor_sfr[indices] / signal.square(-1.0 * freq[indices] * phase[indices])
        evolved_sfr = ancestor_sfr[indices] + \
                ancestor_delta_sfr[indices] * (signal.square(freq[indices] * (t_cosmic - phase[indices])) - 1.0)

        return evolved_sfr

def logsfr_sfduty_fluct(t0, tf, t_q=None, delta_sfr=None, sfrevol_param=None, **sfrevol_prop): 
    '''
    log(SFR) contribution from SF duty cycle fluctuation 

    log(SFR)_sfdutyfluctuation = Delta log(SFR) * squarewave( freq * (t + phase) ) 

    '''
    
    if sfrevol_prop['name'] == 'notperiodic': 

        return delta_sfr 

    elif sfrevol_prop['name'] == 'squarewave': 

        start_time = time.time()
        freq, phase = sfrevol_param 
        #print 'sfduty ', time.time()-start_time

        start_time = time.time()
        if t_q is None: 
            t_cosmic = tf - t0
        else: 
            t_cosmic = t_q - t0
            notqing = np.where(tf <= t_q)
            t_cosmic[notqing] = tf - t0
        #print 'sfduty ', time.time()-start_time

        start_time = time.time()
        sfduty = delta_sfr * signal.square(freq * (t_cosmic - phase))
        #print 'sfduty ', time.time()-start_time

        return sfduty

def logsfr_sfms_evol_t(t0, tf): 
    '''
    log(SFR) from the SFMS evolution 

    log(SFR)_SFMS = -0.76 * (z(t0) - z(tf))

    Takes long
    '''
    return 0.76 * (get_zsnap(tf) - get_zsnap(t0))

def logsfr_sfms_evol(z0, zf, z_q = None): 
    '''
    log(SFR) from the SFMS evolution 

    log(SFR)_SFMS = -0.76 * (z0 - zf)

    '''

    if z_q is None: 
        sfms = 0.76 * (zf - z0)
    else:
        sfms = 0.76 * (z_q - z0)

        notqing = np.where(z_q <= zf)
        sfms[notqing] = 0.76 * (zf - z0)
    
    return sfms 

def logsfr_quenching(tq, tf, tau=None): 
    '''
    log(SFR) contribution from quenching  
                                
    log(SFR)_quenching = np.log10( np.exp( -(t_f - t_Q)/tau) )
    '''
    qing = np.where(tf >= tq)
    logsfrq = np.zeros(len(tq))

    logsfrq[qing] = np.log10( np.exp( (tq[qing] - tf) / tau[qing] ) ) 

    if np.max(logsfrq) > 0.0: 
        i_max = np.argmax(logsfrq)
        print logsfrq[i_max], 'tq', tq[i_max], 'tf', tf, 'tau', tau[i_max]
        raise ValueError

    return logsfrq
