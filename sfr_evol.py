'''

Star formation rate evolution for CenQue and Lineage class objects

'''
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

def logsfr_neverquenched(indices=None, sfrevol_param=None, ancestor_cq=None, **sfrevol_prop):
    '''
    '''
    logsfr_mstar_z, sig_logsfr_mstar_z = get_param_sfr_mstar_z()

    t_ancestor = ancestor_cq.t_cosmic
    
    # below is necessary otherwise, memory gets fucked
    sfrevol_param_indices = []  
    for i_p in xrange(len(sfrevol_param_indices)):
        sfrevol_param_indices.append(sfrevol_param_indices[i_p][indices])

    def logsfr(logmass, t): 

        # avglogsfr = logsfr_mstar_z(logmass, t_ancestor)    
        # this is a simplification. logsfr should depend on the new mass at each time step 
        # but that complicates the SF duty cycle 
        avglogsfr = ancestor_cq.avg_sfr

        logsfr_sfms = logsfr_sfms_evol(t_ancestor, t)

        logsfr_sfduty = logsfr_sfduty_fluct(
                t_ancestor, 
                t, 
                delta_sfr=ancestor_cq.delta_sfr[indices], 
                sfrevol_param=sfr_evol_param_indices, 
                **sfrevol_prop
                )

        return avglogsfr + logsfr_sfms + logsfr_sfduty
   
    return logsfr


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

def logsfr_sfduty_fluct(t0, tf, delta_sfr=None, sfrevol_param=None, **sfrevol_prop): 
    '''
    log(SFR) contribution from SF duty cycle fluctuation 

    log(SFR)_sfdutyfluctuation = Delta log(SFR) * squarewave( freq * (t + phase) ) 

    '''
    
    if sfrevol_prop['name'] == 'notperiodic': 

        return 0.0

    elif sfrevol_prop['name'] == 'squarewave': 

        freq, phase = sfrevol_param 

        t_cosmic = tf - t0

        return delta_sfr * signal.square(freq * (t_cosmic - phase))

def logsfr_sfms_evol(t0, tf): 
    '''
    log(SFR) from the SFMS evolution 

    log(SFR)_SFMS = -0.76 * (z0 - zf)

    '''
    return 0.76 * (get_zsnap(tf) - get_zsnap(t0))

def logsfr_quenching(tq, tf, tau=None): 
    '''
    log(SFR) contribution from quenching  
                                
    log(SFR)_quenching = np.log10( np.exp( -(t_f - t_Q)/tau) )
    '''
    return np.log10( np.exp( (tq - tf) / tau ) ) 
