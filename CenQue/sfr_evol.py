'''

Functions for evolving the SFR in sf_inherit module 

'''
import time
import numpy as np
from scipy import signal

from util.util import get_zsnap

def AverageLogSFR_sfms(mstar, z_in, sfms_prop=None): 
    ''' Model for the average SFR of the SFMS as a function of M* at redshift z_in.
    The model takes the functional form of 

    log(SFR) = A * log M* + B * z + C

    '''
    if sfms_prop is None: 
        raise ValueError

    if sfms_prop['name'] == 'linear': 
        # mass slope
        A_highmass = 0.53
        A_lowmass = 0.53
        mslope = np.repeat(A_highmass, len(mstar))
        # z slope
        zslope = sfms_prop['zslope']            # 0.76, 1.1
        # offset 
        offset = np.repeat(-0.11, len(mstar))

    elif sfms_prop['name'] == 'kinked': # Kinked SFMS 
        # mass slope
        A_highmass = 0.53 
        A_lowmass = sfms_prop['mslope_lowmass'] 
        mslope = np.repeat(A_highmass, len(mstar))
        lowmass = np.where(mstar < 9.5)
        mslope[lowmass] = A_lowmass
        # z slope
        zslope = sfms_prop['zslope']            # 0.76, 1.1
        # offset
        offset = np.repeat(-0.11, len(mstar))
        offset[lowmass] += A_lowmass - A_highmass 

    mu_SFR = mslope * (mstar - 10.5) + zslope * (z_in-0.0502) + offset
    return mu_SFR

def ScatterLogSFR_sfms(mstar, z_in, sfms_prop=None): 
    ''' Scatter of the SFMS logSFR as a function of M* and 
    redshift z_in. Hardcoded at 0.3 
    '''
    if sfms_prop is None: 
        raise ValueError
    return 0.3 

def AverageLogSSFR_q_peak(mstar):  
    ''' Average log(SSFR) of the quiescent peak of the SSFR distribution 
    '''
    #return -0.4 * (mstar - 11.1) - 12.61
    return 0.4 * (mstar - 10.5) - 1.73 - mstar 

def ScatterLogSSFR_q_peak(mstar):  
    ''' Scatter of the log(SSFR) quiescent peak of the SSFR distribution 
    '''
    return 0.18 

def DeltaLogSFR_SF_Q_peak(mstar, z_in, sfms_prop=None): 
    ''' Del log SFR between the average SFMS and average Q
    '''
    return AverageLogSFR_sfms(mstar, z_in, sfms_prop=sfms_prop) - (mstar + AverageLogSSFR_q_peak(mstar))

def dutycycle_param(ngal, dutycycle_prop=None): 
    ''' Get parameters for the star-formation duty cycle contribution to the 
    SFR    
    '''
    if dutycycle_prop['name'] == 'notperiodic':
        # simple  constant
        return dutycycle_prop

    elif dutycycle_prop['name'] == 'squarewave': 
        # simple square wave with fixed amplitude but 
        # uniformly sampled period and phase.
        freq = np.random.uniform(
                dutycycle_prop['freq_range'][0], 
                dutycycle_prop['freq_range'][1], 
                size=ngal)
        phase = np.random.uniform(0.0, 2.*np.pi, size=ngal)
        
        dutycycle_prop['freq'] = freq 
        dutycycle_prop['phase'] = phase 
        return dutycycle_prop 

    elif dutycycle_prop['name'] == 'newamp_squarewave': 
        # square wave with new amplitude assigned every full period  
        freq = np.random.uniform(dutycycle_prop['freq_range'][0], dutycycle_prop['freq_range'][1], size=ngal)
        phase = np.random.uniform(0.0, 2.*np.pi, size=ngal)

        # maximum of number of cycles based on frequency 
        min_period = 2.0 * np.pi / freq.max()
        # overestimated max number of cycles for convenience 
        n_cycle = int((10.0 + phase.max()) // min_period)

        amp = np.random.randn(ngal, n_cycle) * dutycycle_prop['sigma']
        #amp = (np.random.normal(
        #    loc=0., 
        #    scale=dutycycle_prop['sigma'], 
        #    size=n_cycle * ngal)
        #    ).reshape([ngal, n_cycle])
    
        dutycycle_prop['amp'] = amp 
        dutycycle_prop['freq'] = freq 
        dutycycle_prop['phase'] = phase 
        return dutycycle_prop 

    else: 
        raise NotImplementedError("Only squarewave implemented") 

def DeltaLogSFR_dutycycle(t0, tf, t_q=None, dutycycle_prop=None, indices=None): 
    ''' log(SFR) contribution from the SF duty cycle fluctuation 

    log(SFR)_sfdutyfluctuation = Delta log(SFR) * squarewave( freq * (t + phase) ) 

    '''
    if dutycycle_prop['name'] == 'notperiodic': 
        if indices is None: 
            return dutycycle_prop['delta_sfr'] 
        else: 
            return dutycycle_prop['delta_sfr'][indices]

    elif dutycycle_prop['name'] == 'squarewave': 
        freq = dutycycle_prop['freq']
        phase = dutycycle_prop['phase']

        if t_q is None: 
            t_cosmic = tf - t0
        else: 
            t_cosmic = t_q - t0
            notqing = np.where(tf <= t_q)
            t_cosmic[notqing] = tf - t0

        sfduty = dutycycle_prop['delta_sfr'] * signal.square(freq * t_cosmic - phase)

        if indicies is None: 
            return sfduty
        else: 
            return sfduty[indices]
    
    elif dutycycle_prop['name'] == 'newamp_squarewave': 
        
        if indices is not None: 
            freq = dutycycle_prop['freq'][indices]
            phase = dutycycle_prop['phase'][indices]
            amp = dutycycle_prop['amp'][indices]
        else: 
            freq = dutycycle_prop['freq']
            phase = dutycycle_prop['phase']
            amp = dutycycle_prop['amp']

        if t_q is None: 
            t_cosmic = tf - t0
        else: 
            t_cosmic = t_q - t0
            notqing = np.where(tf <= t_q)
            t_cosmic[notqing] = tf - t0
        
        # individual n_cycle 
        n_cycles = ((t_cosmic - phase) // (2.0 * np.pi / freq)).astype(int)
        #print amp.shape, np.max(n_cycles)

        n_cycles += np.arange(amp.shape[0]) * amp.shape[1]

        sfduty = amp.reshape(amp.shape[0] * amp.shape[1])[n_cycles] * \
                signal.square(freq * t_cosmic - phase)

        return sfduty

def logsfr_sfms_evol_t(t0, tf): 
    '''
    log(SFR) from the SFMS evolution 

    log(SFR)_SFMS = -0.76 * (z(t0) - z(tf))

    Takes long
    '''
    return 0.76 * (get_zsnap(tf) - get_zsnap(t0))

def DeltaLogSFR_sfms(zi, zf, z_q=None, sfms_prop=None): 
    ''' The evolution of log(SFR) in the SFMS

    Del log(SFR)_SFMS = -z_slope* (z0 - zf)

    The star formation history of SF galaxies in order to match the 
    overall SFMS evolution. Quenching galaxies still have the SFMS 
    evolution. 
    '''
    z0 = zi #np.min([zi, zf])

    zslope = sfms_prop['zslope']
    sfms = zslope * (zf - z0) 
    #if z_q is None: 
    #    sfms = zslope * (zf - z0)
    #else:
    #    sfms = zslope * (np.maximum(z_q, zf) - z0)
    #    notevol = np.where(zf >= z0)
    #    sfms[notevol] = 0.0
    return sfms 

def DeltaLogSFR_quenching(tq, tf, M_q=None, tau_prop=None): 
    ''' log(SFR) contribution from quenching  
                                
    log(SFR)_quenching = np.log10( np.exp( -(t_f - t_Q)/tau) )
    '''
    qing = np.where(tf >= tq)
    logsfrq = np.zeros(len(tq))

    if len(qing[0]) == 0: 
        return logsfrq
    else: 
        tau = getTauQ(M_q[qing], tau_prop=tau_prop)
        
        if isinstance(tf, float): 
            logsfrq[qing] = np.log10( np.exp( (tq[qing] - tf) / tau ) ) 
        else: 
            logsfrq[qing] = np.log10( np.exp( (tq[qing] - tf[qing]) / tau ) ) 

        if np.max(logsfrq) > 0.0: 
            i_max = np.argmax(logsfrq)
            print logsfrq[i_max], 'tq', tq[i_max], 'tf', tf, 'tau', tau[i_max]
            raise ValueError

        return logsfrq

def getTauQ(mstar, tau_prop={'name': 'instant'}): 
    ''' Return quenching efold based on stellar mass of galaxy, Tau(M*).
    '''
    type = tau_prop['name']

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

    elif type == 'line': 
        # param will give slope and yint of pivoted tau line 
        
        tau = tau_prop['slope'] * (mstar - tau_prop['fid_mass']) + tau_prop['yint']
        
        try: 
            if np.min(tau) < 0.001: 
                tau[np.where( tau < 0.001 )] = 0.001
        except ValueError: 
            pass 

    elif type == 'satellite':   # quenching e-fold of satellite

        tau = -0.57 * ( mstar - 9.78) + 0.8
        if np.min(tau) < 0.001:     
            tau[np.where( tau < 0.001 )] = 0.001

    elif type == 'long':      # long quenching (for qa purposes)

        n_arr = len(mstar) 
        tau = np.array([2.0 for i in xrange(n_arr)]) 

    else: 
        raise NotImplementedError('asdf')

    return tau 


"""
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
        
        elif sfrevol_prop['name'] == 'newamp_squarewave': 

            freq, phase, amp = sfrevol_param  

            #amp = ancestor_sfr[indices] / signal.square(-1.0 * freq[indices] * phase[indices])
            evolved_sfr = ancestor_sfr[indices] + \
                    ancestor_delta_sfr[indices] * (signal.square(freq[indices] * (t_cosmic - phase[indices])) - 1.0)

            return evolved_sfr
"""
