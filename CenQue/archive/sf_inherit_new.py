""" 

Evolve star forming properties of ancestor CenQue objects within the 
Lineage object ancestor for the descendant CenQue objects. New method 

"""
import time
import numpy as np
from scipy import interpolate

from lineage import Lineage

import sfr_evol

from gal_prop import Fq
from util.cenque_utility import intersection_index

import matplotlib.pyplot as plt
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors


# spline between z and t_cosmic
z, t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 

def InheritSF(nsnap_descendant, nsnap_ancestor=20, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = {
            'fq': {'name': 'wetzelsmooth'}, 
            'sfms': {'name': 'linear', 'mslope': 0.55, 'zslope': 1.1}},
        evol_prop = {
            'pq': {'slope': 0.05, 'yint': 0.0}, 
            'tau': {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
            'sfr': {'dutycycle': {'name': 'notperiodic'}}, 
            'mass': {'name': 'sham'}}):
    ''' Evolve star formation properties of ancestor CentralGalaxyPopulation class within 
    the Lineage Class to descendant CentralGalaxlyPopulation object

    Parameters
    ----------
    nsnap_ancestor : int
        Snapshot number of ancestor CGPop object. The ancestor object is constructed 
        with from subhalo catalog with subhalo_prop properties. They are also assigned
        SFRs with sfr_prop properties. 
    subhalo_prop : dict
        Dictionary that describes the subhalo properties. The key 'scatter' corresponds
        to the M*-M_halo relation. The key 'soruce' describes the source SMF used for 
        the SHAM masses.
    sfr_prop : dict
        Dictionary that describes the SFR properties assigned to the ancestor CenQue object. 
        The key 'fq' describes the quiescent fraction used for the ancestor while the key
        'sfr' describes the properties of the SFR assignment. 
    evol_prop : dict
        Dictionary that consists of dictionaries which each describe paramter choices in 
        the model. 
        - evol_prop['pq'] dictates the quenching properties. 
        - evol_prop['tau'] dictates the quenching timescale. 
        - evol_prop['sfr'] dictates the SFR evolution. 
        - evol_prop['mass'] dictates the mass evolution. 
    '''
    # make sure that snapshot = 1 is included among imported descendants
    # and the first element of the list
    if isinstance(nsnap_descendant, list): 
        raise ValueError('nsnap_descendant arg has to be an int')
    
    # evolution properties
    pq_prop = evol_prop['pq']
    tau_prop = evol_prop['tau']
    sfrevol_prop = evol_prop['sfr']
    massevol_prop = evol_prop['mass']
    
    # read in the lineage (< 0.05 seconds for one snapshot)
    read_time = time.time()
    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read(
            range(nsnap_descendant, nsnap_ancestor), 
            sfr_prop='default'
            )
    if 'subhalogrowth' in sfr_prop.keys(): 
        sfr_prop['subhalogrowth']['nsnap_descendant'] = nsnap_descendant
    bloodline.AssignSFR_ancestor(sfr_prop=sfr_prop)
    print 'Lineage Read Time = ', time.time() - read_time 

    ancestor = bloodline.ancestor    # ancestor object
    t_init = ancestor.t_cosmic
    z_init = ancestor.zsnap
    # descendant object
    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant))
    t_final = descendant.t_cosmic 
    z_final = descendant.zsnap
    print 'Evolve until z = ', z_final 

    # initialize SF properties of descendant 
    n_descendant = len(descendant.snap_index)
    descendant.sfr      = np.repeat(-999., n_descendant)
    descendant.ssfr     = np.repeat(-999., n_descendant)
    descendant.min_ssfr = np.repeat(-999., n_descendant)
    descendant.tau      = np.repeat(-999., n_descendant)
    descendant.sfr_class = np.chararray(n_descendant, itemsize=16)
    descendant.sfr_class[:] = ''
    
    succession, will = intersection_index(
            getattr(descendant, 'ancestor'+str(nsnap_ancestor)), 
            ancestor.snap_index)
    if len(succession) != len(descendant.mass): 
        raise ValueError('Something wrong with the lineage')
    q_ancestors = np.where(ancestor.sfr_class[will] == 'quiescent')[0]      # Q ancestors
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]  # SF ancestors
    
    # Evolve queiscent ancestor galaxies 
    q_time = time.time()
    descendant = _QuiescentEvol(ancestor, descendant, 
            succession=succession, will=will)
    print 'Quiescent evolution takes ', time.time()-q_time
    # Evolve Star Forming Galaxies 
    sf_time = time.time()
    _StarformingEvol(ancestor, descendant, succession=succession, will=will, 
            evol_prop=evol_prop, sfr_prop=sfr_prop, lineage=bloodline)
    print 'Star Forming Evolution takes ', time.time() - sf_time

    # Deal with over quenched galaxies 
    overquenched = np.where(
            descendant.min_ssfr[succession[sf_ancestors]] > descendant.ssfr[succession[sf_ancestors]]
            )
    if len(overquenched[0]) > 0: 
        descendant.ssfr[succession[sf_ancestors[overquenched]]] = descendant.min_ssfr[succession[sf_ancestors[overquenched]]]
        descendant.sfr[succession[sf_ancestors[overquenched]]] = descendant.ssfr[succession[sf_ancestors[overquenched]]] \
                + descendant.mass[succession[sf_ancestors[overquenched]]]
        #descendant.tau[succession[sf_ancestors[overquenched]]] = -999.

    descendant.data_columns = list(descendant.data_columns) + ['ssfr', 'sfr', 'min_ssfr', 'sfr_class']#, 'tau'])
    setattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant), descendant)

    return bloodline 

def _QuiescentEvol(ancestor, descendant, succession=None, will=None): 
    ''' Evolve quiescent ancestor galaxy population to descendant. 
    Both ancestor and descendant have to be CGPop objects. Galaxies 
    that are quiescent at z_init remain quiescent until z_final. 
    There are no mechanisms to induce star forming. These galaxies are 
    evolved by keeping the SSFR constant. M*(zf) is the SHAM M* at
    redshift zf. 

    Parameters
    ----------
    ancestor : object
        Central galaxy population object for ancestor galaxy
        population at z_init 
    descendant : object
        Central galaxly population object for descendant galaxy
        population at the final redshit of the evolution. 
    '''     
    # Q ancestors
    q_ancestors = np.where(ancestor.sfr_class[will] == 'quiescent')[0] 
    descendant.sfr_class[succession[q_ancestors]] = 'quiescent'
    descendant.ssfr[succession[q_ancestors]] = ancestor.ssfr[will[q_ancestors]]
    descendant.sfr[succession[q_ancestors]] = \
            descendant.ssfr[succession[q_ancestors]] + descendant.mass[succession[q_ancestors]]
    return descendant

def _StarformingEvol(ancestor, descendant, succession=None, will=None, sfr_prop=None, evol_prop=None, lineage=None): 
    ''' Evolve the stellar mass, star formation rates, and quench galaxies simultaneously. 
    '''
    # SFMS properties
    sfms_prop = sfr_prop['sfms']
    # fq_prop 
    fq_prop = sfr_prop['fq']
    # Quenching timescale properties
    tau_prop = evol_prop['tau']
    # SFR evolution properties 
    sfrevol_prop = evol_prop['sfr']
    # SF dutycycle properties 
    dutycycle_prop = sfrevol_prop['dutycycle']
    # mass evolution properties
    massevol_prop = evol_prop['mass']
    
    # star forming ancestors
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]  # SF ancestors

    # star-formation duty cycle parameters
    dutycycle_time = time.time()
    dutycycle_prop = sfr_evol.dutycycle_param(len(sf_ancestors), dutycycle_prop=dutycycle_prop)
    dutycycle_prop['delta_sfr'] = ancestor.delta_sfr[will[sf_ancestors]]
    print 'SFR dutycycle properties take ', time.time() - dutycycle_time, ' to generate'
    
    # keywords for logSFR(M*,t)
    kwargs_sfr = {
            'dutycycle_prop': dutycycle_prop, 
            'tau_prop': tau_prop, 
            'fq_prop': fq_prop, 
            'sfms_prop': sfms_prop, 
            'massevol_prop': massevol_prop
            }

    sham_mass =  descendant.mass[succession[sf_ancestors]].copy()   # descendant SHAM masses

    mass_time = time.time()
    descendant.mass[succession[sf_ancestors]], descendant.sfr[succession[sf_ancestors]] = \
            MstarSFR_simul_evol(
                    ancestor.mass_genesis[will[sf_ancestors]], 
                    ancestor.tsnap_genesis[will[sf_ancestors]], 
                    descendant.t_cosmic,
                    **kwargs_sfr               
                    )
    if np.min(descendant.mass[succession[sf_ancestors]] - ancestor.mass[will[sf_ancestors]]) < 0.0: 
        raise ValueError("Integrated mass can't reduce the mass")
    '''    
    elif massevol_prop['name'] == 'sham':       # SHAM masses 
        # determine the quenching stellar mass (the SHAM stellar mass at t_Q)
        M_q = _SHAM_Mquenching(kwargs_sfr['t_q'], lineage=lineage, 
                descendant=descendant, 
                descendant_index=succession[sf_ancestors])
        kwargs_sfr['M_q'] = M_q

        descendant.sfr[succession[sf_ancestors]] = \
                logSFR_M_t(sham_mass, descendant.t_cosmic, **kwargs_sfr)
    ''' 
    # calculate SSFR based on evloved SFR and M*
    descendant.ssfr[succession[sf_ancestors]] = \
            descendant.sfr[succession[sf_ancestors]] - descendant.mass[succession[sf_ancestors]]
    return     


def MstarSFR_simul_evol(M0, t0, tf, t_step=0.5, **kwargs): 
    ''' Evolve stellar mass, SFR, and quench galaxies simultaneously. 

    Notes
    -----
    * SF galaxies are quenched based on the following prescription:
        The number of galaxies that start quenching between t0 and t0+tstep 
        is determined by, 

        N_quenching = N_sf(t0) * dFq/dt(t0 + 0.5 tstep) * tstep

    '''
    dutycycle_prop = kwargs['dutycycle_prop']
    tau_prop = kwargs['tau_prop']
    fq_prop = kwargs['fq_prop']
    sfms_prop = kwargs['sfms_prop']
    massevol_prop = kwargs['massevol_prop']

    Mstar = M0 
    SFR = np.repeat(-999., len(M0))
    tQ = np.repeat(-999., len(M0))
    Mq = np.repeat(-999., len(M0))
    
    t00 = t0.min()      # earliest starting time (cosmic time of ancestor snapshot) 
    t_evol = np.arange(t00, tf+t_step, t_step) 
    t_evol[-1] = tf 

    for tt in t_evol: 
        within = np.where(t0 <= tt)    # tsnap_genesis <= t
        sf_within = np.where((t0 <= tt) & (tQ == -999.))
        Nsf_0 = len(sf_within[0])          # Number of SF galaxies  
        
        # SF galaxies that begin quenching between tt and tt + t_step 
        P_sf = np.random.uniform(low=0., high=1., size=Nsf_0)
        P_q = dFqdt(Mstar[within], tt + 0.5 * t_step, lit=fq_prop['name']) * t_step    
        q_ing = np.where(P_q > P_sf)
        # assign them quenching times
        tQ[sf_within[0][q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=len(q_ing[0]))
        
        kwargs_sfr = {
                't_init': t0[within], 
                't_q': tQ[within], 
                'Mq': Mq[within], 
                'dutycycle_prop': dutycycle_prop, 
                'tau_prop': tau_prop, 
                'sfms_prop': sfms_prop
                }
        M_evol, sfr_evol = M_integrate(Mstar[within], tt, tt+t_step, massevol_prop=massevol_prop, kwargs_sfr=kwargs_sfr)
        Mstar[within] = M_evol
        SFR[within] = sfr_evol
    
    if SFR.min() == -999.: 
        raise ValueError

    return Mstar, SFR


def logSFR_M_t(logmass, t_input, t_init=None, t_q=None, M_q=None, dutycycle_prop=None, tau_prop=None, sfms_prop=None): 
    ''' log(SFR) as a function of M* and t_cosmic.
    '''
    # SFR evolution based on solving an ODE of SFR
    logsfr_time = time.time()
    z_init = z_of_t(t_init) # initial redshift
    z_q = z_of_t(t_q)
    
    # update quenched M* ( this is the stellar mass roughly at the time of quenching)
    just_quenched = np.where(
            (t_input > t_q) & 
            (M_q == -999.))
    if len(just_quenched[0]) > 0: 
        M_q[just_quenched] = logmass[just_quenched]

    # average SFR of SFMS at M* and z_init
    quenched = np.where(t_input > t_q)
    tmp_M = logmass
    tmp_M[quenched] = M_q[quenched]
    avglogsfr = sfr_evol.AverageLogSFR_sfms(tmp_M, z_init, sfms_prop=sfms_prop)

    # log(SFR)_SFMS evolutionfrom t0 to tQ
    logsfr_sfms = sfr_evol.DeltaLogSFR_sfms(
            z_init, 
            z_of_t(t_input),
            z_q=z_q)

    # log(SFR)_duty cycle evolution from t0 to tQ
    logsfr_sfduty = sfr_evol.DeltaLogSFR_dutycycle(
            t_init, 
            t_input, 
            t_q=t_q, 
            dutycycle_prop=dutycycle_prop)

    logsfr_quench = sfr_evol.DeltaLogSFR_quenching(
            t_q, 
            t_input, 
            M_q=M_q,
            tau_prop=tau_prop)

    logsfr_tot = avglogsfr + logsfr_sfms + logsfr_sfduty + logsfr_quench
     
    return [logsfr_tot, M_q]

def M_integrate(logmass0, t0, tf, massevol_prop=None, kwargs_sfr=None):
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

    logsfr = logSFR_M_t
    
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

"""
    def _SHAM_Mquenching(tq, lineage=None, descendant=None, descendant_index=None): 
        ''' Based on the quenching time, calculate the SHAM stellar mass at the time 
        the galaxy starts to quench. Code cycles through all the descendants and 
        then obtains the SHAM mass at the snapshots closest to t_q.

        Note
        ----
        * SHAM masses sometimes decrease as a function of cosmic time...

        '''
        # get closest snapshots of galaxies that will quench
        willquench = np.where(tq < 999)[0]
        snap_q = np.repeat(-999, len(tq))
        snap_q[willquench] = get_nsnap_t(tq[willquench])
        #print tq[willquench][:10]
        #print snap_q[willquench][:10]

        M_q = np.repeat(-999., len(tq))

        nsnap_ancestor = lineage.nsnap_ancestor

        ancestor_index = getattr(descendant, 'ancestor'+str(nsnap_ancestor))[descendant_index]

        for snap in np.unique(snap_q[willquench]): 
            # galaxies that quench at this snapshot 
            quenches_here = np.where(snap_q == snap)[0]

            if snap < descendant.nsnap: 
                continue
            elif snap < nsnap_ancestor: 
                descend_snap = getattr(lineage, 'descendant_snapshot'+str(snap))
            
                f_match, snap_match = intersection_index(
                        ancestor_index[quenches_here], 
                        getattr(descend_snap, 'ancestor'+str(nsnap_ancestor))
                        )
                M_q[quenches_here[f_match]] = descend_snap.mass[snap_match] 

            elif snap == nsnap_ancestor: 
                f_match, snap_match = intersection_index(
                        ancestor_index[quenches_here], 
                        lineage.ancestor.snap_index
                        )
                M_q[quenches_here[f_match]] = lineage.ancestor.mass[snap_match]
            else: 
                raise ValueError
            #delM = (descendant.mass[descendant_index])[quenches_here[f_match]] - M_q[quenches_here[f_match]]
        #print M_q[:100]
        return M_q
"""
