"""

Evolve star forming properties of ancestor CenQue objects within the 
Lineage object ancestor for the descendant CenQue objects

"""
import time
import numpy as np
from scipy import interpolate

from ssfr import Ssfr
from cenque import CenQue
from lineage import Lineage

import sfr_evol
import mass_evol

from quiescent_fraction import get_fq
from util.cenque_utility import get_zsnap
from sfms.fitting import get_param_sfr_mstar_z
from util.cenque_utility import get_q_ssfr_mean
from util.tau_quenching import get_quenching_efold

def sf_inherit(nsnap_descendants, 
        nsnap_ancestor = 20, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
        evol_prop = {
            'pq': {'slope': 0.05, 'yint': 0.0}, 
            'tau': {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
            'sfr': {'name': 'notperiodic'}, 
            'mass': {'name': 'sham'}
            }, 
        quiet = True, 
        qaplot = False,
        qaplotname = None, 
        clobber = False):
    ''' 
    Evolve star forming properties of ancestor CenQue object in Lineage Class 
    to descendant CenQue object

    Parameters
    ----------
    nsnap_ancestor : 
        Snapshot number of ancestor CenQue object attribute that constructed Lineage
    subhalo_prop : 
        Dictionary that describes the subhalo properties. The key 'scatter' corresponds
        to the M*-M_halo relation. The key 'soruce' describes the source SMF used for 
        the SHAM masses.
    sfr_prop : 
        Dictionary that describes the SFR properties assigned to the ancestor CenQue object. 
        The key 'fq' describes the quiescent fraction used for the ancestor while the key
        'sfr' describes the properties of the SFR assignment. 
    evol_prop : 
        Dictionary that consists of dictionaries which each describe paramter choices in 
        the model. 
        - evol_prop['pq'] dictates the quenching properties. 
        - evol_prop['tau'] dictates the quenching timescale. 
        - evol_prop['sfr'] dictates the SFR evolution. 
        - evol_prop['mass'] dictates the mass evolution. 

    Notes
    -----

    '''
    # make sure that snapshot = 1 is included among imported descendants
    # and the first element of the list
    if not isinstance(nsnap_descendants, list): 
        raise ValueError
    if 1 not in nsnap_descendants: 
        nsnap_descendants = [1] + nsnap_descendants
    elif nsnap_descendants[0] != 1: 
        nsnap_descendants.pop(1)
        nsnap_descendants = [1] + nsnap_descendants
    
    # spline between z and t_cosmic
    z, t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
    z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 
    
    logsfr_mstar_z, sig_logsfr_mstar_z = get_param_sfr_mstar_z()

    # evolution properties
    pq_prop = evol_prop['pq']
    tau_prop = evol_prop['tau']
    sfrevol_prop = evol_prop['sfr']
    massevol_prop = evol_prop['mass']

    # read in the lineage (< 0.05 seconds for one snapshot)
    read_time = time.time()
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.ancestor(subhalo_prop = subhalo_prop, sfr_prop = sfr_prop)
    bloodline.readin(nsnap_descendants)
    print 'Lineage Read Time = ', time.time() - read_time 

    ancestor = bloodline.ancestor_cq    # ancestor CenQue object
    t_ancestor = ancestor.t_cosmic
    z_ancestor = ancestor.zsnap
    q_ancestors = np.where(ancestor.gal_type == 'quiescent')[0]
    sf_ancestors = np.where(ancestor.gal_type == 'star-forming')[0]
        
    for nsnap_descendant in nsnap_descendants:      # Descendants 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
        t_descendant = descendant.t_cosmic 
        z_descendant = descendant.zsnap
        # initialize SF properties of descendant 
        n_descendant = len(descendant.snap_index)
        descendant.sfr    = np.repeat(-999., n_descendant)
        descendant.ssfr   = np.repeat(-999., n_descendant)
        descendant.q_ssfr = np.repeat(-999., n_descendant)
        descendant.gal_type = np.chararray(n_descendants, itemsize=16)
        descendant.gal_type[:] = ''
        # --------------------------------------------------------------------------------
        # quiescent galaxies at z_ancestor remain quiescent 
        # --------------------------------------------------------------------------------
        descendant.ssfr[q_ancestors] = ancestor.ssfr[q_ancestors]
        descendant.gal_type[q_ancestors] = 'quiescent'

        def logsfr_quiescent(logmass, t_input):
            return ancestor.ssfr[q_ancestors] + logmass
    
        if massevol_prop['name'] == 'integrated': 
            descendant.mass[q_ancestors], descendant.sfr[q_ancestors] = \
                    mass_evol.integrated(
                            massevol_prop['type'], 
                            logsfr_quiescent, 
                            ancestor.mass[q_ancestors], 
                            t_ancestor, 
                            t_descendant, 
                            f_retain=massevol_prop['f_retain'], 
                            delt = massevol_prop['t_step']
                            )
        elif massevol_prop['name'] == 'sham': 
            descendant.sfr[q_ancestors] = logsfr_quiescent(descendant.mass[q_ancestors], t_descendant)

        if nsnap_descendant == 1:       # Final Snapshot -------------------------------
            # quenching probability (just quiescent fraction)
            # P_Q = ( f_Q(zf) - f_Q(z0) ) / (1 - f_Q(z0))
            P_q_offset = pq_prop['slope'] * (descendant.mass[sf_ancestors] - 9.5) + pq_prop['yint']
            P_q = (
                    get_fq(
                        descendant.mass[sf_ancestors], 
                        descendant.zsnap, 
                        lit = ancestor.fq_prop['name']
                        ) - 
                    get_fq(
                        ancestor.mass[sf_ancestors], 
                        ancestor.zsnap, 
                        lit = ancestor.fq_prop['name']
                        )
                    ) / (1.0 -  get_fq( ancestor.mass[sf_ancestors], ancestor.zsnap, lit = ancestor.fq_prop['name']) ) + P_q_offset
            # determine which SF galaxies are quenching or not 
            # based on quenching probability 
            np.random.seed()
            randoms = np.random.uniform(0, 1, len(sf_ancestors)) 
            is_qing = np.where(P_q > randoms) 
            is_notqing = np.where(P_q <= randoms)
            
            # initialize SFR evolution parameters
            sfrevol_param = sfr_evol.get_sfrevol_param(len(descendant.mass), sf_ancestors, **sfrevol_prop)

            # cosmic time/redshift at which the SF galaxy is quenched 
            # is sampled uniformly between t_ancestor and t_nsnap=1
            np.random.seed()
            t_q = np.random.uniform(t_descendant, t_ancestor, len(sf_ancestors))
            z_q = get_zsnap(t_q)
            t_q[is_notqing] = 999.0     # galaxies that will never quench
            z_q[is_notqing] = -999.0
        
            # quenching e-fold timescale
            tau_q = get_quenching_efold(descendant.mass[sf_ancestors], tau_param=tau_prop)
        
            # final quenched SSFR. This is specified due to the fact that 
            # there's a lower bound on the SSFR. These values are effectively 
            # hardcoded in order to reproduce the quiescent peak of the SSFR 
            # distribution
            q_ssfr_mean = get_q_ssfr_mean(descendant.mass[sf_ancestors])
            final_q_ssfr = 0.18 * np.random.randn(len(sf_ancestors)) + q_ssfr_mean 
        
        # star forming decendants that will quench by Snapshot 1
        descendant.q_ssfr[sf_ancestors] = final_q_ssfr
        q_started = np.where(t_q <= t_descendant)   # SF galaxies that have started quenching
        q_notstarted = np.where(t_q > t_descendant) # SF galaxies taht have NOT started quenching
        descendant.gal_type[sf_ancestors[q_started]] = 'quiescent'
        descendant.gal_type[sf_ancestors[q_notstarted]] = 'star-forming'
        
        # ---------------------------------------------------------------------------------------------------
        # star forming evolution 
        # ---------------------------------------------------------------------------------------------------
        # SFR evolution parameters
        sfrevol_param_sf = []
        for i_p in xrange(len(sfrevol_param)):
            sfrevol_param_sf.append(sfrevol_param[i_p][sf_ancestors])

        def logsfr_m_t(logmass, t_input): 

            # SFR evolution based on solving an ODE of SFR
            avglogsfr = logsfr_mstar_z(logmass, z_ancestor)

            # log(SFR)_SFMS evolutionfrom t0 to tQ
            logsfr_sfms = sfr_evol.logsfr_sfms_evol(
                    z_ancestor, 
                    z_of_t(t_input),
                    z_q = z_q
                    )
            
            # log(SFR)_duty evolution from t0 to tQ
            logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                    t_ancestor, 
                    t_input, 
                    t_q = t_q, 
                    delta_sfr=ancestor.delta_sfr[sf_ancestors],
                    sfrevol_param=sfrevol_param_sf, 
                    **sfrevol_prop
                    )

            logsfr_quench = sfr_evol.logsfr_quenching(
                    t_q, 
                    t_input, 
                    tau=tau_q
                    )

            return avglogsfr + logsfr_sfms + logsfr_sfduty + logsfr_quench
         
        #print 'SHAM Masses : ', descendant.mass[sf_ancestors]

        if massevol_prop['name'] == 'integrated': 
            #start_time = time.time()
            if massevol_prop['type'] == 'rk4': 
                descendant.mass[sf_ancestors], descendant.sfr[sf_ancestors] = \
                        mass_evol.integrated_rk4(
                                logsfr_m_t, 
                                ancestor.mass[sf_ancestors], 
                                t_ancestor, 
                                t_descendant, 
                                f_retain = massevol_prop['f_retain'], 
                                delt = massevol_prop['t_step']
                                )
            elif massevol_prop['type'] == 'euler': 
                descendant.mass[sf_ancestors], descendant.sfr[sf_ancestors] = \
                        mass_evol.integrated_euler(
                                logsfr_m_t, 
                                ancestor.mass[sf_ancestors], 
                                t_ancestor, 
                                t_descendant, 
                                f_retain = massevol_prop['f_retain'], 
                                delt = massevol_prop['t_step']
                                )

            if np.min(descendant.mass[sf_ancestors] - ancestor.mass[sf_ancestors]) < 0.0: 
                raise ValueError("Integrated mass can't reduce the mass")
            #print 'integration time ', (time.time() - start_time)
        
        elif massevol_prop['name'] == 'sham': 
            descendant.sfr[sf_ancestors] = logsfr_m_t(
                    descendant.mass[sf_ancestors], 
                    t_descendant
                    )

        descendant.ssfr[sf_ancestors] = descendant.sfr[sf_ancestors] - descendant.mass[sf_ancestors]
        
        # ---------------------------------------------------------------------------------------------------
        # Account for over quenching  
        # ---------------------------------------------------------------------------------------------------
        overquenched = np.where(
                descendant.q_ssfr[sf_ancestors] > descendant.ssfr[sf_ancestors]
                )
        if len(overquenched[0]) > 0: 
            descendant.ssfr[sf_ancestors[overquenched]] = descendant.q_ssfr[sf_ancestors[overquenched]]
            descendant.sfr[sf_ancestors[overquenched]] = descendant.ssfr[sf_ancestors[overquenched]] \
                    + descendant.mass[sf_ancestors[overquenched]]
            descendant.gal_type[sf_ancestors[overquenched]] = 'quiescent'

        descendant.data_columns = np.array(list(descendant.data_columns) + ['ssfr', 'sfr', 'q_ssfr', 'gal_type'])
        setattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant), descendant)

    return bloodline 
