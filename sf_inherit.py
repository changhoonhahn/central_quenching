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
from sfms.fitting import get_param_sfr_mstar_z
from sfms.fitting import get_quiescent_mean_ssfr 
from util.cenque_utility import intersection_index
from util.tau_quenching import get_quenching_efold

import matplotlib.pyplot as plt
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors

def sf_inherit(nsnap_descendant, nsnap_ancestor = 20, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
        evol_prop = {
            'pq': {'slope': 0.05, 'yint': 0.0}, 
            'tau': {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
            'sfr': {'name': 'notperiodic'}, 
            'mass': {'name': 'sham'}}):
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
    if isinstance(nsnap_descendant, list): 
        raise ValueError('nsnap_descendant arg has to be an int')
    
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
    bloodline.readin([nsnap_descendant], subhalo_prop = subhalo_prop, sfr_prop = sfr_prop)
    print 'Lineage Read Time = ', time.time() - read_time 

    ancestor = bloodline.ancestor_cq    # ancestor CenQue object
    t_ancestor = ancestor.t_cosmic
    z_ancestor = ancestor.zsnap

    #for nsnap_descendant in nsnap_descendant: 
    descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
    t_descendant = descendant.t_cosmic 
    z_descendant = descendant.zsnap
    # initialize SF properties of descendant 
    n_descendant = len(descendant.snap_index)
    descendant.sfr    = np.repeat(-999., n_descendant)
    descendant.ssfr   = np.repeat(-999., n_descendant)
    descendant.q_ssfr = np.repeat(-999., n_descendant)
    descendant.gal_type = np.chararray(n_descendant, itemsize=16)
    descendant.gal_type[:] = ''
    descendant.tau    = np.repeat(-999., n_descendant)
    
    succession, will = intersection_index(
            getattr(descendant, 'ancestor'+str(nsnap_ancestor)), 
            ancestor.snap_index)
    if len(succession) != len(descendant.mass): 
        raise ValueError('Something wrong with the lineage')
    q_ancestors = np.where(ancestor.gal_type[will] == 'quiescent')[0]
    sf_ancestors = np.where(ancestor.gal_type[will] == 'star-forming')[0]
        
    q_time = time.time()
    # --------------------------------------------------------------------------------
    # quiescent galaxies at z_ancestor remain quiescent 
    # --------------------------------------------------------------------------------
    descendant.gal_type[succession[q_ancestors]] = 'quiescent'
    descendant.ssfr[succession[q_ancestors]] = ancestor.ssfr[will[q_ancestors]]
    #if massevol_prop['name'] != 'sham': 
    #    descendant.mass[succession[q_ancestors]] = ancestor.mass_genesis[will[q_ancestors]]
    descendant.sfr[succession[q_ancestors]] = descendant.ssfr[succession[q_ancestors]] + \
            descendant.mass[succession[q_ancestors]]

    # quenching probability (just for star-forming ancestors)
    # P_Q = ( f_Q(Mf,zf) - f_Q(M0,z0) ) / (1 - f_Q(M0,z0))
    # Mf and M0 are SHAM masses. Using SHAM probably underestimates the 
    # quenching probabilities but there is a fudge component
    P_q_offset = pq_prop['slope'] * (descendant.mass[succession[sf_ancestors]] - 9.5) + pq_prop['yint']
    lowmass = np.where(ancestor.mass_genesis[will[sf_ancestors]] < 9.5)
    fqf = get_fq(descendant.mass[succession[sf_ancestors]], descendant.zsnap, lit = sfr_prop['fq']['name'])
    fq0 = get_fq(ancestor.mass_genesis[will[sf_ancestors]], 
            ancestor.zsnap_genesis[will[sf_ancestors]], 
            lit=sfr_prop['fq']['name'])
    notallquench = np.where(fq0 < 1.0)
    #P_q = np.repeat(0.0, len(sf_ancestors))
    P_q = np.repeat(1.0, len(sf_ancestors))
    P_q[notallquench] = (fqf[notallquench] - fq0[notallquench] ) / (1.0 - fq0[notallquench]) 
    P_q += P_q_offset
    #P_q[noquench] = 0.0
    print np.min(P_q[lowmass]), np.max(P_q[lowmass]), P_q[lowmass][:10]
    print 'fqf', (fqf[lowmass])[np.argmax(fqf[lowmass])]
    print 'mass_f', descendant.mass[succession[sf_ancestors[lowmass]]][np.argmax(fqf[lowmass])]
    print 'mass_0', ancestor.mass_genesis[will[sf_ancestors[lowmass]]][np.argmax(fqf[lowmass])]
    print 'zsnap_f', descendant.zsnap 
    print 'zsnap_0', ancestor.zsnap_genesis[will[sf_ancestors[lowmass]]][np.argmax(fqf[lowmass])]
    print 'fq0', (fq0[lowmass])[np.argmax(fq0[lowmass])]
    print 'mass_0', ancestor.mass_genesis[will[sf_ancestors[lowmass]]][np.argmax(fq0[lowmass])]
    print 'mass_f', descendant.mass[succession[sf_ancestors[lowmass]]][np.argmax(fq0[lowmass])]
    print 'zsnap', ancestor.zsnap_genesis[will[sf_ancestors[lowmass]]][np.argmax(fq0[lowmass])]

    # determine which SF galaxies are quenching or not 
    # based on quenching probability 
    np.random.seed()
    randoms = np.random.uniform(0., 1., len(sf_ancestors)) 
    is_qing = np.where(P_q > randoms) 
    is_notqing = np.where(P_q <= randoms)
    #print len(is_qing[0]), ' is quenching'
    #print len(is_notqing[0]), ' is not quenching'
            
    # initialize SFR evolution parameters
    sfrevol_time = time.time()
    sfrevol_param = sfr_evol.get_sfrevol_param(len(descendant.mass), sf_ancestors, **sfrevol_prop)
    print 'sfr evolution param took ', time.time() - sfrevol_time     

    # cosmic time/redshift at which the SF galaxy is quenched 
    # is sampled uniformly between t_ancestor and t_nsnap=1
    np.random.seed()
    t_q = (t_descendant - ancestor.tsnap_genesis[will[sf_ancestors]]) * np.random.uniform(0., 1., len(sf_ancestors))
    t_q += ancestor.tsnap_genesis[will[sf_ancestors]]
    z_q = z_of_t(t_q)
    t_q[is_notqing] = 999.0     # galaxies that will never quench
    z_q[is_notqing] = -999.0
    # quenching e-fold timescale
    tau_q = get_quenching_efold(descendant.mass[succession[sf_ancestors]], tau_param=tau_prop)
    descendant.tau[succession[sf_ancestors[is_qing]]] = tau_q[is_qing]
    descendant.gal_type[succession[sf_ancestors[is_qing]]] = 'quiescent'
    descendant.gal_type[succession[sf_ancestors[is_notqing]]] = 'star-forming'
    print 'quiescent evolution profiling : ', time.time() - q_time     
    # -----------------------------------------------------------------------------------------
    # star forming evolution 
    # -----------------------------------------------------------------------------------------
    # SFR evolution parameters
    sfrevol_param_sf = []
    for i_p in xrange(len(sfrevol_param)):
        sfrevol_param_sf.append(sfrevol_param[i_p][sf_ancestors])
    
    masses, sfres = [], [] 
    def logsfr_m_t(logmass, t_input): 
        # SFR evolution based on solving an ODE of SFR
        logsfr_time = time.time()
        avglogsfr = logsfr_mstar_z(logmass, ancestor.zsnap_genesis[will[sf_ancestors]])
        # log(SFR)_SFMS evolutionfrom t0 to tQ
        logsfr_sfms = sfr_evol.logsfr_sfms_evol(
                ancestor.zsnap_genesis[will[sf_ancestors]],
                z_of_t(t_input),
                z_q = z_q)
        # log(SFR)_duty evolution from t0 to tQ
        logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                ancestor.tsnap_genesis[will[sf_ancestors]],
                t_input, 
                t_q = t_q, 
                delta_sfr=ancestor.delta_sfr[will[sf_ancestors]],
                sfrevol_param=sfrevol_param_sf, 
                **sfrevol_prop)
        logsfr_quench = sfr_evol.logsfr_quenching(
                t_q, 
                t_input, 
                tau=tau_q)
        logsfr_tot = avglogsfr + logsfr_sfms + logsfr_sfduty + logsfr_quench
        masses.append(logmass[:20])
        sfres.append(logsfr_tot[:20])
        return logsfr_tot 
         
    if massevol_prop['name'] == 'integrated': 

        sham_mass =  descendant.mass[succession[sf_ancestors]].copy()

        mass_time = time.time()
        descendant.mass[succession[sf_ancestors]], descendant.sfr[succession[sf_ancestors]] = \
                mass_evol.integrated(
                        massevol_prop['type'],
                        logsfr_m_t, 
                        ancestor.mass_genesis[will[sf_ancestors]], 
                        ancestor.tsnap_genesis[will[sf_ancestors]], 
                        t_descendant, 
                        f_retain = massevol_prop['f_retain'], 
                        delt = massevol_prop['t_step']
                        )
        if np.min(descendant.mass[succession[sf_ancestors]] - ancestor.mass[will[sf_ancestors]]) < 0.0: 
            raise ValueError("Integrated mass can't reduce the mass")
        print 'Integrated Masses takes ', time.time() - mass_time
        
        if nsnap_descendant == 1:  
            fig = plt.figure(1)
            sub = fig.add_subplot(111)
            prettyplot()
            pretty_colors = prettycolors()

            masses = np.vstack(masses)
            sfres = np.vstack(sfres)
            colors = [] 
            for i in xrange(np.shape(masses)[1]): #sf_masses.shape[1]):
                sub.plot(
                        masses[:,i], 
                        sfres[:,i],
                        color=pretty_colors[i % 20],
                        lw=2
                        )
                colors.append(pretty_colors[i % 20])
            sub.scatter(sham_mass[:np.shape(masses)[1]], sfres[-1,:], c=colors, marker='^', s=20, lw=0)
            fig.savefig('sf_inherit_galaxytracking.png', bbox_inches='tight')
            plt.close()
    elif massevol_prop['name'] == 'sham': 
        descendant.sfr[succession[sf_ancestors]] = logsfr_m_t(descendant.mass[succession[sf_ancestors]], t_descendant)

    descendant.ssfr[succession[sf_ancestors]] = \
            descendant.sfr[succession[sf_ancestors]] - descendant.mass[succession[sf_ancestors]]
    
    # final quenched SSFR. This is specified due to the fact that 
    # there's a lower bound on the SSFR. These values are effectively 
    # hardcoded in order to reproduce the quiescent peak of the SSFR 
    # distribution
    #q_ssfr_mean = get_q_ssfr_mean(descendant.mass[succession[sf_ancestors]])
    q_ssfr_mean = get_quiescent_mean_ssfr(descendant.mass[succession[sf_ancestors]])
    final_q_ssfr = 0.18 * np.random.randn(len(sf_ancestors)) + q_ssfr_mean 
    descendant.q_ssfr[sf_ancestors] = final_q_ssfr
    # ---------------------------------------------------------------------------------------------------
    # Account for over quenching  
    # ---------------------------------------------------------------------------------------------------
    overquenched = np.where(
            descendant.q_ssfr[succession[sf_ancestors]] > descendant.ssfr[succession[sf_ancestors]]
            )
    if len(overquenched[0]) > 0: 
        descendant.ssfr[succession[sf_ancestors[overquenched]]] = descendant.q_ssfr[succession[sf_ancestors[overquenched]]]
        descendant.sfr[succession[sf_ancestors[overquenched]]] = descendant.ssfr[succession[sf_ancestors[overquenched]]] \
                + descendant.mass[succession[sf_ancestors[overquenched]]]
        descendant.tau[succession[sf_ancestors[overquenched]]] = -999.

    descendant.data_columns = np.array(list(descendant.data_columns) + ['ssfr', 'sfr', 'q_ssfr', 'gal_type', 'tau'])
    setattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant), descendant)

    return bloodline 



"""
    if __name__=='__main__': 
        sf_inherit(1, 
                nsnap_ancestor = 20, 
                subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
                sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
                evol_prop = {
                    'pq': {'slope': 0.05, 'yint': 0.0}, 
                    'tau': {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
                    'sfr': {'name': 'notperiodic'}, 
                    'mass': {'name': 'sham'}
                    })
"""
