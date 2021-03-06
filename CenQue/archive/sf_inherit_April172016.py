""" 

Evolve star forming properties of ancestor CenQue objects within the 
Lineage object for the descendant CenQue objects

"""
import time
import numpy as np
from scipy import interpolate

from lineage import Lineage

import sfr_evol

from gal_prop import Fq
from gal_prop import dFqdt
import util.util as Util
from util.util import get_nsnap_t
from util.util import intersection_index

import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


# spline between z and t_cosmic
z_snap, t_snap = np.loadtxt(Util.snapshottable(), unpack=True, usecols=[2, 3]) 
z_of_t = interpolate.interp1d(list(reversed(t_snap)), list(reversed(z_snap)), kind='cubic') 


def InheritSF(nsnap_descendant, nsnap_ancestor=20, subhalo_prop=None, sfr_prop=None, evol_prop=None, quiet=True):
    ''' Evolve star formation properties of 'ancestor' CentralGalaxyPopulation class at 
    redshift specified by nsnap_ancestor to descendant CentralGalaxlyPopulation object
    at redshift specified by nsnap_descendant. Both ancestor and descendant objects
    are attributes in the Lineage class, which contains all the 'lineage' information 
    i.e. all the halo tracking information. 

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
    if isinstance(nsnap_descendant, list): 
        d_list = True 
    else: 
        d_list = False
    # read in the lineage (< 0.05 seconds for one snapshot)
    #read_time = time.time()
    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop, quiet=quiet)
    if not d_list: 
        bloodline.Read(range(nsnap_descendant, nsnap_ancestor), quiet=quiet)
    else: 
        bloodline.Read(range(np.min(nsnap_descendant), nsnap_ancestor), quiet=quiet)
    if 'subhalogrowth' in sfr_prop.keys():  # depending on whether SFR assign includes subhalo growth AM
        sfr_prop['subhalogrowth']['nsnap_descendant'] = nsnap_descendant
    bloodline.AssignSFR_ancestor(sfr_prop=sfr_prop, quiet=quiet)
    #print 'Lineage Read Time = ', time.time() - read_time 

    ancestor = bloodline.ancestor    # ancestor object
    t_init = ancestor.t_cosmic
    z_init = ancestor.zsnap
    ancestor.tQ[np.where(ancestor.tQ != 999.)] = 0.

    if not isinstance(nsnap_descendant, list): 
        nsnap_descendant = [nsnap_descendant]

    descendants = [] 
    sf_ancestors, q_ancestors = [], [] 
    successions, wills = [], [] 
    for nd in nsnap_descendant: 
        des = getattr(bloodline, 'descendant_snapshot'+str(nd))
        des._clean_initialize()   # initialize SF properties
        des.sfr_prop = ancestor.sfr_prop 
        # match indices up with each other 
        succession, will = intersection_index(
                getattr(des, 'ancestor'+str(nsnap_ancestor)), 
                ancestor.snap_index
                )
        if len(succession) != len(des.mass): 
            raise ValueError('Something wrong with the lineage')
        q_ancestor = np.where(ancestor.sfr_class[will] == 'quiescent')[0]      # Q ancestors
        sf_ancestor = np.where(ancestor.sfr_class[will] == 'star-forming')[0]  # SF ancestors
        if not quiet: 
            print "nsnap_descendant = ", nd, "Ancestors: Nq = ", len(q_ancestor), ', Nsf = ', len(sf_ancestor)

        descendants.append(des)
        wills.append(will)
        successions.append(succession)             
        q_ancestors.append(q_ancestor)
        sf_ancestors.append(sf_ancestor)

    # Evolve queiscent ancestor galaxies 
    if not quiet: 
        q_time = time.time()
    descendants = _QuiescentEvol(ancestor, descendants, successions=successions, wills=wills, q_ancestors=q_ancestors)
    if not quiet: 
        print 'Quiescent evolution takes ', time.time()-q_time

    # Evolve Star Forming Galaxies 
    if not quiet: 
        sf_time = time.time()
    ################## PQ BASED DROPPED ##############################
    # if evol_prop['type'] == 'pq_based':
    #     _StarformingEvol_Pq(ancestor, descendant, succession=succession, will=will, 
    #             evol_prop=evol_prop, sfr_prop=sfr_prop, 
    #             lineage=bloodline)
    # elif evol_prop['type'] == 'simult': 
    ################## PQ BASED DROPPED ##############################
    descendants = _StarformingEvol_SimulEvo(ancestor, descendants, successions=successions, wills=wills, 
            evol_prop=evol_prop, sfr_prop=sfr_prop, lineage=bloodline, quiet=quiet)
    if not quiet: 
        print 'Star Forming Evolution takes ', time.time() - sf_time
    
    descendants = _Overquenching(descendants, successions=successions, wills=wills, sf_ancestors=sf_ancestors)
    
    for descendant in descendants: 
        descendant.data_columns = list(descendant.data_columns) + ['ssfr', 'sfr', 'min_ssfr', 'sfr_class']
        setattr(bloodline, 'descendant_snapshot'+str(descendant.nsnap), descendant)
    return bloodline 


# Quiescent population evolution
def _QuiescentEvol(ancestor, descendants, successions=None, wills=None, q_ancestors=None): 
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
    for i_des in range(len(descendants)):
        # Q ancestors
        a_index = (wills[i_des])[q_ancestors[i_des]]
        d_index = (successions[i_des])[q_ancestors[i_des]]

        descendants[i_des].sfr_class[d_index] = 'quiescent'
        descendants[i_des].ssfr[d_index] = ancestor.ssfr[a_index]
        descendants[i_des].sfr[d_index] = \
                descendants[i_des].ssfr[d_index] + descendants[i_des].mass[d_index]
    return descendants


# SF Evolution (Simultaneous evolution) 
def _StarformingEvol_SimulEvo(ancestor, descendants, successions=None, wills=None, 
        sfr_prop=None, evol_prop=None, lineage=None, quiet=True): 
    ''' Evolve the stellar mass, star formation rates, and quench galaxies simultaneously. 
    '''
    sfms_prop = sfr_prop['sfms']        # SFMS properties
    fq_prop = sfr_prop['fq']            # fq_prop 
    tau_prop = evol_prop['tau']         # Quenching timescale properties
    fudge_prop = evol_prop['fudge']
    dutycycle_prop = evol_prop['sfr']['dutycycle']   # SF dutycycle properties 
    massevol_prop = evol_prop['mass']   # mass evolution properties
    
    if len(descendants) > 1: 
        allwill = np.array(list(set(list(np.concatenate(wills)))))   # combine will of unique elements
        t_output = [] 
        for des in descendants: 
            t_output.append(des.t_cosmic)
    else: 
        t_output = [descendants[0].t_cosmic]
        allwill = wills[0]
    sf_ancestors = np.where(ancestor.sfr_class[allwill] == 'star-forming')[0]  # SF ancestors
    q_ancestors = np.where(ancestor.sfr_class[allwill] == 'quiescent')[0]  # Q ancestors

    # star-formation duty cycle parameters
    if not quiet: 
        dutycycle_time = time.time()
    dutycycle_prop = sfr_evol.dutycycle_param(len(sf_ancestors), dutycycle_prop=dutycycle_prop)
    dutycycle_prop['delta_sfr'] = ancestor.delta_sfr[allwill[sf_ancestors]]

    if not quiet: 
        print 'SFR dutycycle properties take ', time.time() - dutycycle_time, ' to generate'
    
    MshamEvol = ancestor.Msham_evol[0]
    # keywords for logSFR(M*,t)
    kwargs_sfr = {'dutycycle_prop': dutycycle_prop, 'tau_prop': tau_prop, 'fudge_prop': fudge_prop,
            'fq_prop': fq_prop, 'sfms_prop': sfms_prop, 'massevol_prop': massevol_prop}
    if not quiet: 
        mass_time = time.time()

    # DEFUNCT FOR NOW. needs to be updated
    #if massevol_prop['name'] == 'integrated': 
    #    M_evolved, SFR_evolved, tQ = \
    #            Evol_IntegMstarSFR(
    #                    ancestor.mass_genesis[will[sf_ancestors]], 
    #                    ancestor.tsnap_genesis[will[sf_ancestors]], 
    #                    descendant.t_cosmic,
    #                    tQ0=ancestor.tQ[will[sf_ancestors]],
    #                    MQ0=ancestor.MQ[will[sf_ancestors]],
    #                    ancestor_Mq=MshamEvol[will[q_ancestors],:], 
    #                    **kwargs_sfr               
    #                    )
    #elif massevol_prop['name'] == 'sham': 
    bef_tsnap = ancestor.tsnap_genesis[allwill[sf_ancestors]]
    bef_tQ = ancestor.tQ[allwill[sf_ancestors]] 
    bef_MQ = ancestor.MQ[allwill[sf_ancestors]] 
    bef_sfr = ancestor.sfr[allwill[sf_ancestors]]
    bef_MshamEvol = MshamEvol[allwill]
    SFRs_evolved, tQ = \
            Evol_shamMstarSFR(
                    ancestor.tsnap_genesis[allwill[sf_ancestors]], 
                    t_output,
                    tQ0=ancestor.tQ[allwill[sf_ancestors]],
                    MQ0=ancestor.MQ[allwill[sf_ancestors]],
                    SFR0=ancestor.sfr[allwill[sf_ancestors]], 
                    q_Mshams=MshamEvol[allwill[q_ancestors],:], 
                    sf_Mshams=MshamEvol[allwill[sf_ancestors],:],
                    **kwargs_sfr               
                    )
    print 'tsnap stays the same', np.array_equal(bef_tsnap, ancestor.tsnap_genesis[allwill[sf_ancestors]])
    print 'tQ stays the same', np.array_equal(bef_tQ, ancestor.tQ[allwill[sf_ancestors]])
    print 'MQ stays the same', np.array_equal(bef_MQ, ancestor.MQ[allwill[sf_ancestors]]) 
    print 'SFR stays the same', np.array_equal(bef_sfr, ancestor.sfr[allwill[sf_ancestors]])
    print 'MshamEvol stays the same', np.array_equal(bef_MshamEvol, MshamEvol[allwill])

    print 'recalculating for just t=', t_output[0]
    SFRs_evolved1, tQ1 = \
            Evol_shamMstarSFR(
                    ancestor.tsnap_genesis[allwill[sf_ancestors]], 
                    [t_output[0]],
                    tQ0=ancestor.tQ[allwill[sf_ancestors]],
                    MQ0=ancestor.MQ[allwill[sf_ancestors]],
                    SFR0=ancestor.sfr[allwill[sf_ancestors]], 
                    q_Mshams=MshamEvol[allwill[q_ancestors],:], 
                    sf_Mshams=MshamEvol[allwill[sf_ancestors],:],
                    **kwargs_sfr               
                    )

    print 'tsnap stays the same', np.array_equal(bef_tsnap, ancestor.tsnap_genesis[allwill[sf_ancestors]])
    print 'tQ stays the same', np.array_equal(bef_tQ, ancestor.tQ[allwill[sf_ancestors]])
    print 'MQ stays the same', np.array_equal(bef_MQ, ancestor.MQ[allwill[sf_ancestors]]) 
    print 'SFR stays the same', np.array_equal(bef_sfr, ancestor.sfr[allwill[sf_ancestors]])
    print 'MshamEvol stays the same', np.array_equal(bef_MshamEvol, MshamEvol[allwill])

    print np.array_equal(SFRs_evolved[0], SFRs_evolved1) 
    
    dist, bins = np.histogram(SFRs_evolved[0], range=[-5, 2])
    plt.plot(0.5 * (bins[:-1] + bins[1:]), dist)
    dist, bins = np.histogram(SFRs_evolved1, range=[-5, 2])
    plt.plot(0.5 * (bins[:-1] + bins[1:]), dist)

    plt.show()

    for i_d in range(len(t_output)): 
        sf_anc = np.where(ancestor.sfr_class[wills[i_d]] == 'star-forming')[0] 
        sf_succession = (successions[i_d])[sf_anc]
        sf_will = (wills[i_d])[sf_anc]

        sf_allwill = allwill[sf_ancestors]

        ss, ww = intersection_index(sf_allwill, sf_will)
        if len(ww) != len(sf_will): 
            raise ValueError

        #descendants[i_d].mass[sf_succession[ww]] = (Ms_evolved[i_d])[ss]
        descendants[i_d].sfr[sf_succession[ww]] = (SFRs_evolved[i_d])[ss]

        is_qing = np.where(tQ[ss] < t_output[i_d]) 
        is_notqing = np.where(tQ[ss] >= t_output[i_d]) 
        descendants[i_d].sfr_class[sf_succession[is_qing]] = 'quiescent'
        descendants[i_d].sfr_class[sf_succession[is_notqing]] = 'star-forming'

        # calculate SSFR based on evloved SFR and M*
        descendants[i_d].ssfr[sf_succession[ww]] = \
                descendants[i_d].sfr[sf_succession[ww]] - descendants[i_d].mass[sf_succession[ww]]
    
        #if massevol_prop['name'] == 'integrated': 
        #    if np.min(descendant.mass[succession[sf_ancestors]] - ancestor.mass[will[sf_ancestors]]) < 0.0: 
        #        raise ValueError("Integrated mass can't reduce the mass")
    return descendants 

def Evol_IntegMstarSFR(M0, t0, tf, t_step=0.5, tQ0=None, MQ0=None, ancestor_Mq=None, **kwargs): 
    ''' Evolve stellar mass and SFR of the SF galaxies while quenching some of them
    at the same time. 

    Notes
    -----
    * SF galaxies are quenched based on the following prescription:
        The number of galaxies that start quenching between t0 and t0+tstep 
        is determined by first guessing, 

        N_quenching = N_sf(t0) * (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep
    '''
    dutycycle_prop = kwargs['dutycycle_prop']
    tau_prop = kwargs['tau_prop']
    fq_prop = kwargs['fq_prop']
    sfms_prop = kwargs['sfms_prop']
    massevol_prop = kwargs['massevol_prop']
    fudge_prop = kwargs['fudge_prop']

    Mstar = M0 
    SFR = np.repeat(-999., len(M0))
    tQ = tQ0
    Mq = MQ0
    
    t00 = t0.min()      # earliest starting time (cosmic time of ancestor snapshot) 
    t_evol = np.arange(t00, tf+t_step, t_step) 
    t_evol[-1] = tf 

    qf = Fq()
    Fq_anal = qf.model

    # Mass bins
    M_bins = np.arange(6.0, 13., 0.5) 
    M_mid = 0.5 * (M_bins[:-1] + M_bins[1:]) 

    for tt in t_evol: 
        print 't_cosmic = ', tt
        t_one_tstep = time.time()
        within = np.where(t0 <= tt)    # tsnap_genesis <= t
        sf_within = np.where((t0 <= tt) & (tQ == 999.))[0] # Not quenching SF galaxies 
        Nsf_0 = len(sf_within)          

        # M_sham of the quiescent galaxies at the closest snapshot. 
        closest_t0_snap = np.abs(t_snap-tt).argmin() - 1
        Msham_q0 = ancestor_Mq[:, closest_t0_snap]
        M_t0_sample = np.concatenate([Mstar[within], Msham_q0[np.where(Msham_q0 > 0.)]])
        
        P_sf = np.random.uniform(0., 1., Nsf_0)
        
        # Initial estimate of quenching probability given by dFq/dt 
        # P_q = (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep
        Ng0, dum = np.histogram(M_t0_sample, bins=M_bins)
        Nsf0, dum = np.histogram(Mstar[sf_within], bins=M_bins)
        P_q_arr = Ng0/Nsf0 * dFqdt(M_mid, tt + 0.5 * t_step, lit=fq_prop['name']) * t_step
        P_q_arr[np.where(Nsf0 == 0)] = 0.
        print 'M* : ', M_mid[-6:]
        print 'P_q : ', P_q_arr[-6:]
        Pq_M = interpolate.interp1d(M_mid, P_q_arr, kind='linear') 
        P_q = Pq_M(Mstar[sf_within]) 
        q_ing = np.where(P_q > P_sf)
        Nqing0 = len(q_ing[0])
        # assign them quenching times
        tQ_tmp = tQ
        tQ_tmp[sf_within[q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=Nqing0)

        kwargs_sfr = {
                't_init': t0[within], 
                't_q': tQ_tmp[within], 
                'M_q': Mq[within], 
                'dutycycle_prop': dutycycle_prop, 
                'tau_prop': tau_prop, 
                'sfms_prop': sfms_prop, 
                'indices': within
                }
        M_evol, sfr_evol, Mq_evol = M_integrate(Mstar[within], tt, tt+t_step, massevol_prop=massevol_prop, kwargs_sfr=kwargs_sfr)

        # M_sham of the quiescent galaxies at the closest snapshot. 
        closest_t_snap = np.abs(t_snap - tt - t_step).argmin() - 1 
        Msham_qf = ancestor_Mq[:, closest_t_snap]
        M_tf_sample = np.concatenate([M_evol, Msham_qf])
        
        # dPQ correction to the quenching to account for change in Ng(M*,t)
        Ngp, dum = np.histogram(M_tf_sample, bins=M_bins)
        dPq = (Ngp - Ng0)/Nsf0.astype('float') * Fq_anal(M_mid, z_of_t(tt + t_step), lit=fq_prop['name'])
        dPq[np.where(dPq < 0.)] = 0.
        dPq[np.where(Nsf0 == 0)] = 0.
        print 'dPq : ', dPq[-6:]
        dPq_M = interpolate.interp1d(M_mid, dPq, kind='linear') 
        P_q += dPq_M(Mstar[sf_within])
        
        fudge_factor = fudge_prop['slope'] * (Mstar[sf_within] - fudge_prop['fidmass']) + fudge_prop['offset']
        fudge_factor[np.where(fudge_factor < 1.)] = 1.
        P_q *= fudge_factor 
        q_ing = np.where(P_sf < P_q)
        print 'Initial guess ', Nqing0, ' after correction: ', len(q_ing[0]), ' SF galaxies out of ', Nsf_0, ' galaxies  start quenching'
        print time.time() - t_one_tstep
        print '----------------'
        # assign them quenching times then actually evolve their stellar masses  
        tQ[sf_within[q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=len(q_ing[0]))
        kwargs_sfr = {
                't_init': t0[within], 
                't_q': tQ[within], 
                'M_q': Mq[within], 
                'dutycycle_prop': dutycycle_prop, 
                'tau_prop': tau_prop, 
                'sfms_prop': sfms_prop, 
                'indices': within
                }
        M_evol, sfr_evol, Mq_evol = M_integrate(Mstar[within], tt, tt+t_step, 
                massevol_prop=massevol_prop, kwargs_sfr=kwargs_sfr)

        Mstar[within] = M_evol
        SFR[within] = sfr_evol
        Mq[within] = Mq_evol

    if SFR.min() == -999.: 
        raise ValueError

    return Mstar, SFR, tQ

def Evol_shamMstarSFR(t0, tfs, t_step=0.5, tQ0=None, MQ0=None, SFR0=None, q_Mshams=None, 
        sf_Mshams=None, quiet=True, **kwargs): 
    ''' Evolve SFR of the SF galaxies while quenching some of them at the same time. 

    Notes
    -----
    * SF galaxies are quenched based on the following prescription:
        The number of galaxies that start quenching between t0 and t0+tstep 
        is determined by first guessing, 

        N_quenching = N_sf(t0) * (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep
    '''
    dutycycle_prop = kwargs['dutycycle_prop']
    tau_prop = kwargs['tau_prop']
    fq_prop = kwargs['fq_prop']
    sfms_prop = kwargs['sfms_prop']
    massevol_prop = kwargs['massevol_prop']
    fudge_prop = kwargs['fudge_prop']

    tQ = tQ0.copy()
    Mq = MQ0.copy()
    
    t00 = t0.min()      # earliest starting time (cosmic time of ancestor snapshot) 
    t_evol = t_snap[np.abs(t_snap-np.max(tfs)).argmin()-1: np.abs(t_snap-t00).argmin()+1][::-1]

    qf = Fq()
    Fq_anal = qf.model

    # Mass bins
    M_bins = np.arange(6.0, 13., 0.5) 
    M_mid = 0.5 * (M_bins[:-1] + M_bins[1:]) 

    for i_t, tt in enumerate(t_evol[:-1]): 
        if not quiet: 
            print 't_cosmic = ', tt
            t_one_tstep = time.time()
        
        # M_sham of the quiescent galaxies at the closest snapshot. 
        closest_t0_snap = np.abs(t_snap-tt).argmin() - 1
        Msham_q0 = q_Mshams[:, closest_t0_snap].copy()
        Msham_sf0 = sf_Mshams[:, closest_t0_snap].copy()

        within = np.where((Msham_sf0 > 0.) & (t0 <= tt))    # tsnap_genesis <= t
        sf_within = np.where((Msham_sf0 > 0.) & (t0 <= tt) & (tQ == 999.))[0] # Not quenching SF galaxies 
        Nsf_0 = len(sf_within)          
        M_t0_sample = np.concatenate([Msham_sf0[within], Msham_q0[np.where(Msham_q0 > 0.)]])
        
        P_sf = np.random.uniform(0., 1., Nsf_0)
        
        # Initial estimate of quenching probability given by dFq/dt 
        # P_q = (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep
        Ng0, dum = np.histogram(M_t0_sample, bins=M_bins)
        Nsf0, dum = np.histogram(Msham_sf0[sf_within], bins=M_bins)
        P_q_arr = Ng0/Nsf0 * dFqdt(M_mid, tt + 0.5 * t_step, lit=fq_prop['name']) * t_step
        P_q_arr[np.where(Nsf0 == 0)] = 0.
        if not quiet: 
            print 'M* : ', M_mid[-6:]
            print 'P_q : ', P_q_arr[-6:]
        Pq_M = interpolate.interp1d(M_mid, P_q_arr, kind='linear') 
        P_q = Pq_M(Msham_sf0[sf_within]) 
        q_ing0 = np.where(P_q > P_sf)
        Nqing0 = len(q_ing0[0])

        # M_sham of the quiescent galaxies at the closest snapshot. 
        closest_t_snap = np.abs(t_snap - t_evol[i_t+1]).argmin() - 1 
        Msham_q1 = q_Mshams[:, closest_t_snap]
        Msham_sf1 = sf_Mshams[:, closest_t_snap]
        M_t1_sample = np.concatenate([Msham_sf1[within], Msham_q1])
        
        # dPQ correction to the quenching to account for change in Ng(M*,t)
        Ngp, dum = np.histogram(M_t1_sample, bins=M_bins)
        dPq = (Ngp - Ng0)/Nsf0.astype('float') * Fq_anal(M_mid, z_of_t(tt + t_step), lit=fq_prop['name'])
        dPq[np.where(dPq < 0.)] = 0.
        dPq[np.where(Nsf0 == 0)] = 0.
        if not quiet: 
            print 'dPq : ', dPq[-6:]
        dPq_M = interpolate.interp1d(M_mid, dPq, kind='linear') 
        P_q += dPq_M(Msham_sf0[sf_within])
        # fudge factor for P_Q
        fudge_factor = fudge_prop['slope'] * (Msham_sf0[sf_within] - fudge_prop['fidmass']) + fudge_prop['offset']
        fudge_factor[np.where(fudge_factor < 1.)] = 1.
        P_q *= fudge_factor 
        
        q_ing = np.where(P_sf < P_q)
        if not quiet: 
            print 'Initial guess ', Nqing0, ' final: ', len(q_ing[0]), ' SF gal out of ', Nsf_0, '  start quenching'
            print time.time() - t_one_tstep
            print '----------------'
        # assign them quenching times then actually evolve their stellar masses  
        tQ[sf_within[q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=len(q_ing[0]))
    
    quenching = np.where((Mq == -999.) & (tQ < 999.) & (tQ > 0))[0] 
    Nquenching = len(quenching)
    closest_tQ_index = np.abs(
            np.tile(t_snap, (Nquenching,1)) - np.tile(tQ[quenching].reshape(Nquenching,1), (1, len(t_snap)))
            ).argmin(axis=1) - 1
    Mq[quenching] = sf_Mshams[quenching, closest_tQ_index].copy()

    SFR_list = []
    for tf in tfs: 
        closest_tf_snap = np.abs(t_snap-tf).argmin() - 1
        Msham_f = sf_Mshams[:, closest_tf_snap].copy()
        if closest_tf_snap == 0: 
            before = Msham_f[:10]

        kwargs_sfr = {
                't_init': t0.copy(), 
                't_q': tQ.copy(), 
                'M_q': Mq.copy(), 
                'dutycycle_prop': dutycycle_prop, 
                'tau_prop': tau_prop, 
                'sfms_prop': sfms_prop
                }
        SFR, M_qq = logSFR_M_t(Msham_f.copy(), tf, **kwargs_sfr) 
        if closest_tf_snap == 0: 
            print 'is same? ', np.array_equal(before, Msham_f[:10])
        print tf, SFR.min(), SFR.max()

        # correct for the galaxies that were originally in the green valley. 
        # ultimately doesn't make a difference in the outcome, but to be meticulous 
        gv0 = np.where(tQ0 == 0.)   # galaxies that were initiall quenching 
        SFR[gv0] -= sfr_evol.AverageLogSFR_sfms(Mq[gv0], z_of_t(t0[gv0]), sfms_prop=sfms_prop)
        SFR[gv0] -= sfr_evol.DeltaLogSFR_quenching(tQ[gv0], t0[gv0], M_q=Mq[gv0], tau_prop=tau_prop)
        SFR[gv0] += SFR0[gv0].copy()
        SFR_list.append(SFR.copy())
        del SFR

    return SFR_list, tQ.copy()

# Deal with over quenching 
def _Overquenching(descendants, successions=None, wills=None, sf_ancestors=None): 
    for i_d, descendant in enumerate(descendants): 
        # Assign final quenched SSFR to star-forming galaxies. This is specified 
        # due to the fact that there's a lower bound on the SSFR. These values are effectively 
        # hardcoded in order to reproduce the quiescent peak of the SSFR 
        # distribution, which is more of a lower bound. 
        #q_ssfr_mean = get_q_ssfr_mean(descendant.mass[succession[sf_ancestors]])

        sf_succession = (successions[i_d])[sf_ancestors[i_d]]
        avg_q_ssfr = sfr_evol.AverageLogSSFR_q_peak(descendant.mass[sf_succession])
        sigma_q_ssfr = sfr_evol.ScatterLogSSFR_q_peak(descendant.mass[sf_succession])

        min_q_ssfr = sigma_q_ssfr * np.random.randn(len(sf_succession)) + avg_q_ssfr 
        descendants[i_d].min_ssfr[sf_succession] = min_q_ssfr 

        # Deal with over quenched galaxies 
        overquenched = np.where(
                descendants[i_d].min_ssfr[sf_succession] > descendants[i_d].ssfr[sf_succession]
                )
        if len(overquenched[0]) > 0: 
            descendants[i_d].ssfr[sf_succession[overquenched]] = \
                    descendant.min_ssfr[sf_succession[overquenched]]
            descendants[i_d].sfr[sf_succession[overquenched]] = \
                    descendant.ssfr[sf_succession[overquenched]] \
                    + descendant.mass[sf_succession[overquenched]]
    return descendants 


# SF evolution (Quenching probability based) 
def _StarformingEvol_Pq(ancestor, descendant, succession=None, will=None, 
        sfr_prop=None, evol_prop=None, lineage=None): 
    ''' Evolve the stellar mass and star formation rates simultaneously. 
    The stellar mass evolution is dictated by the SFR.
    '''
    sfms_prop = sfr_prop['sfms']
    pq_prop = evol_prop['pq'] 
    tau_prop = evol_prop['tau']
    sfrevol_prop = evol_prop['sfr']
    dutycycle_prop = sfrevol_prop['dutycycle']
    massevol_prop = evol_prop['mass']

    # quenching probability 
    pq_time = time.time() 
    P_q = _getPq(ancestor, descendant, 
            succession=succession, will=will, 
            pq_prop=pq_prop, sfr_prop=sfr_prop)
    print 'Quenching Probabilities takes : ', time.time() - pq_time
    
    # quenching time/redshift
    tq_time = time.time()
    descendant, t_q, z_q = _tQuench(ancestor, descendant, 
            P_q, succession=succession, will=will)
    print 'tQuench takes: ', time.time() - tq_time     
    
    # star forming ancestors
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]  # SF ancestors

    # star-formation duty cycle parameters
    dutycycle_time = time.time()
    dutycycle_prop = sfr_evol.dutycycle_param(len(sf_ancestors), dutycycle_prop=dutycycle_prop)
    dutycycle_prop['delta_sfr'] = ancestor.delta_sfr[will[sf_ancestors]]
    print 'SFR dutycycle properties take ', time.time() - dutycycle_time, ' to generate'
    
    M_q = np.repeat(-999., len(sf_ancestors))
    # keywords for logSFR(M*,t)
    kwargs_sfr = {
            't_init': descendant.tsnap_genesis[succession[sf_ancestors]], 
            't_q': t_q, 
            'z_q': z_q, 
            'M_q': M_q, 
            'dutycycle_prop': dutycycle_prop, 
            'tau_prop': tau_prop, 
            'sfms_prop': sfms_prop
            }

    sham_mass =  descendant.mass[succession[sf_ancestors]].copy()   # descendant SHAM masses

    if massevol_prop['name'] == 'integrated':   # integrated stellar mass 
        mass_time = time.time()
        descendant.mass[succession[sf_ancestors]], descendant.sfr[succession[sf_ancestors]], blah = \
                M_integrate(
                        ancestor.mass_genesis[will[sf_ancestors]], 
                        ancestor.tsnap_genesis[will[sf_ancestors]], 
                        descendant.t_cosmic,
                        massevol_prop=massevol_prop,
                        kwargs_sfr=kwargs_sfr               # kwargs for logSFR(M*,t) function
                        )
        if np.min(descendant.mass[succession[sf_ancestors]] - ancestor.mass[will[sf_ancestors]]) < 0.0: 
            raise ValueError("Integrated mass can't reduce the mass")
        print 'Integrated Masses takes ', time.time() - mass_time
        
    elif massevol_prop['name'] == 'sham':       # SHAM masses 
        # determine the quenching stellar mass (the SHAM stellar mass at t_Q)
        M_q = _SHAM_Mquenching(kwargs_sfr['t_q'], lineage=lineage, 
                descendant=descendant, 
                descendant_index=succession[sf_ancestors])
        kwargs_sfr['M_q'] = M_q

        descendant.sfr[succession[sf_ancestors]] = \
                logSFR_M_t(sham_mass, descendant.t_cosmic, **kwargs_sfr)
    
    # calculate SSFR based on evloved SFR and M*
    descendant.ssfr[succession[sf_ancestors]] = \
            descendant.sfr[succession[sf_ancestors]] - descendant.mass[succession[sf_ancestors]]
    
    return 

def _getPq(ancestor, descendant, succession=None, will=None, pq_prop=None, sfr_prop=None):
    '''' Calculate the quenching probability for Star-Forming ancestor galaxies
    in the ancestor CGPop object. The quenching probability is calculated using 

    P_Q = ( f_Q(Mf,zf) - f_Q(M0,z0) ) / (1 - f_Q(M0,z0)) + Offset_P_Q
    
    where Offset_P_Q is a fudge factor dictated by nuisance parameters described 
    by pq_prop. 
    
    Parameters
    ----------
    pq_prop : dict
        Quenching properties dictionary. Keywords are 'slope' and 'yint' which
        describes the offset to the quenchign probability. The quenchign probability 
        is fudged around to match the observed quenching fraction at z_final: 
        fQ(M*, z_final). 

    Returns
    -------
    P_q : array
        Array that specifies the quenching probabilities for the SF ancestor 
        galaxies. 
    '''
    z_final = descendant.zsnap  # z_final of evolution
    # SF ancestors
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]      

    P_q_offset = pq_prop['slope'] * (descendant.mass[succession[sf_ancestors]] - 9.5) + pq_prop['yint']
    
    # fQ(M*, z_final)
    fq_obj = Fq()
    fqf = fq_obj.model(descendant.mass[succession[sf_ancestors]], 
            z_final, 
            lit = sfr_prop['fq']['name'])
    # fQ(M*, z_initial), which is calculated using mass and redshifts when 
    # the host subhalo passes the M* threshold. 
    fq0 = fq_obj.model(descendant.mass_genesis[succession[sf_ancestors]], 
            descendant.zsnap_genesis[succession[sf_ancestors]], 
            lit=sfr_prop['fq']['name'])

    notallquench = np.where(fq0 < 1.0)
    P_q = np.repeat(1.0, len(sf_ancestors))
    P_q[notallquench] = (fqf[notallquench] - fq0[notallquench] ) / (1.0 - fq0[notallquench]) 
    P_q += P_q_offset

    return P_q

def _tQuench(ancestor, descendant, P_q, succession=None, will=None):
    ''' Given ancestor and descendent galaxy populations AND the quenching
    probabilities, calculate the cosmic time when the star-forming 
    ancestor galaxies, begins to quench. The the quenching time is 
    sampled uniformly between t_final and t_snap_gensis.
    '''
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]      # SF ancestors

    np.random.seed()
    randoms = np.random.uniform(0., 1., len(sf_ancestors)) 
    is_qing = np.where(P_q > randoms)           # quenching
    is_notqing = np.where(P_q <= randoms)       # not quenching 

    descendant.sfr_class[succession[sf_ancestors[is_qing]]] = 'quiescent'
    descendant.sfr_class[succession[sf_ancestors[is_notqing]]] = 'star-forming'
    #print len(is_qing[0]), ' is quenching'
    #print len(is_notqing[0]), ' is not quenching'

    np.random.seed()
    t_final = descendant.t_cosmic       # final t_cosmic 
    t_q = (t_final - ancestor.tsnap_genesis[will[sf_ancestors]]) * np.random.uniform(0., 1., len(sf_ancestors))
    t_q += ancestor.tsnap_genesis[will[sf_ancestors]]
    z_q = z_of_t(t_q)
    t_q[is_notqing] = 999.0     # galaxies that will never quench
    z_q[is_notqing] = -999.0

    return descendant, t_q, z_q



def MstarSFR_simul_evol_OBSOLETE(M0, t0, tf, t_step=0.2, ancestor_Mq=None, quiet=True, **kwargs): 
    ''' Evolve stellar mass, SFR, and quench galaxies simultaneously. 

    Notes
    -----
    * SF galaxies are quenched based on the following prescription:
        The number of galaxies that start quenching between t0 and t0+tstep 
        is determined by, 

        N_quenching = N_sf(t0) * (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep

    '''
    dutycycle_prop = kwargs['dutycycle_prop']
    tau_prop = kwargs['tau_prop']
    fq_prop = kwargs['fq_prop']
    sfms_prop = kwargs['sfms_prop']
    massevol_prop = kwargs['massevol_prop']

    Mstar = M0 
    SFR = np.repeat(-999., len(M0))
    tQ = np.repeat(999., len(M0))
    Mq = np.repeat(-999., len(M0))
    
    t00 = t0.min()      # earliest starting time (cosmic time of ancestor snapshot) 
    t_evol = np.arange(t00, tf+t_step, t_step) 
    t_evol[-1] = tf 

    qf = Fq()
    Fq_anal = qf.model

    for tt in t_evol: 
        if not quiet: 
            print 't_cosmic = ', tt
        within = np.where(t0 <= tt)    # tsnap_genesis <= t
        sf_within = np.where((t0 <= tt) & (tQ == 999.))[0] # Not quenching SF galaxies 
        Nsf_0 = len(sf_within)          
        
        P_sf = np.random.uniform(0., 1., Nsf_0)
    
        # Initial estimate of P_q = (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep
        t_offset = np.zeros(len(sf_within))
        #t_offset = -2.7 * (Mstar[sf_within] - 11.5)
        #t_offset[np.where(t_offset < 0.)] = 0.
        fudge = 1.  # fudge factor
        f_Ng = (1./(1.- Fq_anal(Mstar[sf_within], z_of_t(tt), lit=fq_prop['name'])))
        whereinf = np.where(f_Ng == np.inf) 
        P_q = fudge * f_Ng * dFqdt(Mstar[sf_within], tt + t_offset + 0.5 * t_step, lit=fq_prop['name']) * t_step    
        P_q[whereinf] = 1.0
        q_ing = np.where(P_q > P_sf)
        if not quiet: 
            print 'P_q', P_q.min(), P_q.max(), P_q.mean()
            print 'Initial guess ', len(q_ing[0]), ' SF galaxies out of ', Nsf_0, ' galaxies  start quenching'

        # assign them quenching times
        tQ[sf_within[q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=len(q_ing[0]))

        kwargs_sfr = {
                't_init': t0[within], 
                't_q': tQ[within], 
                'M_q': Mq[within], 
                'dutycycle_prop': dutycycle_prop, 
                'tau_prop': tau_prop, 
                'sfms_prop': sfms_prop, 
                'indices': within
                }
        M_evol, sfr_evol, Mq_evol = M_integrate(Mstar[within], tt, tt+t_step, massevol_prop=massevol_prop, kwargs_sfr=kwargs_sfr)

        Mstar[within] = M_evol
        SFR[within] = sfr_evol
        Mq[within] = Mq_evol

    if SFR.min() == -999.: 
        raise ValueError

    return Mstar, SFR, tQ


# Universal functions 
def logSFR_M_t(logmass, t_input, t_init=None, t_q=None, M_q=None, 
        dutycycle_prop=None, tau_prop=None, sfms_prop=None, indices=None): 
    ''' log(SFR) as a function of M* and t_cosmic.

    Notes
    -----
    * kwarg indices is for dutycycle_prop 
    '''
    # SFR evolution based on solving an ODE of SFR
    logsfr_time = time.time()
    z_init = z_of_t(t_init) # initial redshift
    
    # update quenched M* ( this is the stellar mass roughly at the time of quenching)
    just_quenched = np.where(
            (t_input > t_q) & 
            (M_q == -999.))
    if len(just_quenched[0]) > 0: 
        M_q[just_quenched] = logmass[just_quenched]

    # average SFR of SFMS at M* and z_init
    quenched = np.where(t_input > t_q)
    tmp_M = logmass.copy() 
    tmp_M[quenched] = M_q[quenched]
    avglogsfr = sfr_evol.AverageLogSFR_sfms(tmp_M, z_init, sfms_prop=sfms_prop)

    # log(SFR)_SFMS evolutionfrom t0
    logsfr_sfms = sfr_evol.DeltaLogSFR_sfms(
            z_init, 
            z_of_t(t_input),
            sfms_prop=sfms_prop)

    # log(SFR)_duty cycle evolution from t0 to tQ
    logsfr_sfduty = sfr_evol.DeltaLogSFR_dutycycle(
            t_init, 
            t_input, 
            t_q=t_q, 
            dutycycle_prop=dutycycle_prop, 
            indices=indices)

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
        #print niter, ' ', type, ' iterations'
        #print 'f_reatin = ', f_retain, 'delta_t = ', delt
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

    return logM_n_1, logSFR_n_1, M_q_new

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
