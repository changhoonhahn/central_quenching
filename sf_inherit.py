""" 

Evolve star forming properties of ancestor CenQue objects within the 
Lineage object ancestor for the descendant CenQue objects

"""
import time
import numpy as np
from scipy import interpolate

from lineage import Lineage

import sfr_evol

from gal_prop import Fq
from gal_prop import dFqdt
from util.util import get_nsnap_t
from util.util import intersection_index

import matplotlib.pyplot as plt
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors


# spline between z and t_cosmic
z_snap, t_snap = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
z_of_t = interpolate.interp1d(list(reversed(t_snap)), list(reversed(z_snap)), kind='cubic') 

def InheritSF(nsnap_descendant, nsnap_ancestor=20, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = {
            'fq': {'name': 'wetzelsmooth'}, 
            'sfms': {'name': 'lineage', 'mslope': 0.55, 'zslope': 1.1}},
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
    
    # read in the lineage (< 0.05 seconds for one snapshot)
    read_time = time.time()
    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read(range(nsnap_descendant, nsnap_ancestor))
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
    descendant.sfr_prop = ancestor.sfr_prop 
    
    succession, will = intersection_index(
            getattr(descendant, 'ancestor'+str(nsnap_ancestor)), 
            ancestor.snap_index)
    if len(succession) != len(descendant.mass): 
        raise ValueError('Something wrong with the lineage')
    q_ancestors = np.where(ancestor.sfr_class[will] == 'quiescent')[0]      # Q ancestors
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]  # SF ancestors
    print "nUmber of ancestors"
    print len(q_ancestors), len(sf_ancestors)
    #ancestor.plotFq(model=True, savefig='figure/test/test_fq_anc.png')      # testing Fq
    #ancestor.plotSFMS(sfqcut=True, allgal=True, bovyplot=False, sigSFR=False, model=False, savefig='figure/test/test_sfms_anc.png')      # testing SFMS 
    
    # Evolve queiscent ancestor galaxies 
    q_time = time.time()
    descendant = _QuiescentEvol(ancestor, descendant, succession=succession, will=will)
    print 'Quiescent evolution takes ', time.time()-q_time
    print np.sum(descendant.sfr_class == 'quiescent') 
    #descendant.plotFq(model=True, savefig='figure/test/test_fq_desc_afterQevol.png')      # testing Fq
    #descendant.plotSFMS(sfqcut=True, allgal=True, bovyplot=False, sigSFR=False, model=False, savefig='figure/test/test_sfms_desc_afterQevol.png')      # testing Fq

    # Evolve Star Forming Galaxies 
    sf_time = time.time()
    if evol_prop['type'] == 'pq_based':
        _StarformingEvol_Pq(ancestor, descendant, succession=succession, will=will, 
                evol_prop=evol_prop, sfr_prop=sfr_prop, 
                lineage=bloodline)
    elif evol_prop['type'] == 'simult': 
        _StarformingEvol_SimulEvo(ancestor, descendant, succession=succession, will=will, 
                evol_prop=evol_prop, sfr_prop=sfr_prop, 
                lineage=bloodline)
    elif evol_prop['type'] == 'predict': 
        _StarformingEvol_Predict(ancestor, descendant, succession=succession, will=will, 
                evol_prop=evol_prop, sfr_prop=sfr_prop, 
                lineage=bloodline)
    print 'Star Forming Evolution takes ', time.time() - sf_time
    #descendant.plotFq(model=True, savefig='figure/test/test_fq_desc_afterSFevol.png')      # testing Fq
    #descendant.plotSFMS(sfqcut=True, allgal=True, bovyplot=False, sigSFR=False, model=False, savefig='figure/test/test_sfms_desc_afterSFevol.png')      # testing Fq

    # Assign final quenched SSFR to star-forming galaxies. This is specified 
    # due to the fact that there's a lower bound on the SSFR. These values are effectively 
    # hardcoded in order to reproduce the quiescent peak of the SSFR 
    # distribution, which is more of a lower bound. 
    #q_ssfr_mean = get_q_ssfr_mean(descendant.mass[succession[sf_ancestors]])
    avg_q_ssfr = sfr_evol.AverageLogSSFR_q_peak(descendant.mass[succession[sf_ancestors]])
    sigma_q_ssfr = sfr_evol.ScatterLogSSFR_q_peak(descendant.mass[succession[sf_ancestors]])

    min_q_ssfr = sigma_q_ssfr * np.random.randn(len(sf_ancestors)) + avg_q_ssfr 
    descendant.min_ssfr[succession[sf_ancestors]] = min_q_ssfr 
    #print descendant.min_ssfr[succession[sf_ancestors]]

    # Deal with over quenched galaxies 
    overquenched = np.where(
            descendant.min_ssfr[succession[sf_ancestors]] > descendant.ssfr[succession[sf_ancestors]]
            )
    if len(overquenched[0]) > 0: 
        descendant.ssfr[succession[sf_ancestors[overquenched]]] = descendant.min_ssfr[succession[sf_ancestors[overquenched]]]
        descendant.sfr[succession[sf_ancestors[overquenched]]] = descendant.ssfr[succession[sf_ancestors[overquenched]]] \
                + descendant.mass[succession[sf_ancestors[overquenched]]]

    descendant.data_columns = list(descendant.data_columns) + ['ssfr', 'sfr', 'min_ssfr', 'sfr_class']
    setattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant), descendant)

    #descendant.plotFq(model=True, savefig='figure/test/test_fq_desc_afterOverQuenching.png')      # testing Fq
    #descendant.plotSFMS(sfqcut=True, allgal=True, bovyplot=False, model=False, sigSFR=False, savefig='figure/test/test_sfms_desc_afterOverQuenching.png')      # testing Fq

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

# SF Evolution (Simultaneous evolution) 
def _StarformingEvol_SimulEvo(ancestor, descendant, succession=None, will=None, 
        sfr_prop=None, evol_prop=None, lineage=None): 
    ''' Evolve the stellar mass, star formation rates, and quench galaxies simultaneously. 
    '''
    sfms_prop = sfr_prop['sfms']        # SFMS properties
    fq_prop = sfr_prop['fq']            # fq_prop 
    tau_prop = evol_prop['tau']         # Quenching timescale properties
    dutycycle_prop = evol_prop['sfr']['dutycycle']   # SF dutycycle properties 
    massevol_prop = evol_prop['mass']   # mass evolution properties
    if massevol_prop['name'] != 'integrated': 
        raise ValueError('Only supported for integrated') 
    
    # star forming ancestors
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]  # SF ancestors
    q_ancestors = np.where(ancestor.sfr_class[will] == 'quiescent')[0]  # SF ancestors

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
    MshamEvol = ancestor.Msham_evol[0]

    mass_time = time.time()
    M_evolved, SFR_evolved, tQ = \
            MstarSFR_simul_evol(
                    ancestor.mass_genesis[will[sf_ancestors]], 
                    ancestor.tsnap_genesis[will[sf_ancestors]], 
                    descendant.t_cosmic,
                    ancestor_Mq=MshamEvol[will[q_ancestors],:], 
                    **kwargs_sfr               
                    )
    descendant.mass[succession[sf_ancestors]] = M_evolved
    descendant.sfr[succession[sf_ancestors]] = SFR_evolved
    is_notqing = np.where(tQ == 999.)
    is_qing = np.where(tQ < 999.) 
    descendant.sfr_class[succession[sf_ancestors[is_qing]]] = 'quiescent'
    descendant.sfr_class[succession[sf_ancestors[is_notqing]]] = 'star-forming'

    #print descendant.mass[succession[sf_ancestors]][:100]
    #print descendant.sfr[succession[sf_ancestors]][:100]
    if np.min(descendant.mass[succession[sf_ancestors]] - ancestor.mass[will[sf_ancestors]]) < 0.0: 
        raise ValueError("Integrated mass can't reduce the mass")
    # calculate SSFR based on evloved SFR and M*
    descendant.ssfr[succession[sf_ancestors]] = \
            descendant.sfr[succession[sf_ancestors]] - descendant.mass[succession[sf_ancestors]]
    return     

def MstarSFR_simul_evol(M0, t0, tf, t_step=0.2, ancestor_Mq=None, **kwargs): 
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
        print 'P_q', P_q.min(), P_q.max(), P_q.mean()
        q_ing = np.where(P_q > P_sf)
        print 'Initial guess ', len(q_ing[0]), ' SF galaxies out of ', Nsf_0, ' galaxies  start quenching'

        # assign them quenching times
        tQ_tmp = tQ
        tQ_tmp[sf_within[q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=len(q_ing[0]))

        # M_sham of the quiescent galaxies at the closest snapshot. 
        closest_t0_snap = list(np.abs(t_snap-tt)).index(np.abs(t_snap - tt).min()) - 1
        Msham_q0 = ancestor_Mq[:, closest_t0_snap]
        M_t0_sample = np.concatenate([Mstar[within], Msham_q0[np.where(Msham_q0 > 0.)]])

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
        closest_t_snap = list(np.abs(t_snap-tt-t_step)).index(np.abs(t_snap - tt-t_step).min()) - 1
        Msham_qf = ancestor_Mq[:, closest_t_snap]
        M_tf_sample = np.concatenate([M_evol, Msham_qf[np.where(Msham_qf > 0.)]])
        
        # dPQ correction to the quenching to account for change in Ng(M*,t)
        M_bins = np.arange(
                np.min([M_t0_sample.min(), M_tf_sample.min()])-0.5, 
                np.max([M_t0_sample.max(), M_tf_sample.max()])+1.0, 0.5) 
        M_low = M_bins[:-1]
        M_high = M_bins[1:]
        dPq = np.zeros(len(M_low))
        for i_m in range(len(M_low)): 
            gal_0 = np.where((M_t0_sample >= M_low[i_m]) & (M_t0_sample < M_high[i_m]))
            gal_p = np.where((M_tf_sample >= M_low[i_m]) & (M_tf_sample < M_high[i_m]))
            gal_sf_0 = np.where((t0 <= tt) & (tQ == 999.) & (Mstar >= M_low[i_m]) & (Mstar < M_high[i_m])) 

            Ng0 = np.float(len(gal_0[0]))
            Ngp = np.float(len(gal_p[0]))
            Nsf0 = np.float(len(gal_sf_0[0]))
            if Nsf0 == 0: 
                continue

            #t_offset = -2.7 * (0.5*(M_high[i_m]+M_low[i_m]) - 11.5)
            #if t_offset < 0.: 
            #    t_offset = 0.
            #print t_offset, z_of_t(tt+t_step+t_offset)
            fq_tf = Fq_anal(0.5*(M_high[i_m]+M_low[i_m]), z_of_t(tt+t_step), lit=fq_prop['name'])
            dPq[i_m] = (Ngp - Ng0) * (1. - fq_tf)/Nsf0

        dPq[np.where(dPq < 0.)] = 0.
        dPq_M = interpolate.interp1d(0.5*(M_low + M_high), dPq, kind='linear') 
        print 'dPq', dPq_M(Mstar[sf_within]).min(), dPq_M(Mstar[sf_within]).max(), np.mean(dPq_M(Mstar[sf_within]))
        P_q += dPq_M(Mstar[sf_within])
        
        q_ing = np.where(P_q > P_sf)
        print 'After correction: ', len(q_ing[0]), ' SF galaxies out of ', Nsf_0, ' galaxies  start quenching'
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

# SF Evolution (Predictive) 
def _StarformingEvol_Predict(ancestor, descendant, succession=None, will=None, 
        sfr_prop=None, evol_prop=None, lineage=None): 
    ''' Evolve the stellar mass, star formation rates, and quench galaxies simultaneously. 
    '''
    sfms_prop = sfr_prop['sfms']        # SFMS properties
    fq_prop = sfr_prop['fq']            # fq_prop 
    tau_prop = evol_prop['tau']         # Quenching timescale properties
    dutycycle_prop = evol_prop['sfr']['dutycycle']   # SF dutycycle properties 
    massevol_prop = evol_prop['mass']   # mass evolution properties
    
    # star forming ancestors
    sf_ancestors = np.where(ancestor.sfr_class[will] == 'star-forming')[0]  # SF ancestors
    q_ancestors = np.where(ancestor.sfr_class[will] == 'quiescent')[0]  # Q ancestor 

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
    MshamEvol = ancestor.Msham_evol[0]

    mass_time = time.time()
    M_evolved, SFR_evolved, tQ = \
            MstarSFR_predict_evol(
                    ancestor.mass_genesis[will[sf_ancestors]], 
                    ancestor.tsnap_genesis[will[sf_ancestors]], 
                    descendant.t_cosmic,
                    SFR0=ancestor.sfr[will[sf_ancestors]],
                    ancestor_Mq=MshamEvol[will[q_ancestors],:], 
                    **kwargs_sfr               
                    )
    descendant.mass[succession[sf_ancestors]] = M_evolved
    descendant.sfr[succession[sf_ancestors]] = SFR_evolved
    is_notqing = np.where(tQ == 999.)
    is_qing = np.where(tQ < 999.) 
    descendant.sfr_class[succession[sf_ancestors[is_qing]]] = 'quiescent'
    descendant.sfr_class[succession[sf_ancestors[is_notqing]]] = 'star-forming'

    #print descendant.mass[succession[sf_ancestors]][:100]
    #print descendant.sfr[succession[sf_ancestors]][:100]
    if np.min(descendant.mass[succession[sf_ancestors]] - ancestor.mass[will[sf_ancestors]]) < 0.0: 
        raise ValueError("Integrated mass can't reduce the mass")
    # calculate SSFR based on evloved SFR and M*
    descendant.ssfr[succession[sf_ancestors]] = \
            descendant.sfr[succession[sf_ancestors]] - descendant.mass[succession[sf_ancestors]]
    return     

def MstarSFR_predict_evol(M0, t0, tf, t_step=0.2, SFR0=None, ancestor_Mq=None, **kwargs): 
    ''' Evolve stellar mass, SFR, and quench galaxies predictively. 
    In other words, determine how many galaxies need to be quenched 
    in order to match the quiescent fraction at t_step timesteps.  

    Parameters
    ----------
    M0 : (np.ndarray)
        Numpy nd array that specifies the M* at "genesis"
    t0 : (np.ndarray)
        Numpy nd array that specifies the cosmic time at "genesis"
    tf : (float)
        Float that specifies the cosmic time of the descendant. 
        Final t_cosmic. 
    t_step : (float) 
        time step of the evolution  
    ancestor_Mq : (np.ndarray)
        Numpy nd array with dimensions q_ancestors x ancestor_snapshot-1 that 
        specifies the M_sham for all snapshots (used to measure fQ_obvs). 

    Notes
    -----
    * SF galaxies are quenched based on the following prescription:
        The number of galaxies that start quenching between t0 and t0+tstep 
        is determined by, 

        N_quenching = N_sf(t0) * ( Fq(t0 + 0.5 tstep) - Fq_obv(t0) )/ 0.5 tstep * tstep 
    
    The main problem with this scheme is that if the quenching timescale 
    is significantly longer than the timesteps, the change in Fq_obv
    will be very small, causing a lot of galaxies to be quenched. 

    This may be solved by considering both Fq_obv from SFR_cut and 
    Fq_obv based on sfr_class, that way galaxies that are in the process
    of quenching is reflected.
    '''
    dutycycle_prop = kwargs['dutycycle_prop']
    tau_prop = kwargs['tau_prop']
    fq_prop = kwargs['fq_prop']
    sfms_prop = kwargs['sfms_prop']
    massevol_prop = kwargs['massevol_prop']

    Mstar = M0 
    SFR = SFR0
    tQ = np.repeat(999., len(M0))
    Mq = np.repeat(-999., len(M0))
    
    t00 = t0.min()      # earliest starting time (cosmic time of ancestor snapshot) 
    t_evol = np.arange(t00, tf+t_step, t_step) 
    t_evol[-1] = tf 

    qf = Fq()
    Fq_anal = qf.model
    Fq_cut = qf.SFRcut

    for tt in t_evol: 
        print 't_cosmic = ', tt
        within = np.where(t0 <= tt)[0]    # tsnap_genesis <= t
        sf_within = np.where((t0 <= tt) & (tQ == 999.))[0] # Not quenching SF galaxies 
        Nsf_0 = len(sf_within)          
        
        P_sf = np.random.uniform(0., 1., Nsf_0)

        # M_sham of the quiescent galaxies at the closest snapshot. 
        closest_t0_snap = list(np.abs(t_snap-tt)).index(np.abs(t_snap - tt).min()) - 1
        Msham_q0 = ancestor_Mq[:, closest_t0_snap]
        M_t0_sample = np.concatenate([Mstar[within], Msham_q0[np.where(Msham_q0 > 0.)]])

        # calculate the observed quiescent fraction fQ_obvs
        M_bins = np.arange(M_t0_sample.min()-0.5, M_t0_sample.max()+1.0, 0.5) 
        M_low = M_bins[:-1]
        M_high = M_bins[1:]

        fq_obvs = np.zeros(len(M_low)) 
        for i_m in range(len(M_low)):  
            mbin = np.where((Mstar[within] >= M_low[i_m]) & (Mstar[within] < M_high[i_m]))
            q_mbin = np.where((Msham_q0 >= M_low[i_m]) & (Msham_q0 < M_high[i_m]))
            sf_mbin = np.where((Mstar[sf_within] >= M_low[i_m]) & (Mstar[sf_within] < M_high[i_m]))

            gal_sf = np.where(SFR[within[mbin]] > Fq_cut(Mstar[within[mbin]], z_of_t(tt), sfms_prop=sfms_prop))
            Nsf = np.float(len(gal_sf[0])) 

            if np.float(len(mbin[0])+len(q_mbin[0])) == 0: 
                continue 
            #fq_tmp = 1 - np.float(len(sf_mbin[0]))/np.float(len(mbin[0])+len(q_mbin[0]))
            fq_tmp = 1 - (0.5 * (Nsf+np.float(len(sf_mbin[0]))))/np.float(len(mbin[0])+len(q_mbin[0]))
            fq_obvs[i_m] = fq_tmp #np.float(len(q_mbin[0]))/np.float(len(sf_mbin[0])+len(q_mbin[0])))
        print 0.5 * (M_low + M_high) 
        print fq_obvs

        # Fq_obvs function  
        fq_obvs_M = interpolate.interp1d(0.5 * (M_low + M_high), fq_obvs, kind='linear') 

        f_Ng = (1./(1.- Fq_anal(Mstar[sf_within], z_of_t(tt), lit=fq_prop['name'])))
        whereinf = np.where(f_Ng == np.inf) 

        dfqdt_obvs = (Fq_anal(Mstar[sf_within], z_of_t(tt+t_step), lit=fq_prop['name']) - fq_obvs_M(Mstar[sf_within]))/t_step
        dfqdt_theo = dFqdt(Mstar[sf_within], tt + t_step, lit=fq_prop['name'])
        dfqdt_obvs[np.where(dfqdt_obvs < dfqdt_theo)] = dfqdt_theo[np.where(dfqdt_obvs < dfqdt_theo)]
        #print 'M*', Mstar[sf_within][:5]
        #print 'dFq/dt observed: ', dfqdt_obvs[:5]
        #print 'dFq/dt theoreti: ', dFqdt(Mstar[sf_within][:5], tt + t_step, lit=fq_prop['name'])
        P_q = f_Ng * dfqdt_obvs * t_step    

        P_q[whereinf] = 1.0
        print 'P_q', P_q.min(), P_q.max(), P_q.mean()
        q_ing = np.where(P_q > P_sf)
        print 'Initial guess ', len(q_ing[0]), ' SF galaxies out of ', Nsf_0, ' galaxies  start quenching'

        # assign them quenching times
        tQ_tmp = tQ
        tQ_tmp[sf_within[q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=len(q_ing[0]))

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
        closest_t_snap = list(np.abs(t_snap-tt-t_step)).index(np.abs(t_snap - tt-t_step).min()) - 1
        Msham_qf = ancestor_Mq[:, closest_t_snap]
        M_tf_sample = np.concatenate([M_evol, Msham_qf[np.where(Msham_qf > 0.)]])

        # dPQ correction to the quenching to account for change in Ng(M*,t)
        M_bins = np.arange(
                np.min([M_t0_sample.min(), M_tf_sample.min()])-0.5, 
                np.max([M_t0_sample.max(), M_tf_sample.max()])+1.0, 0.5) 
        M_low = M_bins[:-1]
        M_high = M_bins[1:]
        dPq = np.zeros(len(M_low))
        for i_m in range(len(M_low)): 
            gal_0 = np.where((M_t0_sample >= M_low[i_m]) & (M_t0_sample < M_high[i_m]))
            gal_p = np.where((M_tf_sample >= M_low[i_m]) & (M_tf_sample < M_high[i_m]))
            gal_sf_0 = np.where((t0 <= tt) & (tQ == 999.) & (Mstar >= M_low[i_m]) & (Mstar < M_high[i_m])) 

            Ng0 = np.float(len(gal_0[0]))
            Ngp = np.float(len(gal_p[0]))
            Nsf0 = np.float(len(gal_sf_0[0]))
            if Nsf0 == 0: 
                continue

            t_offset = -2.7 * (0.5*(M_high[i_m]+M_low[i_m]) - 11.5)
            if t_offset < 0.: 
                t_offset = 0.

            fq_tf = Fq_anal(0.5*(M_high[i_m]+M_low[i_m]), z_of_t(tt+t_step), lit=fq_prop['name'])
            dPq[i_m] = (Ngp - Ng0) * (1. - fq_tf)/Nsf0

        dPq[np.where(dPq < 0.)] = 0.
        dPq_M = interpolate.interp1d(0.5*(M_low + M_high), dPq, kind='linear') 
        print 'dPq', dPq_M(Mstar[sf_within]).min(), dPq_M(Mstar[sf_within]).max(), np.mean(dPq_M(Mstar[sf_within]))
        P_q += dPq_M(Mstar[sf_within])
        
        q_ing = np.where(P_q > P_sf)
        print 'After correction: ', len(q_ing[0]), ' SF galaxies out of ', Nsf_0, ' galaxies  start quenching'
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

def logSFR_M_t(logmass, t_input, t_init=None, t_q=None, z_q=None, M_q=None, 
        dutycycle_prop=None, tau_prop=None, sfms_prop=None, indices=None): 
    ''' log(SFR) as a function of M* and t_cosmic.

    Notes
    -----
    * kwarg indices is for dutycycle_prop 
    '''
    # SFR evolution based on solving an ODE of SFR
    logsfr_time = time.time()
    z_init = z_of_t(t_init) # initial redshift
    if z_q is None: 
        z_q = np.repeat(-999., len(t_q))
        qing = np.where(t_q != 999.)
        z_q[qing] = z_of_t(t_q[qing]) 
    
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
