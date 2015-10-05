"""

Evolve CenQue object over cosmic time 

Author(s): ChangHoon Hahn

"""
import time
import random 
import numpy as np
import warnings

from cenque import CenQue
from util import cenque_utility as util
from assign_sfr import assign_sfr
from quiescent_fraction import get_fq
from quiescent_fraction import get_fq_nsnap
from util.gal_classify import sfr_cut 
from util.gal_classify import sfq_classify
from sfms.fitting import get_param_sfr_mstar_z
from util.tau_quenching import get_quenching_efold

def evolve_cq(
        cenque, 
        final_nsnap = 1, 
        sf_prop = {'name': 'average'}, 
        fq_prop = {'name': 'wetzelsmooth'}, 
        tau_prop = {'name': 'instant'}, 
        mass_evol = 'sham', 
        quiet = False, 
        **kwargs
        ): 
    """ Evolve SF properties and stellar mass properties of galaxies 
    in the CenQue class object over cosmic time until final_nsnap.
    
    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    cenque : CenQue class object
    final_nsnap : final snapshot number 

    ----------------------------------------------------------------
    Notes
    ----------------------------------------------------------------
    * Quenching of SF galaxies now involve quenching-efold *and* 
    overall SF MS SFR decrease

    """
    evo_start_time = time.time()
    
    if cenque.cenque_type != 'sf_assigned': 
        raise ValueError()

    if not quiet: 
        print 'Quiescent Fraction', 
        print np.float(len(cenque.gal_type[cenque.gal_type == 'quiescent']))/np.float(len(cenque.gal_type)) 
    
    start_nsnap = cenque.nsnap 
    parent_cq = cenque

    for i_snap in range(final_nsnap, start_nsnap)[::-1]:    

        print ''
        print '---------------------------------------'
        print 'Evolving to ', str(i_snap) 

    
        # Import halo and SHAM properties from TreePM 
        # catalogs and run time evolution on the SF
        # properties of the parent CenQue object.
        child_cq = CenQue() 
        child_cq.import_treepm(i_snap) 
        child_cq.cenque_type = ''.join(['evol_from', str(start_nsnap)])
        child_cq.sf_prop = sf_prop
        child_cq.fq_prop = fq_prop
        child_cq.tau_prop = tau_prop
        child_cq.mass_evol = mass_evol

        child_cq = evolve_onestep(parent_cq, child_cq, quiet=quiet)

        child_cq.writeout()        

        parent_cq = child_cq
    
    print 'Total Evolution takes ', time.time() - evo_start_time

    return child_cq 

def evolve_onestep(parent_cq, child_cq, quiet=False):
    """ Evolve CenQue class object by one time step 
    """
    evo_start = time.time()

    mass_bins = child_cq.mass_bins
    
    # remove galaxies below the min and max mass
    within_massbin = np.where(
            (child_cq.mass > min(mass_bins.mass_low)) & 
            (child_cq.mass <= max(mass_bins.mass_high))
            ) 
    child_cq.sample_trim(within_massbin)                   

    n_child = len(within_massbin[0])     
    n_parent = len(parent_cq.mass)

    # SF and stellar mass properties are added as attributes 
    # to the CenQue object
    child_cq.gal_type         = np.array(['' for i in xrange(n_child)], dtype='|S16') 
    child_cq.sfr              = np.array([-999. for i in xrange(n_child)]) 
    child_cq.ssfr             = np.array([-999. for i in xrange(n_child)]) 
    child_cq.q_ssfr           = np.array([-999. for i in xrange(n_child)]) 
    child_cq.tau              = np.array([-999. for i in xrange(n_child)]) 
    child_cq.parent_sfr       = np.array([-999. for i in xrange(n_child)]) 
    child_cq.parent_mass      = np.array([-999. for i in xrange(n_child)]) 
    child_cq.parent_halo_mass = np.array([-999. for i in xrange(n_child)]) 

    for attrib in ['gal_type', 'sfr', 'ssfr', 'q_ssfr', 'tau', 
            'parent_sfr', 'parent_mass', 'parent_halo_mass']: 
        if attrib not in child_cq.data_columns: 
            child_cq.data_columns.append(attrib)

    if child_cq.sf_prop['name'] == 'average': 
        child_cq.delta_sfr = np.array([-999. for i in range(n_child)]) 
        extra_attr = ['delta_sfr']
    else: 
        raise NotImplementedError() 

    for attrib in extra_attr: 
        if attrib not in child_cq.data_columns: 
            child_cq.data_columns.append(attrib)

    # parent children match which assigns indices into dictionaries and then get dictionary values
    parents, children = parent_children_match(parent_cq.snap_index, child_cq.parent)

    # Deal with parent to children inheritance (children inherit the parents' attributes)
    # This step needs to be thought out more. 
    inheritance_time = time.time()
    child_cq.parent_sfr[children]       = parent_cq.sfr[parents]
    child_cq.parent_mass[children]      = parent_cq.mass[parents]
    child_cq.parent_halo_mass[children] = parent_cq.halo_mass[parents]
    
    for attr in child_cq.data_columns: 
        # Inherit everything but the TreePM and Halo information
        if attr not in ('mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index', 'pos', 
                'parent_sfr', 'parent_mass', 'parent_halo_mass'): 
            try: 
                parent_attr = getattr(parent_cq, attr)[parents]
                getattr(child_cq, attr)[children] = parent_attr
            except AttributeError: 
                if 'evol_from' not in parent_cq.cenque_type:
                    pass 
                else: 
                    raise ValueError()
    print 'Inheritance takes ', time.time() - inheritance_time
    
    # Stellar mass evolution of star forming galaxies. Either integrated or 
    # SHAM masses assigned from TreePM 
    if child_cq.mass_evol == 'integrated':
        # integrated SFR mass evolution of star-forming galaxies 
        # (only star-forming galaxies) 

        if kwargs['sfr'] == 'sfr_func': 

            integrated_mass = util.integrated_mass_rk4(
                    util.sfr_squarewave, (child_cq.parent_mass)[sf_child_indx], 
                    parent_cq.t_cosmic, child_cq.t_cosmic, 
                    amp = (child_cq.sfr_amp)[sf_child_indx], 
                    freq = (child_cq.sfr_freq)[sf_child_indx], 
                    phase = (child_cq.sfr_phase)[sf_child_indx])
            (child_cq.mass)[sf_child_indx] = integrated_mass

        elif kwargs['sfr'] == 'sfr_avg': 

            (child_cq.mass)[sf_child_indx] = util.integrated_mass_rk4(
                    util.sfr_avg_residual, (child_cq.parent_mass)[sf_child_indx], 
                    parent_cq.t_cosmic, child_cq.t_cosmic, 
                    resid = (child_cq.sfr_resid)[sf_child_indx])
            if not silent: 
                print 'Integrated Mass vs SHAM mass' 
                print (child_cq.mass)[sf_child_indx] - (child_cq.sham_mass)[sf_child_indx]

    # Evolve Quiescent Children. 
    # Quiecsent galaxies keep the same sSFR while the quiescent galaxy 
    # masses are assigned by SHAM masses
    quiescent_children = np.where(
            (child_cq.gal_type[children] == 'quiescent') & 
            (child_cq.tau[children] == -999.0)
            )
    q_children = children[quiescent_children]
    
    # correct the SFR in order to preserve SSFR
    # this is because the mass evolves through SHAM  
    child_cq.sfr[q_children] = child_cq.ssfr[q_children] + child_cq.mass[q_children]
    
    # Evolve Star-forming Children 
    # include doucmentation details 
    # include doucmentation details 
    # include doucmentation details 
    # include doucmentation details 
    # include doucmentation details 
    starforming_children = np.where(
            child_cq.gal_type[children] == 'star-forming'
            )
    sf_children = children[starforming_children] 

    if child_cq.sf_prop['name'] == 'average': 
        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()
        child_sfr = sfr_mstar_z(
                child_cq.mass[sf_children], 
                child_cq.zsnap
                )

        child_cq.sfr[sf_children] = child_sfr + child_cq.delta_sfr[sf_children]
    else: 
        raise NotImplementedError() 

    child_cq.ssfr[sf_children] = child_cq.sfr[sf_children] - child_cq.mass[sf_children]

    # Evolve quenching children by exp(- delta t_cosmic / tau)
    quenching_time = time.time()
    still_quenching = np.where(child_cq.tau > 0.0)

    if len(still_quenching[0]) > 0: 
        # for galaxies with tau, keep them quenching!
        tau_quench = np.log10( 
                np.exp( 
                    -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[still_quenching]
                    )
                )                                              # SFR quenching amount  

        child_cq.sfr[still_quenching] = child_cq.parent_sfr[still_quenching] + tau_quench 
        child_cq.ssfr[still_quenching] = child_cq.sfr[still_quenching] - child_cq.mass[still_quenching]
    else: 
        if 'evol_from' in parent_cq.cenque_type:
            raise ValueError()

    print 'Quenching SF galaxies takes ', time.time() - quenching_time
    
    # numpy.where of mass bins 
    mass_bin_indices = np.array([
        np.where(
            (child_cq.mass > mass_bins.mass_low[i_m]) & 
            (child_cq.mass <= mass_bins.mass_high[i_m]) & 
            (child_cq.gal_type != '') 
            )
        for i_m in xrange(mass_bins.nbins)
        ])
    
    # number of galaxies in each of the mass bins 
    ngal_mass_bin = np.array([
        len(mass_bin_indices[i_m][0])
        for i_m in xrange(mass_bins.nbins)
        ])
        
    # f_Q of child CenQue as a function of mass
    # f_Q(M_bin)
    child_fq = get_fq(
            mass_bins.mass_mid, 
            child_cq.zsnap, 
            lit = child_cq.fq_prop['name']
            ) 
    
    # f_Q of descendant 
    descendant_fq = get_fq_nsnap(
            mass_bins.mass_mid, 
            child_cq.nsnap - 2, 
            lit = child_cq.fq_prop['name']
            )
    descendant_tcosmic = util.get_t_nsnap(child_cq.nsnap - 2)
    
    # expected number of quiescent galaxies = Ngal * f_Q
    # N_gal,exp * f_Q(M_bin)
    exp_ngal_q = np.rint( child_fq * ngal_mass_bin ).astype(int)

    child_galtype_list = [
            child_cq.gal_type[mass_bin_indices[i_m][0]]
            for i_m in xrange(mass_bins.nbins)
            ]
    sf_massbin_indices = [ 
            (mass_bin_indices[i_m][0])[np.where(child_galtype_list[i_m] == 'star-forming')]
            for i_m in xrange(mass_bins.nbins)
            ]
    ngal_sf_mass_bin = [ 
            len(sf_massbin_indices[i_m]) 
            for i_m in xrange(mass_bins.nbins)
            ]
    ngal2quench = exp_ngal_q - ngal_mass_bin + ngal_sf_mass_bin 
    
    exp_exp_ngal_q = np.rint( descendant_fq * ngal_mass_bin ).astype(int)
    ngal2quenchfuture = exp_exp_ngal_q - ngal_mass_bin + ngal_sf_mass_bin 

    quench_index = [] 
    quenching_start_time = time.time()
    # Quench a number of star-forming galaxies in order to match the 
    # quiescent fraction at the next cosmic time step 

    for i_m in xrange(mass_bins.nbins):              

        quench_index += quenching_galaxies_massbin(
                child_cq, 
                child_cq.t_cosmic - parent_cq.t_cosmic, 
                descendant_tcosmic - parent_cq.t_cosmic, 
                mass_bins.mass_mid[i_m], 
                mass_bin_indices[i_m], 
                sf_massbin_indices[i_m], 
                ngal2quench[i_m], 
                ngal2quenchfuture[i_m], 
                quiet=quiet
                )

    quench_index = np.array(quench_index)

    if len(quench_index) != 0:

        child_cq.gal_type[quench_index] = 'quiescent'  # boom quenched 

        # assign quenching e-folds for quenching galaxies
        child_cq.tau[quench_index] = get_quenching_efold(
                child_cq.parent_mass[quench_index], 
                tau_param = child_cq.tau_prop
                )

        tau_quench = np.log10( np.exp( 
            -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[quench_index]
            )) 

        child_cq.sfr[quench_index] = child_cq.parent_sfr[quench_index] + tau_quench 
        child_cq.ssfr[quench_index] = child_cq.sfr[quench_index] - child_cq.mass[quench_index]

        q_ssfr_mean = util.get_q_ssfr_mean(child_cq.mass[quench_index]) 
        child_cq.q_ssfr[quench_index] = 0.18 * np.random.randn(len(quench_index)) + q_ssfr_mean 
        
    print 'Quenching takes ', time.time() - quenching_start_time 
    
    if not quiet: 
        print "Child's redshift = ", child_cq.zsnap

    # deal with orphans ------------------------------------------------------------
    orphan_time = time.time()
    print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
    child_cq = assign_sfr(child_cq, quiet=quiet)

    if len(child_cq.gal_type[child_cq.gal_type == '']) != 0:
        print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
        raise ValueError() 
    print 'Orphan care takes ', time.time() - orphan_time

    # deal with galaxies that are being quenched beyond the designated final quenching  
    overquench_time = time.time()
    over_quenched = np.where(child_cq.ssfr < child_cq.q_ssfr)
    child_cq.ssfr[over_quenched] = child_cq.q_ssfr[over_quenched]
    child_cq.tau[over_quenched] = -999.0            # done quenching 
    print 'Overquenching takes ', time.time() - overquench_time

    print 'Evolution for one time step takes ', time.time() - evo_start

    if not quiet: 
        print 'Quiescent Fraction = ', np.float(len(parent_cq.gal_type[parent_cq.gal_type == 'quiescent']))/np.float(len(parent_cq.gal_type)) 

    return child_cq

def quenching_galaxies_massbin(cenque, delta_t_cosmic, descendant_t, m_bin_mid, indices, sf_indices, ngal2quench, ngal2quench_future, quiet=False):
    """ Select which star-forming galaxies to quench in a specified mass bin 
    """

    ngal = len(indices[0])
    ngal_sf = len(sf_indices)

    if ngal_sf > ngal: 
        raise ValueError()

    if (ngal == 0) or (ngal_sf == 0):
        return [] 

    # number of SF galaxies that need to be quenched 
    if ngal2quench <= 0: 
        # quenching, overshot in previous snapshot  move onto next mass bin 
        #if not quiet: 
        print "No SF galaxies need to be quenched"
        return [] 

    if not quiet: 
        print '--------------------------------------------------------------------------------'
        print 'Ngal,sf = ', ngal_sf, ' Ngal,q = ', ngal - ngal_sf
        print 'current fq = ', 1.0 - np.float(ngal_sf)/np.float(ngal)
        print ngal2quench, ' SF galaxies need to be quenched'

    # number of SF galaxies to quench will be determined by first determining how many 
    # galaxies will be quenched if the entire SF population was quenched, taking the 
    # ratio of that number over the 'exp_ngal_q_massbin'. We use that factor as the 
    # "quenching fraction" (the fraction of SF galaxies that will be quenched for this mass
    # bin and redshift) 
    pred_taus = get_quenching_efold(
            cenque.parent_mass[sf_indices], 
            tau_param = cenque.tau_prop
            )
    
    pred_tau_quench = np.log10(np.exp( 
        -1.0 * delta_t_cosmic / pred_taus 
        ))      

    pred_sfqs = sfq_classify(
            cenque.mass[sf_indices], 
            cenque.sfr[sf_indices] + pred_tau_quench, 
            cenque.zsnap
            )
    
    pred_ngal_q = np.float(np.sum(pred_sfqs == 'quiescent'))
    
    predpred_tau_quench = np.log10(np.exp( 
        -1.0 * descendant_t / pred_taus 
        ))      

    predpred_sfqs = sfq_classify(
            cenque.mass[sf_indices], 
            cenque.sfr[sf_indices] + predpred_tau_quench, 
            util.get_z_nsnap(cenque.nsnap - 2) 
            )
    predpred_ngal_q = np.float(np.sum(predpred_sfqs == 'quiescent'))

    if pred_ngal_q == 0: 
        warnings.warn('Quenching time-scale is too long')
        #raise ValueError("What the fuck") 
        
    # Quenching fraction evalutation
    #alpha = 1.5
    #fqing1 = 0.025 * (m_bin_mid - 9.5) * (1.8 - cenque.zsnap)**alpha + 0.15

    if ngal2quench == ngal_sf: 
        pred_fqing = 1.0
        predpred_fqing = 1.0
    else: 
        try: 
            pred_fqing = np.float(ngal2quench)/np.float(pred_ngal_q)
        except ZeroDivisionError: 
            pred_fqing = 1.0

        try: 
            predpred_fqing = np.float(ngal2quench_future)/np.float(predpred_ngal_q)
        except ZeroDivisionError: 
            predpred_fqing = 1.0

    f_quenching = min(pred_fqing, predpred_fqing)

    if f_quenching <= 0.0: 
        f_quenching = 0.0
        return []
    elif f_quenching > 1.0: 
        f_quenching = 1.0 
    
    quench_index = random.sample( 
            sf_indices, 
            int(np.rint(f_quenching * np.float(ngal_sf))) 
            )
    return quench_index

def parent_children_match(parent_snap_index, child_parent_snap_index):
    """ Match snapshot index of parents to the parent snapshot index of
    childrens by matching (index, parent), (index, child) key pairs 
    Takes approximately < 1 second
    """
    parent_child_start_time = time.time()
    parent_index_dict = dict( 
            (k, i) 
            for i, k in enumerate(parent_snap_index) 
            )
    child_index_dict = dict( 
            (k, i) 
            for i, k in enumerate(child_parent_snap_index) 
            )

    parent_child_meet = np.intersect1d( parent_snap_index, child_parent_snap_index ) 
    parents = np.array([
            parent_index_dict[k] for k in parent_child_meet
            ])
    children = np.array([
            child_index_dict[k] for k in parent_child_meet
            ])
    print 'Parent Child match takes ', time.time() - parent_child_start_time, ' seconds'

    return parents, children 

if __name__=='__main__': 
    blah = CenQue(n_snap=13, cenque_type='sf_assigned')
    #blah.import_treepm(20)
    #blah.writeout()
    #blah = assign_sfr(blah)
    #blah.writeout()
    blah.readin()
    blah = evolve_cq(blah, tau_prop = {'name': 'satellite'}, quiet=True)
