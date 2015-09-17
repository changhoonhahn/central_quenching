"""

Evolve CenQue object over cosmic time 

Author(s): ChangHoon Hahn

"""
import time
import numpy as np

from cenque import CenQue
from quiescent_fraction import get_fq

def evolve_CenQue(
        cenque, 
        final_nsnap = 1, 
        sf_prop = {'name': 'average'}, 
        fq_prop = {'name': 'wetzelsmooth'}, 
        tau_prop = {'name': 'instant'}, 
        mass_evol = 'sham', 
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
    
    if cenque.cenque_type != 'sf_assigned': 
        raise ValueError()

    if not quiet: 
        print 'Quiescent Fraction', 
        print np.float(len(cenque.gal_type[cenque.gal_type == 'quiescent']))/np.float(len(cenque.gal_type)) 
    
    if kw_sfr == 'average': 
        sfr_mstar_z, sig_sfr_mstar_z = get_bestfit_sfr_mstar_z(
                Mrcut = 18, 
                fid_mass = 10.5
                )
    else: 
        raise NotImplementedError() 

    start_nsnap = cenque.nsnap 
    parent_cq = cenque
    
    for i_snap in xrange(start_nsnap, final_nsnap):    

        child_cq = CenQue() 
        child_cq.nsnap = i_snap - 1  
        child_cq.import_treepm(child_cq.nsnap) 
        child_cq.cenque_type = ''.join(['evol_from', str(start_nsnap)])
        child_cq.sf_prop = sf_prop
        child_cq.fq_prop = fq_prop
        child_cq.tau_prop = tau_prop
        child_cq.mass_evol = mass_evol

        mass_bins = child_cq.mass_bins
        
        # remove galaxies below the min and max mass
        within_massbin = np.where(
                (child_cq.mass > min(mass_bins.mass_low)) & 
                (child_cq.mass <= max(mass_bins.mass_high)) &
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
            if attrib not in cenque.data_columns: 
                child_cq.data_columns.append(attrib)
            else: 
                raise ValueError()
    
        if kw_sfr == 'average': 
            child_cq.delta_sfr = np.array([-999. for i in range(n_child)]) 
            extra_attr = ['delta_sfr']
        else: 
            raise NotImplementedError() 

        for attrib in extra_attr: 
            if attrib not in cenque.data_columns: 
                cenque.data_columns.append(attrib)

        # parent children match
        # assign indices into dictionaries and then get dictionary values
        # not very good -----------------------------------------------------------------
        parent_child_start_time = time.time()
        parent_index_dict = dict( 
                (k, i) 
                for i, k in enumerate(parent_cq.snap_index) 
                )
        child_index_dict = dict( 
                (k, i) 
                for i, k in enumerate(child_cq.parent) 
                )

        parent_child_meet = np.intersect1d( parent_cq.snap_index, child_cq.parent ) 
        parents = np.array([
                parent_index_dict[k] for k in parent_child_meet
                ])
        children = np.array([
                child_index_dict[k] for k in parent_child_meet
                ])
        print 'Parent Child match takes ', time.time() - parent_child_start_time, ' seconds'
        # --------------------------------------------------------------------------------
        # Deal with parent to children inheritance (children inherit the parents' attributes)
        child_cq.gal_type[children]         = parent_cq.gal_type[parents]
        child_cq.parent_sfr[children]       = parent_cq.sfr[parents]
        child_cq.parent_mass[children]      = parent_cq.mass[parents]
        child_cq.parent_halo_mass[children] = parent_cq.halo_mass[parents]
        
        for attr in extra_attr: 
            parent_attr = getattr(parent_cq, attr)[parents]
            setattr(child_cq, attr, parent_attr)
        
        # Evolve Quiescent Children -------------------------------------------------------
        quiescent_children = np.where(
                (child_cq.gal_type[children] == 'quiescent') & 
                (child_cq.tau[children] == -999.0)
                )
        q_children = (children[0])[quiescent_children]
        q_children_parents = (parents[0])[quiescent_children]

        # keep SSFR same 
        child_cq.ssfr[q_children] = parent_cq.ssfr[q_children_parents]
        child_cq.sfr[q_children] = child_cq.ssfr[q_children] + child_cq.mass[q_children]

        # --------------------------------------------------------------------------------
        # Evolve star-forming children using the specified sfr perscription 
        starforming_children = np.where(
                child_cq.gal_type[children] == 'star-forming'
                )
        sf_children = (children[0])[starforming_children] 
        sf_children_parents = (parents[0])[starforming_children]
    
        if kw_sfr == 'average': 
            child_sfr, child_sig_sfr = sfr_mstar_z(
                    child_cq.mass[sf_children], 
                    child_cq.zsnap
                    )

            child_cq.sfr[sf_children] = child_sfr + child_cq.delta_sfr[sf_children]
        
        else: 
            #elif kwargs['sfr'] == 'sfr_func': 

            #    child_cq.sfr[sf_child_indx] = util.sfr_squarewave(
            #            child_cq.mass[sf_child_indx], child_cq.t_cosmic, 
            #            amp = child_cq.sfr_amp[sf_child_indx], 
            #            freq = child_cq.sfr_freq[sf_child_indx], 
            #            phase = child_cq.sfr_phase[sf_child_indx]) 
            raise NotImplementedError() 

        child_cq.ssfr[sf_children] = child_cq.sfr[sf_children] - child_cq.mass[sf_children]

        # --------------------------------------------------------------------------------
        # Stellar mass evolution of star forming galaxies 
        if mass_evol == 'integrated':
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

        # Evolve quenching children by exp(- delta t_cosmic / tau)
        try: 
            # children inherit parent tau 
            child_cq.tau[children] = parent_cq.tau[parents]
            has_tau = np.where(
                    child_cq.tau[children] > 0.0
                    )
            still_quenching = (children[0])[has_tau]
            still_quenching_parent = (parents[0])[has_tau]
            
            # for galaxies with tau, keep them quenching!
            tau_quench = np.log10( 
                    np.exp( 
                        -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[still_quenching]
                        )
                    )                                              # SFR quenching amount  

            child_cq.sfr[still_quenching] = child_cq.parent_sfr[still_quenching] + tau_quench 
            child_cq.ssfr[still_quenching] = child_cq.sfr[still_quenching] - child_cq.mass[still_quenching]
            
            # children inherit final quenched SSFR 
            child_cq.q_ssfr[children] = parent_cq.q_ssfr[parents]

        except AttributeError:
            pass

        quenching_fractions = []
        quenching_fractionss = [] 

        # Quench a number of star-forming galaxies in order to match the 
        # quiescent fraction at the next cosmic time step 
        for i_m in xrange(mass_bins.nbins):              

            mass_bin = np.where(
                    (child_cq.mass > mass_bins.mass_low[i_m]) & 
                    (child_cq.mass <= mass_bins.mass_high[i_m]) & 
                    (child_cq.gal_type != '') 
                    )
            ngal_massbin = len(mass_bin[0])

            if not quiet: 
                print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]

            if ngal_massbin == 0:              
                print 'No Galaxies in mass bin'
                continue 
            
            fq_massbin = get_fq(
                    mass_bins.mass_mid[i_m], 
                    child_cq.zsnap, 
                    lit = fq_prop['name']
                    ) 
            
            if not quiet: 
                print 'nsnap = ', child_cq.nsnap, ' z = ', child_cq.zsnap
                print ' M* = ', mass_bins.mass_mid[i_m], ' fq = ', fq_massbin 
                
            # Expected Ngal_q,mass
            exp_ngal_q_massbin = int(np.rint(
                        fq_massbin * np.float(ngal_massbin)
                        ))    

            if ngal_massbin < exp_ngal_q_massbin: 
                raise ValueError()
            
            child_gal_type = child_cq.gal_type[mass_bin] 
            ngal_sf_massbin = np.sum(child_gal_type == 'star-forming') 
            ngal_q_massbin  = np.sum(child_gal_type == 'quiescent') 
            if ngal_sf_massbin == 0: 
                print "There are no SF galaxies to quench"
                continue

            # number of SF galaxies that need to be quenched 
            ngal_2quench = exp_ngal_q_massbin - ngal_q_massbin 
            if ngal_2quench <= 0: 
                # quenching, overshot in previous snapshot  move onto next mass bin 
                if not quiet: 
                    print 'Quenching overshot in previous snapshot' 
                continue 

            if not quiet: 
                print 'Ngal,sf = ', ngal_sf_massbin, ' Ngal,q = ', ngal_q_massbin 
                print exp_ngal_q_massbin, ' expected quiescent galaxies' 
                print 'current fq = ', np.float(ngal_q_massbin)/np.float(ngal_massbin)
                print ngal_2quench, ' SF galaxies need to be quenched'

            # Star-forming children 
            sf_massbin = (mass_bin[0])[
                    np.where(child_gal_type == 'star-forming')
                    ] 
                
            # number of SF galaxies to quench will be determined by first determining how many galaxies will be quenched
            # if the entire SF population was quenched, taking the ratio of that number over the 'exp_ngal_q_massbin'. 
            # we use that factor as the "quenching fraction" (the fraction of SF galaxies that will be quenched for this mass
            # bin and redshift) 
            pred_taus = util.get_quenching_efold(
                    child_cq.parent_mass[sf_massbin], 
                    tau_param = tau_prop
                    )
            
            pred_tau_quench = np.log10(np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / pred_taus 
                ))      

            pred_sfqs = util.sfq_classify(
                    child_cq.mass[sf_massbin], 
                    child_cq.sfr[sf_massbin] + tau_quench, 
                    child_cq.zsnap
                    )
            pred_ngal_q = np.float(np.sum(pred_sfqs == 'quiescent'))
    
            if pred_ngal_q == 0: 
                raise ValueError("What the fuck") 

            alpha = 1.5
            quenching_fraction = 0.025 * (mass_bins.mass_mid[i_m] - 9.5) * (1.8 - child_cq.zsnap)**alpha + 0.15
            quenching_fractions.append(quenching_fraction)

            quenching_fraction = np.float(ngal_2quench)/np.float(ngal_totalq)
            quenching_fractionss.append(quenching_fraction)
            
            alpha = 1.5
            fqing1 = 0.025 * (mass_bins.mass_mid[i_m] - 9.5) * (1.8 - child_cq.zsnap)**alpha + 0.15
                
            fqing2 = np.float(ngal_2quench)/np.float(ngal_totalq)

            if (mass_bins.mass_mid[i_m] < 11.0) and (fqing2 > fqing1): 
                quenching_fraction = fqing1
                if not silent: 
                    print '####################################### FQING 2 > FQING 1'
            else: 
                quenching_fraction = fqing2

            if quenching_fraction < 0.0: 
                quenching_fraction = 0.0
            elif quenching_fraction > 1.0: 
                quenching_fraction = 1.0
            #print 'quenching fraction ', quenching_fraction
            if quenching_fraction == 0.0: 
                continue

            quench_index = random.sample(mbin_sf_index, int( np.rint(quenching_fraction * np.float(mbin_sf_ngal) ) ) )
            
        child_cq.gal_type[quench_index] = 'quiescent'  # boom quenched 
        
        # assign quenching e-folds for quenching galaxies
        child_cq.tau[quench_index] = util.get_quenching_efold(
                child_cq.parent_mass[sf_massbin], 
                tau_param = tau_prop
                )
        tau_quench = np.log10( np.exp( 
            -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[quench_index]
            ))      

        child_cq.sfr[quench_index] = child_cq.parent_sfr[quench_index] + tau_quench 
        child_cq.ssfr[quench_index] = child_cq.sfr[quench_index] - child_cq.mass[quench_index]

        q_ssfr_mean = util.get_q_ssfr_mean(child_cq.mass[quench_index]) 
        child_cq.q_ssfr[quench_index] = 0.18 * np.random.randn(len(quench_index)) + q_ssfr_mean 
        
        if not quiet: 
            print child_cq.zsnap
            print mass_bins.mass_mid
            print quenching_fractions
            print quenching_fractionss

        # deal with orphans ------------------------------------------------------------

        sfing = np.where(child_cq.gal_type == 'star-forming') 
        print max(child_cq.mass[sfing]) 
        print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
        child_cq.AssignSFR(child_cq.nsnap, **kwargs) 

        if len(child_cq.gal_type[child_cq.gal_type == '']) != 0:
            print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
            raise NameError('asdflkjasdlfj') 
        
        # deal with galaxies that are being quenched beyond the designated final quenching  
        over_quenched = child_cq.ssfr < child_cq.q_ssfr 
        child_cq.ssfr[over_quenched] = child_cq.q_ssfr[over_quenched]
        child_cq.tau[over_quenched] = -999.0            # done quenching 

        #if child_cq.nsnap == final_nsnap: 
    
        if kwargs['stellmass'] == 'sham': 
            if kwargs['sfr'] == 'sfr_avg': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_resid', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                        **kwargs)  
            elif kwargs['sfr'] == 'sfr_func': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_amp', 'sfr_freq', 'sfr_phase', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                    **kwargs)  

        elif kwargs['stellmass'] == 'integrated': 
            if kwargs['sfr'] == 'sfr_avg': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_resid', 'sham_mass', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                        **kwargs)  
            elif kwargs['sfr'] == 'sfr_func': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_amp', 'sfr_freq', 'sfr_phase', 'sham_mass', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                    **kwargs)  
        else: 
            raise NotImplementedError('asdflkjasdkf') 

        parent_cq = child_cq

        if not silent: 
            print 'Quiescent Fraction = ', np.float(len(parent_cq.gal_type[parent_cq.gal_type == 'quiescent']))/np.float(len(parent_cq.gal_type)) 

