"""

Evolve CenQue object over cosmic time 

Author(s): ChangHoon Hahn

"""

def EvolveCenQue(origin_nsnap, final_nsnap, mass_bin=None, silent=True, **kwargs): 
    ''' Evolve SF properties from origin_nsnap to final_nsnap 
    
    Parameters
    ----------

    Notes
    -----
    * Quenching of SF galaxies now involve quenching-efold *and* overall SF MS SFR decrease

    '''

    if mass_bin is None:            # mass bins
        mass_bins = util.simple_mass_bin()  # use simplest mass bins
    else: 
        raise NotImplementedError("not yet coded") 

    # SF-MS fits 
    if 'sfms_slope' in kwargs.keys(): 
        sfms_sfr_fit = sfms.get_sfmsfit_sfr(kwargs['sfms_slope'], kwargs['sfms_yint'])
    else: 
        pass
        #groupcat_slope, groupcat_yint = sfms.get_bestfit_groupcat_sfms(Mrcut=18)
        #sfms_sfr_fit = sfms.get_sfmsfit_sfr(groupcat_slope, groupcat_yint)
  
    # import original snap SF prop 
    parent_cq = CenQue()
    parent_cq.readin(nsnap=origin_nsnap, file_type='sf assign', **kwargs)   
    
    if not silent: 
        print 'Quiescent Fraction = ', np.float(len(parent_cq.gal_type[parent_cq.gal_type == 'quiescent']))/np.float(len(parent_cq.gal_type)) 
    
    # evolve snapshot by snapshot ------------------------------------------------------
    for i_step in range(0, origin_nsnap - final_nsnap):    

        i_snap = origin_nsnap - i_step  # current Nsnap 
        child_snap = i_snap-1

        child_cq = CenQue() 
        child_cq.readin(nsnap=child_snap)         # read in snapshot SHAM masses
        
        # remove galaxies beyond mass bins
        child_mass_limit = (child_cq.mass > min(mass_bins.mass_low)) & \
                (child_cq.mass <= max(mass_bins.mass_high))
        child_cq.sample_select( child_mass_limit, 
                columns = ['mass', 'sham_mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index']  
                )                   
        n_child = len(child_cq.mass)     # number of children left 

        # set up columns for assignment 
        child_cq.gal_type = np.array(['' for i in range(n_child)], dtype='|S16') 
        child_cq.sfr = np.array([-999. for i in range(n_child)]) 
        child_cq.ssfr = np.array([-999. for i in range(n_child)]) 
        child_cq.q_ssfr = np.array([-999. for i in range(n_child)]) 
        child_cq.tau = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_sfr = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_mass = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_halo_mass = np.array([-999. for i in range(n_child)]) 

        if kwargs['sfr'] == 'sfr_func': 
            child_cq.sfr_amp = np.array([-999. for i in range(n_child)]) 
            child_cq.sfr_freq = np.array([-999. for i in range(n_child)]) 
            child_cq.sfr_phase = np.array([-999. for i in range(n_child)]) 

        elif kwargs['sfr'] == 'sfr_avg': 
            child_cq.sfr_resid = np.array([-999. for i in range(n_child)]) 

        else: 
            raise NotImplementedError('asdlkfj') 

        # parent children match
        # assign indices into dictionaries and then get dictionary values
        # not very good -----------------------------------------------------------------
        parent_index_dict = dict( (k, i) for i, k in enumerate(parent_cq.snap_index) )
        child_index_dict = dict( (k, i) for i, k in enumerate(child_cq.parent) )

        parent_child_intersect = np.intersect1d( parent_cq.snap_index, child_cq.parent ) 
        parent_indx = [ parent_index_dict[x] for x in parent_child_intersect ] 
        child_indx = [ child_index_dict[x] for x in parent_child_intersect ]  
        # --------------------------------------------------------------------------------
        
        # children inherit galaxy type and store parent SFR and mass 
        (child_cq.gal_type)[child_indx] = [(parent_cq.gal_type)[i] for i in parent_indx]
        (child_cq.parent_sfr)[child_indx] = [(parent_cq.sfr)[i] for i in parent_indx]
        (child_cq.parent_mass)[child_indx] = [(parent_cq.mass)[i] for i in parent_indx]
        (child_cq.parent_halo_mass)[child_indx] = [(parent_cq.halo_mass)[i] for i in parent_indx]
       
        if kwargs['sfr'] == 'sfr_func': 
            (child_cq.sfr_amp)[child_indx] = [(parent_cq.sfr_amp)[i] for i in parent_indx]
            (child_cq.sfr_freq)[child_indx] = [(parent_cq.sfr_freq)[i] for i in parent_indx]
            (child_cq.sfr_phase)[child_indx] = [(parent_cq.sfr_phase)[i] for i in parent_indx]

        elif kwargs['sfr'] == 'sfr_avg': 
            (child_cq.sfr_resid)[child_indx] = [(parent_cq.sfr_resid)[i] for i in parent_indx]

        else: 
            raise NotImplementedError('lkajsdfkj')
        
        # Star-forming Children ----------------------------------------
        sf_child = (child_cq.gal_type[child_indx] == 'star-forming')
        sf_child_indx = np.array(child_indx)[sf_child]                  # index stuff
        sf_child_parent_indx = np.array(parent_indx)[sf_child]
        
        # SFR evolution amount
        #child_sfr, child_sig_sfr = util.get_sfr_mstar_z(child_cq.mass[sf_child_indx], 
        #        child_cq.zsnap, lit='primusfit')
        #parent_sfr, parent_sig_sfr = util.get_sfr_mstar_z( parent_cq.mass[sf_child_parent_indx], 
        #        parent_cq.zsnap, lit='primusfit') 
        #child_sfr, child_sig_sfr = util.get_sfr_mstar_z_flex(child_cq.mass[sf_child_indx], 
        #        child_cq.zsnap, sfms_sfr_fit)
        #parent_sfr, parent_sig_sfr = util.get_sfr_mstar_z_flex(
        #        parent_cq.mass[sf_child_parent_indx], parent_cq.zsnap, sfms_sfr_fit) 
    
        if kwargs['sfr'] == 'sfr_avg': 
            child_sfr, child_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    child_cq.mass[sf_child_indx], child_cq.zsnap, Mrcut=18)
            parent_sfr, parent_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    parent_cq.mass[sf_child_parent_indx], parent_cq.zsnap, Mrcut=18) 
        
            ''' SFR evolution is assumed to be equal to the overall change in SFR  
            '''
            dSFRt = child_sfr - parent_sfr          # overall change in SFR
            
            # evolve sf children SFR
            #child_cq.sfr[sf_child_indx] = child_cq.parent_sfr[sf_child_indx] + dSFRt
            child_cq.sfr[sf_child_indx] = util.sfr_avg_residual(
                    child_cq.mass[sf_child_indx], child_cq.zsnap, 
                    resid = child_cq.sfr_resid[sf_child_indx] ) 
        
        elif kwargs['sfr'] == 'sfr_func': 

            child_cq.sfr[sf_child_indx] = util.sfr_squarewave(
                    child_cq.mass[sf_child_indx], child_cq.t_cosmic, 
                    amp = child_cq.sfr_amp[sf_child_indx], 
                    freq = child_cq.sfr_freq[sf_child_indx], 
                    phase = child_cq.sfr_phase[sf_child_indx]) 

        else: 
            raise NotImplementedError('asdflkjadf') 

        child_cq.ssfr[sf_child_indx] = child_cq.sfr[sf_child_indx] - child_cq.mass[sf_child_indx]
        
        # --------------------------------------------------------------------------------
        # mass evolution of star-forming galaxies (only star-forming galaxies) 
        # should quiescent galaxies remain the same mass or adopt SHAM masses? 

        if kwargs['stellmass'].lower() == 'integrated':

            # integrated stellar mass
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

        elif kwargs['stellmass'].lower() == 'sham': 
            # leave stellar mass 
            pass
        else: 
            raise NotImplementedError('asdfasdflkjasdflkj;') 

        # Quenching Children ------------------------------------------------        
        try: 
            parent_cq.tau
        except AttributeError:
            pass
        else: 
            # children inherit parent tau 
            (child_cq.tau)[child_indx] = [(parent_cq.tau)[i] for i in parent_indx]
            has_tau = (child_cq.tau)[child_indx] > 0.0  # has tau 
            still_quenching_indx = np.array(child_indx)[has_tau]
            still_quenching_parent_indx = np.array(parent_indx)[has_tau]
            
            # for galaxies with tau, keep them quenching!
            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[still_quenching_indx]
                ))                                              # SFR quenching amount  

            c_sfr, c_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    child_cq.mass[still_quenching_indx], child_cq.zsnap, Mrcut=18)
            p_sfr, p_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    child_cq.parent_mass[still_quenching_indx], parent_cq.zsnap, Mrcut=18) 

            child_cq.sfr[still_quenching_indx] = \
                    child_cq.parent_sfr[still_quenching_indx] + tau_quench 
            child_cq.ssfr[still_quenching_indx] = \
                    child_cq.sfr[still_quenching_indx] - child_cq.mass[still_quenching_indx]
            
            # children inherit final quenched SSFR 
            (child_cq.q_ssfr)[child_indx] = [(parent_cq.q_ssfr)[i] for i in parent_indx]

        # Quiescent Children --------------------------------------------
        q_child = (child_cq.gal_type[child_indx] == 'quiescent') & (child_cq.tau[child_indx] == -999.0) # not quenching
        q_child_indx = np.array(child_indx)[q_child]
        q_child_parent_indx = np.array(parent_indx)[q_child]

        # keep SSFR same 
        child_cq.ssfr[q_child_indx] = parent_cq.ssfr[q_child_parent_indx]
        child_cq.sfr[q_child_indx] = child_cq.ssfr[q_child_indx] + \
                child_cq.mass[q_child_indx]
        
        quenching_fractions = [] 
        quenching_fractionss = [] 

        for i_m in range(mass_bins.nbins):              
            # boolean list for mass range
            mass_bin_bool = (child_cq.mass > mass_bins.mass_low[i_m]) & \
                    (child_cq.mass <= mass_bins.mass_high[i_m]) & \
                    (child_cq.gal_type != '') 
            if not silent: 
                print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]

            # indices of galaxies within mass range
            mass_bin_index = child_cq.get_index(mass_bin_bool) 
        
            mbin_ngal = np.sum(mass_bin_bool)   # Ngal in mass bin 

            if mbin_ngal == 0:              
                # if there are no galaxies within mass bin
                print 'No Galaxies in ', mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]
                continue 
            
            # observed quiescent fraction at M* and z
            mbin_qf = util.get_fq(mass_bins.mass_mid[i_m], child_cq.zsnap, lit='wetzelsmooth') 
            
            if not silent: 
                print 'nsnap = ', child_cq.nsnap, ' z = ', child_cq.zsnap, ' M* = ', mass_bins.mass_mid[i_m], ' fq = ', mbin_qf
                
            # Number of expected quiescent galaxies (Ngal,Q_exp)
            mbin_exp_n_q = int( np.rint(mbin_qf * np.float(mbin_ngal)) )    

            if mbin_ngal < mbin_exp_n_q: 
                # if for some reason expected number of quiescent galaxies
                # exceed the number of galaxies in the bin, all of them are 
                # quiescent
                mbin_exp_n_q = mbin_ngal 
            
            child_gal_type = child_cq.gal_type[mass_bin_bool] 
            mbin_sf_ngal = np.sum(child_gal_type == 'star-forming') 
            mbin_q_ngal = np.sum(child_gal_type == 'quiescent')  

            # number of SF galaxies that need to be quenched 
            ngal_2quench = mbin_exp_n_q - mbin_q_ngal 

            if not silent: 
                print 'sf', mbin_sf_ngal, 'q', mbin_q_ngal
                print mbin_exp_n_q, ' expected quiescent galaxies' 
                print 'current fq = ', np.float(mbin_q_ngal)/np.float(mbin_ngal)
                print ngal_2quench, ' SF galaxies need to be quenched'

            if ngal_2quench <= 0: 
                # quenching, overshot in previous snapshot  
                # move onto next mass bin 
                if not silent: 
                    print 'Quenching overshot in previous snapshot' 
                pass
                #continue 

            if mbin_sf_ngal == 0: 
                continue
            
            # Star-forming children 
            mbin_sf_index = child_cq.get_index( 
                    [ mass_bin_bool & (child_cq.gal_type == 'star-forming')] 
                    ) 
                
            # number of SF galaxies to quench will be determined by first determining how many galaxies will be quenched
            # if the entire SF population was quenched, taking the ratio of that number over the 'mbin_exp_n_q'. 
            # we use that factor as the "quenching fraction" (the fraction of SF galaxies that will be quenched for this mass
            # bin and redshift) 
            if 'tau' in kwargs.keys(): 
                if 'tau_param' in kwargs.keys(): 
                    taus = util.get_quenching_efold(
                            child_cq.parent_mass[mbin_sf_index], 
                            type=kwargs['tau'], param=kwargs['tau_param'])
                else: 
                    taus = util.get_quenching_efold(
                            child_cq.parent_mass[mbin_sf_index], 
                            type=kwargs['tau'])
            else: 
                raise TypeError('specify quenching e-fold: tau = instant, constant, linear') 
            
            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / taus 
                ))      
            #print child_cq.parent_sfr[mbin_sf_index] + tau_quench - child_cq.mass[mbin_sf_index]
            #print tau_quench
            #print min(child_cq.sfr[mbin_sf_index] + tau_quench - child_cq.mass[mbin_sf_index])
            #print max(child_cq.sfr[mbin_sf_index] + tau_quench - child_cq.mass[mbin_sf_index])

            sfqs = util.sfq_classify(child_cq.mass[mbin_sf_index], child_cq.sfr[mbin_sf_index] + tau_quench, child_cq.zsnap)
            ngal_totalq = np.float(np.sum(sfqs == 'quiescent'))
            #print min(child_cq.sfr[mbin_sf_index[sfqs == 'quiescent']] + tau_quench[sfqs == 'quiescent'] - child_cq.mass[mbin_sf_index[sfqs == 'quiescent']])
            #print max(child_cq.sfr[mbin_sf_index[sfqs == 'quiescent']] + tau_quench[sfqs == 'quiescent'] - child_cq.mass[mbin_sf_index[sfqs == 'quiescent']])
    
            if ngal_totalq == 0: 
                raise NameError("What the fuck") 
            #print 'ngal_totalq', ngal_totalq

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

            #if quenching_fraction > 1.0: 
            #    #raise NameError('asdfasdfasdf')
            #    #quench_index = mbin_sf_index
            #    print 'QUENCHING FRACTION TOO LARGE ******************************************************************'
            #    quench_index = random.sample(mbin_sf_index, ngal_2quench) 
            #else: 
            quench_index = random.sample(mbin_sf_index, int( np.rint(quenching_fraction * np.float(mbin_sf_ngal) ) ) )
            
            #quench_index = random.sample(mbin_sf_index, ngal_2quench) 
            
            child_cq.gal_type[quench_index] = 'quiescent'  # boom quenched 
            
            # assign quenching e-folds for quenching galaxies
            if 'tau' in kwargs.keys(): 
                if 'tau_param' in kwargs.keys(): 
                    child_cq.tau[quench_index] = util.get_quenching_efold(
                            child_cq.parent_mass[quench_index], 
                            type=kwargs['tau'], param=kwargs['tau_param'])
                else: 
                    child_cq.tau[quench_index] = util.get_quenching_efold(
                            child_cq.parent_mass[quench_index], 
                            type=kwargs['tau'])
            else: 
                raise TypeError('specify quenching e-fold: tau = instant, constant, linear') 
    
            # SFR quenching amount
            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[quench_index]
                ))      

            child_cq.sfr[quench_index] = child_cq.parent_sfr[quench_index] + tau_quench 
            child_cq.ssfr[quench_index] = child_cq.sfr[quench_index] - child_cq.mass[quench_index]

            q_ssfr_mean = util.get_q_ssfr_mean(child_cq.mass[quench_index]) 
            child_cq.q_ssfr[quench_index] = 0.18 * np.random.randn(len(quench_index)) + q_ssfr_mean 
        
        if not silent: 
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

