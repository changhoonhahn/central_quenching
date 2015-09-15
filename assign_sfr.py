"""

Assign StarFomation properties to CenQue object

Author(s): ChangHoon Hahn

"""

import sf_mainseq as sfms

def assign_sfr(cenque, quiet=True, kw_sfr='average', **kwargs):
    """ Assign star-formation properties to CenQue object
    
    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    cenque : CenQue class object 
    ----------------------------------------------------------------
    Notes 
    ----------------------------------------------------------------
    * Re imagine how to set up SFR 
    """

    start_time = time.time()

    if cenque.mass == None: 
        raise ValueError()

    if cenque.zsnap == None: 
        raise ValueError()
    
    mass_bins = cenque.mass_bins
    mass_bin_low = round(mass_bins.mass_low, 2) 
    mass_bin_mid = round(mass_bins.mass_mid, 2) 
    mass_bin_high = round(mass_bins.mass_high, 2) 

    # remove galaxies below the minimum and maximum mass and galaxies without children 
    if (min(cenque.mass) < mass_bins.mass_low.min()) or (max(cenque.mass) > mass_bins.mass_high.max()): 
        
        within_massbin_with_child = np.where(
                (cenque.mass > mass_bins.mass_low.min()) & 
                (cenque.mass <= mass_bins.mass_high.max()) & 
                (cenque.child >= 0) 
                ) 
    
        self.sample_trim(within_massbin_with_child)

    ngal_tot = len(within_massbin_with_child[0])

    if cenque.gal_type is None:         
        cenque.gal_type = np.array(['' for i in xrange(ngal_tot)], dtype='|S16') 
        cenque.sfr = np.array([-999. for i in xrange(ngal_tot)]) 
        cenque.ssfr = np.array([-999. for i in xrange(ngal_tot)]) 

    for attrib in ['gal_type', 'sfr', 'ssfr']: 
        if attrib not in cenque.data_columns: 
            cenque.data_columns.append(attrib)
    
    # simplest SFR assignment for starforming galaxies. Use mu_SFR(M*, z)
    # and randomly sampled normal delta_SFR. 
    if kw_sfr == 'average': 

        self.delta_sfr = np.array([-999. for i in xrange(ngal_tot)])

        for attrib in ['delta_sfr']: 
            if attrib not in cenque.data_columns: 
                cenque.data_columns.append(attrib)
    else: 
        raise NotImplementedError() 

    for i_m in xrange(mass_bins.nbins):              

            massbin_unassigned = np.where(
                    (cenque.mass > mass_bin_low[i_m]) &
                    (cenque.mass <= mass_bin_high[i_m]) &
                    (cenque.gal_type == '')
                    )

            if ngal_massbin == 0: 
                continue 
    
            # quiescent fraction for mass bin at z_snapshot 
            qf_massbin = util.get_fq(
                    mass_bin_mid[i_m], 
                    cenque.zsnap, 
                    lit = 'wetzelsmooth'
                    ) 
            """
            mass_bin_qf = util.get_fquenching(mass_bins.mass_mid[i_m], self.zsnap, 
                    yint=kwargs['fqing_yint'])
            """

            # Ngal(M_mid), Ngal,Q(M_mid) and Ngal,SF(M_mid)
            ngal_massbin = len(massbin_unassigned[0])
            ngal_q_massbin= int( 
                    np.rint(qf_massbin * np.float(ngal_massbin)) 
                    ) 
            ngal_sf_massbin = ngal_massbin - ngal_q_massbin

            if not quiet: 
                print mass_bin_low[i_m], ' < M < ', mass_bin_high[i_m]
                print 'fQ = ', qf_massbin, ' Ngal = ', ngal_massbin
                print 'Ngal,Q = ', ngal_q_massbin, ' Ngal,SF = ', ngal_sf_massbin 
            
            # Randomly select ngal_q_massbin quiescent galaxies from the 
            # massbin. Assign them 'quiescent' gal_type and sSFR and SFR 
            # based on a predetermined gaussian distribution about sSFR.
            if ngal_q_massbin > 0: 

                try: 
                    q_massbin = random.sample(xrange(ngal_massbin), ngal_q_massbin) 
                except ValueError: 
                    q_massbin = xrange(ngal_massbin) 

                i_q_massbin = massbin_unassigned[q_massbin]

                cenque.gal_type[i_q_massbin] = 'quiescent'   

                # sample SSFR from log-normal distribution with a preset
                # 0.18 dex scatter centered about some predetermined mean
                mu_q_ssfr = util.get_q_ssfr_mean(
                        cenque.mass[i_q_massbin]
                        ) 
                cenque.ssfr[i_q_massbin] = 0.18 * np.random.randn(mass_bin_n_q) + q_ssfr_mean
                cenque.sfr[i_q_massbin]  = cenque.ssfr[i_q_massbin] + cenque.mass[i_q_massbin]
            
            # ngal_sf_massbin starforming galaxies from the massbin. Assign 
            # them 'star-forming' gal_type and sSFR and SFR in some manner
            if ngal_sf_massbin > 0: 

                try: 
                    sf_massbin = [x for x in xrange(ngal_massbin) if x not in q_massbin]
                except NameError:       
                    sf_massbin = xrange(ngal_massbin)

                i_sf_massbin = massbin_unassigned[sf_massbin]
                
                cenque.gal_type[i_sf_massbin] = 'star-forming'
                
                # sample SFR 
                if  kw_sfr == 'average': 

                    mu_sf_sfr, sigma_sf_sfr = util.get_sfr_mstar_z_bestfit(
                            mass_bin_mid[i_m], 
                            cenque.zsnap, 
                            Mrcut=18
                            ) 


                    cenque.delta_sfr[i_sf_massbin] = sigma_sf_sfr * np.random.randn(ngal_sf_massbin)

                    cenque.sfr[i_sf_massbin] = mu_sf_sfr + cenque.delta_sfr[i_sf_massbin]

                    if not quiet:
                        print 'Starforming galaxies: '
                        print 'Average(SFR) = ', mu_sf_sfr, ' sigma(SFR) = ', sigma_sf_sfr 
                
                else:
                    raise NotImplementedError()
                    """
                    # Fluctuating SFR with random amplitude, frequency, and phase 

                    # sfr amplitude, frequency and phase respectively 
                    self.sfr_amp[mass_bin_sf_index] = \
                            sf_sig_sfr * np.random.randn(mass_bin_n_sf) 
                    self.sfr_freq[mass_bin_sf_index] = \
                            (2.0 * np.pi)/np.random.uniform(0.01, 0.1, mass_bin_n_sf)  
                    self.sfr_phase[mass_bin_sf_index] = \
                            np.random.uniform(0.0, 1.0, mass_bin_n_sf)

                    self.sfr[mass_bin_sf_index] = util.sfr_squarewave(
                            self.mass[mass_bin_sf_index], self.t_cosmic, 
                            amp = self.sfr_amp[mass_bin_sf_index], 
                            freq = self.sfr_freq[mass_bin_sf_index], 
                            phase = self.sfr_phase[mass_bin_sf_index])
                    """

                cenque.ssfr[i_sf_massbin] = cenque.sfr[i_sf_massbin] - cenque.mass[i_sf_massbin]
    
    # double check that SF assign didn't fail anywhere
    assign_fail = np.where(
            (cenque.sfr == -999.0) | 
            (cenqe.ssfr == -999.0)
            )
    if len(assign_fail[0]) > 0: 
        raise NameError('Function failed!')

    print 'Assign SFR function takes', (start_time - time.time())/60.0, ' minutes'

    return cenque
