"""

Assign StarFomation properties to CenQue object

Author(s): ChangHoon Hahn

"""

import sf_mainseq as sfms

def assign_sfr(cenque, **kwargs):
    """ Assign star-formation properties to CenQue object
    
    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    cenque : CenQue class object 

    """

    if cenque.mass == None: 
        raise ValueError()
    
    mass_bins = cenque.mass_bins

    # if slope/yint of the SF-MS is specified in kwargs, use those values
    # otherwise, just use the fit from the group catalogs (hardcoded)
    # groupcat_slope, groupcat_yint = sfms.get_bestfit_groupcat_sfms(Mrcut=18)
    # sfms_sfr_fit = sfms.get_sfmsfit_sfr(groupcat_slope, groupcat_yint)
    if 'sfms_slope' in kwargs.keys(): 
        sfms_sfr_fit = sfms.get_sfmsfit_sfr(
                kwargs['sfms_slope'], 
                kwargs['sfms_yint']
                )
        
    # remove galaxies below the minimum and maximum mass and galaxies without children 
    if (min(cenque.mass) < mass_bins.mass_low.min()) or (max(cenque.mass) > mass_bins.mass_high.max()): 
        
        within_massbin_with_child = np.where(
                (cenque.mass > mass_bins.mass_low.min()) & 
                (cenque.mass <= mass_bins.mass_high.max()) & 
                (cenque.child >= 0) 
                ) 
    
        self.sample_trim(within_massbin_with_child)

    if cenque.gal_type is None:         
        cenque.gal_type = np.array(['' for i in range(len(self.mass))], dtype='|S16') 
        cenque.sfr = np.array([-999. for i in range(len(self.mass))]) 
        cenque.ssfr = np.array([-999. for i in range(len(self.mass))]) 

    # SFR parameters
    if kwargs['sfr'] == 'sfr_func': 
        # periodic SFR 
        self.sfr_amp = np.array([-999. for i in range(len(self.mass))])
        self.sfr_freq = np.array([-999. for i in range(len(self.mass))])
        self.sfr_phase = np.array([-999. for i in range(len(self.mass))])

    elif kwargs['sfr'] == 'sfr_avg': 
        self.sfr_resid = np.array([-999. for i in range(len(self.mass))])
            
    else: 
        raise NotImplementedError('asdlkfj') 

        
    
    
    
    
    
    
    
    # loop through mass bins and assign SFRs ---------------------------------------------
        for i_m in range(mass_bins.nbins):              
        
            mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
            mass_bin_high = round(mass_bins.mass_high[i_m], 2) 
            # boolean list for mass range
            mass_bin_bool = (self.mass > mass_bin_low) & \
                    (self.mass <= mass_bin_high) &\
                    (self.gal_type == '')           # only assign for galaxies without galtype

            # indices of galaxies within mass range
            mass_bin_index = np.array(range(len(self.mass)))[mass_bin_bool] 
        
            mass_bin_ngal = np.sum(mass_bin_bool)   # Ngal in mass bin 
            if mass_bin_ngal == 0:                  # if there are no galaxies within mass bin
                continue 
    
            # get quiescent fraction for mass bin at z_snapshot 
            #mass_bin_qf = util.get_fquenching(mass_bins.mass_mid[i_m], self.zsnap, 
            #        yint=kwargs['fqing_yint'])
            mass_bin_qf = util.get_fq(mass_bins.mass_mid[i_m], self.zsnap, lit='wetzelsmooth') 

            # Ngal,quiescent in mass bin
            mass_bin_n_q = int( np.rint(mass_bin_qf * np.float(mass_bin_ngal)) ) 

            # Ngal,active in mass bin 
            mass_bin_n_sf = mass_bin_ngal - mass_bin_n_q

            #print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]
            #print 'fQ = ', mass_bin_qf, ' Ngal = ', mass_bin_ngal, \
            #        ' Ngal,Q = ', mass_bin_n_q, ' Ngal,SF = ', mass_bin_n_sf
            
            #### Quiescent ####
            if mass_bin_n_q > 0: 
                # randomly sample mass_bin_n_q galaxies from the mass bin 
                try: 
                    mass_bin_q_index = random.sample(mass_bin_index, mass_bin_n_q) 
                except ValueError: 
                    mass_bin_q_index = mass_bin_index 

                self.gal_type[mass_bin_q_index] = 'quiescent'   # label galaxy type 

                '''
                sample SSFR from log-normal distribution with 0.2dex scatter 
                centered about some determined mean
                '''
                q_ssfr_mean = util.get_q_ssfr_mean(self.mass[mass_bin_q_index]) 
                self.ssfr[mass_bin_q_index] = 0.18 * np.random.randn(mass_bin_n_q) + q_ssfr_mean
                # calculate SFR by SSFR+Mass
                self.sfr[mass_bin_q_index] = self.ssfr[mass_bin_q_index] + self.mass[mass_bin_q_index]

            #### Star-Forming ####
            if mass_bin_n_sf > 0: 
                try: 
                    mass_bin_sf_index = [x for x in mass_bin_index if x not in mass_bin_q_index]
                except NameError:       # if there aren't any quiescent galaxies
                    mass_bin_sf_index = mass_bin_index
                
                self.gal_type[mass_bin_sf_index] = 'star-forming'   # label galaxy type 
            
                if kwargs['sfr'] == 'sfr_avg': 
                    #get average and scatter of SF main sequence 
                    [sf_avg_sfr, sf_sig_sfr] = util.get_sfr_mstar_z_bestfit(
                            mass_bins.mass_mid[i_m], self.zsnap, Mrcut=18) 

                    #print 'SF Average(SFR) = ', sf_avg_sfr, ' sigma_SFR = ', sf_sig_sfr

                    self.sfr_resid[mass_bin_sf_index] = sf_sig_sfr * np.random.randn(mass_bin_n_sf)
                    self.sfr[mass_bin_sf_index] = util.sfr_avg_residual(
                            mass_bins.mass_mid[i_m], self.zsnap, 
                            resid = self.sfr_resid[mass_bin_sf_index]) 
                    #self.sfr[mass_bin_sf_index] = sf_sig_sfr * np.random.randn(mass_bin_n_sf) + sf_avg_sfr 

                elif kwargs['sfr'] == 'sfr_func': 
                    # Fluctuating SFR with random amplitude, frequency, and phase 

                    # average and 1sigma of SF-MS 
                    [sf_avg_sfr, sf_sig_sfr] = util.get_sfr_mstar_z_bestfit(
                            mass_bins.mass_mid[i_m], self.zsnap, Mrcut=18) 

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
                else: 
                    raise NotImplementedError('asdfkjlkjasdf') 

                self.ssfr[mass_bin_sf_index] = \
                        self.sfr[mass_bin_sf_index] - self.mass[mass_bin_sf_index]

        # check for SFR/SSFR assign fails
        assign_fail = (self.sfr == -999.0) | (self.ssfr == -999.0)
        if np.sum(assign_fail) > 0: 
            raise NameError('asdfasdfasdfasdf')

