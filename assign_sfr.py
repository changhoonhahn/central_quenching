"""

Assign StarFomation properties to CenQue object

Author(s): ChangHoon Hahn

"""
import time
import random 
import numpy as np

# --- Local ----
from quiescent_fraction import get_fq
from util import cenque_utility as util
from sfms.fitting import get_param_sfr_mstar_z

def assign_sfr(
        cenque, 
        sf_prop={'name': 'average'}, 
        fq_prop={'name': 'wetzelsmooth'}, 
        quiet=True, 
        **kwargs
        ):
    """ Assign star-formation properties to CenQue object. 

    The function Goes through mass bins and then classifies unassigned 
    galaxies into quiescent/star-forming based on quiescent fraction 
    function. 
    
    Then SF properties are assigned to the galaxies. For quiescent 
    galaxies the sSFR is drawn from a log normal distribution about 
    some predeterimined mu_sSFR. For star-forming galaxies, SFR is
    sampled from a designated SF-MS model. 
    
    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    cenque : CenQue class object 
    ----------------------------------------------------------------
    Notes 
    ----------------------------------------------------------------
    * Re imagine how to set up SFR 

    """

    # time the code 
    start_time = time.time()

    if cenque.mass == None: 
        raise ValueError()

    if cenque.zsnap == None: 
        raise ValueError()
    
    if 'sf_prop' not in cenque.__dict__.keys(): 
        cenque.sf_prop = sf_prop
    else: 
        if cenque.sf_prop != sf_prop: 
            raise ValueError()

    if 'fq_prop' not in cenque.__dict__.keys():
        cenque.fq_prop = fq_prop 
    else: 
        if cenque.fq_prop != fq_prop:
            raise ValueError()
    
    mass_bins = cenque.mass_bins
    mass_bin_low  = mass_bins.mass_low
    mass_bin_mid  = mass_bins.mass_mid
    mass_bin_high = mass_bins.mass_high

    # remove galaxies below the minimum and maximum mass and galaxies without children 
    within_massbin_with_child = np.where(
            (cenque.mass > mass_bins.mass_low.min()) & 
            (cenque.mass <= mass_bins.mass_high.max()) & 
            (cenque.child >= 0) 
            ) 

    if (min(cenque.mass) < mass_bins.mass_low.min()) or (max(cenque.mass) > mass_bins.mass_high.max()): 
        cenque.sample_trim(within_massbin_with_child)

    ngal_tot = len(within_massbin_with_child[0])

    for attrib in ['gal_type', 'sfr', 'ssfr']: 
        if attrib not in cenque.data_columns: 
            cenque.data_columns.append(attrib)

    if cenque.gal_type is None:         
        cenque.gal_type = np.array(['' for i in xrange(ngal_tot)], dtype='|S16') 
        cenque.sfr = np.array([-999. for i in xrange(ngal_tot)]) 
        cenque.ssfr = np.array([-999. for i in xrange(ngal_tot)]) 
    
    # simplest SFR assignment for starforming galaxies. Use mu_SFR(M*, z)
    # and randomly sampled normal delta_SFR. 
    if sf_prop['name'] == 'average': 
        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

        if 'delta_sfr' not in cenque.__dict__.keys(): 
            cenque.delta_sfr = np.array([-999. for i in xrange(ngal_tot)])

        extra_attr = ['delta_sfr']
    else: 
        raise NotImplementedError() 

    for attrib in extra_attr: 
        if attrib not in cenque.data_columns: 
            cenque.data_columns.append(attrib)

    # f_Q(M_*,mid, z_snapshot) 
    qf_massbin = get_fq(
            mass_bin_mid, 
            cenque.zsnap, 
            lit = fq_prop['name']
            ) 

    massbin_unassigned = [
            np.where(
                (cenque.mass > mass_bin_low[i_m]) &
                (cenque.mass <= mass_bin_high[i_m]) &
                (cenque.gal_type == '')
                )
            for i_m in xrange(mass_bins.nbins)
            ]

    # Ngal(M_mid), Ngal,Q(M_mid) and Ngal,SF(M_mid)
    ngal_massbin = np.array(
            [len(x[0]) for x in massbin_unassigned]
            )
    ngal_q_massbin= np.array([
        int(np.rint(qf_massbin[i_m] * np.float(ngal_massbin[i_m])))
        for i_m in xrange(mass_bins.nbins)
        ]) 
    ngal_sf_massbin = ngal_massbin - ngal_q_massbin
    # fail-safe for ngal_q_massbin
    ngal_q_massbin[np.where(ngal_q_massbin > ngal_massbin)] = ngal_massbin[np.where(ngal_q_massbin > ngal_massbin)]

    for i_m in xrange(mass_bins.nbins):             

        if len(massbin_unassigned[i_m][0]) == 0: 
            continue

        begin_loop_time = time.time()
        if not quiet: 
            print mass_bin_low[i_m], ' < M < ', mass_bin_high[i_m]
            print 'fQ = ', qf_massbin[i_m], ' Ngal = ', ngal_massbin[i_m]
            print 'Ngal,Q = ', ngal_q_massbin[i_m], ' Ngal,SF = ', ngal_sf_massbin[i_m]

        shuffled_massbin_index = np.arange(ngal_massbin[i_m])
        np.random.shuffle(shuffled_massbin_index)
        i_q_end = ngal_q_massbin[i_m]
        
        # Randomly select ngal_q_massbin quiescent galaxies from the 
        # massbin. Assign them 'quiescent' gal_type and sSFR and SFR 
        # based on a predetermined gaussian distribution about sSFR.
        if ngal_q_massbin[i_m] > 0: 

            q_massbin = shuffled_massbin_index[:i_q_end]
            i_q_massbin = (massbin_unassigned[i_m][0])[q_massbin]

            cenque.gal_type[i_q_massbin] = 'quiescent'   

            # sample SSFR from log-normal distribution with a preset
            # 0.18 dex scatter centered about some predetermined mean
            mu_q_ssfr = util.get_q_ssfr_mean(
                    cenque.mass[i_q_massbin]
                    ) 
            if len(i_q_massbin) != ngal_q_massbin[i_m]:
                print shuffled_massbin_index
                print q_massbin
                print i_q_massbin
                print i_q_end, ngal_q_massbin[i_m]
            cenque.ssfr[i_q_massbin] = 0.18 * np.random.randn(ngal_q_massbin[i_m]) + mu_q_ssfr 
            cenque.sfr[i_q_massbin]  = cenque.ssfr[i_q_massbin] + cenque.mass[i_q_massbin]
            q_time = time.time()
        
        # ngal_sf_massbin starforming galaxies from the massbin. Assign 
        # them 'star-forming' gal_type and sSFR and SFR in some manner
        if ngal_sf_massbin[i_m] > 0: 

            sf_massbin = shuffled_massbin_index[i_q_end:]
            i_sf_massbin = (massbin_unassigned[i_m][0])[sf_massbin]
            
            cenque.gal_type[i_sf_massbin] = 'star-forming'
            
            # sample SFR 
            if  sf_prop['name'] == 'average': 

                mu_sf_sfr = sfr_mstar_z(mass_bin_mid[i_m], cenque.zsnap)
                sigma_sf_sfr = sig_sfr_mstar_z(mass_bin_mid[i_m], cenque.zsnap)
                
                cenque.delta_sfr[i_sf_massbin] = sigma_sf_sfr * np.random.randn(ngal_sf_massbin[i_m])

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
            (cenque.ssfr == -999.0)
            )
    if len(assign_fail[0]) > 0: 
        raise NameError('Function failed!')
    
    if not quiet: 
        print 'Assign SFR function takes', (time.time()-start_time)/60.0, ' minutes'
    if 'evol_from' in cenque.cenque_type: 
        pass
    else: 
        cenque.cenque_type = 'sf_assigned'

    return cenque
