"""

Assign StarFomation properties to CenQue object

Author(s): ChangHoon Hahn

"""
import time
import random 
import numpy as np
import warnings

# --- Local ----
from quiescent_fraction import get_fq
from util import cenque_utility as util
from util.tau_quenching import get_quenching_efold
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
            warnings.warn("SF properties do not match")

    if 'fq_prop' not in cenque.__dict__.keys():
        cenque.fq_prop = fq_prop 
    else: 
        if cenque.fq_prop != fq_prop:
            warnings.warn("fQ properties do not match")
    
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

    for attrib in ['sf_prop', 'fq_prop']:
        if attrib not in cenque.metadata: 
            cenque.metadata.append(attrib)

    for attrib in ['gal_type', 'sfr', 'ssfr']: 
        if attrib not in cenque.data_columns: 
            cenque.data_columns.append(attrib)

    if cenque.gal_type is None:         
        cenque.gal_type = np.array(['' for i in xrange(ngal_tot)], dtype='|S16') 
        cenque.sfr      = np.array([-999. for i in xrange(ngal_tot)]) 
        cenque.ssfr     = np.array([-999. for i in xrange(ngal_tot)]) 
    
    # simplest SFR assignment for starforming galaxies. Use mu_SFR(M*, z)
    # and randomly sampled normal delta_SFR. 
    if sf_prop['name'] == 'average': 
        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

        if 'delta_sfr' not in cenque.__dict__.keys(): 
            cenque.delta_sfr = np.array([-999. for i in xrange(ngal_tot)])
        
        if 'avg_sfr' not in cenque.__dict__.keys(): 
            cenque.avg_sfr = np.array([-999. for i in xrange(ngal_tot)])

        extra_attr = ['avg_sfr', 'delta_sfr']
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
                )[0]
            for i_m in xrange(mass_bins.nbins)
            ]

    # Ngal(M_mid), Ngal,Q(M_mid) and Ngal,SF(M_mid)
    ngal_massbin = np.array([x.size for x in massbin_unassigned])
    ngal_q_massbin = np.rint(qf_massbin * ngal_massbin.astype(float)).astype(int)
    ngal_sf_massbin = ngal_massbin - ngal_q_massbin

    # fail-safe for ngal_q_massbin
    if len(np.where(ngal_q_massbin > ngal_massbin)[0]) > 0: 
        ngal_q_massbin[np.where(ngal_q_massbin > ngal_massbin)] = ngal_massbin[np.where(ngal_q_massbin > ngal_massbin)]

    for i_m in xrange(mass_bins.nbins):             

        if len(massbin_unassigned[i_m]) == 0: 
            continue

        begin_loop_time = time.time()
        if not quiet: 
            print mass_bin_low[i_m], ' < M < ', mass_bin_high[i_m]
            print 'fQ = ', qf_massbin[i_m], ' Ngal = ', ngal_massbin[i_m]
            print 'Ngal,Q = ', ngal_q_massbin[i_m], ' Ngal,SF = ', ngal_sf_massbin[i_m]

        shuffled_massbin_index = np.arange(ngal_massbin[i_m])
        np.random.seed()
        np.random.shuffle(shuffled_massbin_index)
        i_q_end = ngal_q_massbin[i_m]
        
        # Randomly select ngal_q_massbin quiescent galaxies from the 
        # massbin. Assign them 'quiescent' gal_type and sSFR and SFR 
        # based on a predetermined gaussian distribution about sSFR.
        if ngal_q_massbin[i_m] > 0: 

            q_massbin = shuffled_massbin_index[:i_q_end]
            i_q_massbin = (massbin_unassigned[i_m])[q_massbin]

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
        
        # ngal_sf_massbin starforming galaxies from the massbin. Assign 
        # them 'star-forming' gal_type and sSFR and SFR in some manner
        if ngal_sf_massbin[i_m] > 0: 

            #sf_massbin = shuffled_massbin_index[i_q_end:i_sf_end]
            sf_massbin = shuffled_massbin_index[i_q_end:]
            i_sf_massbin = (massbin_unassigned[i_m])[sf_massbin]
            
            cenque.gal_type[i_sf_massbin] = 'star-forming'
            
            # sample SFR 
            if  sf_prop['name'] == 'average': 

                mu_sf_sfr = sfr_mstar_z(mass_bin_mid[i_m], cenque.zsnap)
                sigma_sf_sfr = sig_sfr_mstar_z(mass_bin_mid[i_m], cenque.zsnap)

                mu_sf_ssfr = mu_sf_sfr - mass_bin_mid[i_m]
                
                cenque.avg_sfr[i_sf_massbin] = mu_sf_sfr
                cenque.delta_sfr[i_sf_massbin] = sigma_sf_sfr * np.random.randn(ngal_sf_massbin[i_m])

                cenque.sfr[i_sf_massbin] = mu_sf_sfr + cenque.delta_sfr[i_sf_massbin]

                if not quiet:
                    print 'Starforming galaxies: '
                    print 'Average(SFR) = ', mu_sf_sfr, ' sigma(SFR) = ', sigma_sf_sfr 
            
            else:
                raise NotImplementedError

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

"""
Attempt at adding the green valley to the assign SFR code. 
----------------------------------------------------------
    if 'tau_prop' not in cenque.__dict__.keys():
        cenque.tau_prop = tau_prop 
    else: 
        if cenque.tau_prop != tau_prop:
            warnings.warn("tau properties do not match")

        tau_prop = cenque.tau_prop
    
    if 'prev_fquenching' not in cenque.__dict__.keys():
        cenque.prev_fquenching = 0.01*(mass_bin_mid - 9.2)
        #0.03*(mass_bin_mid - 9.4)
        #0.001+0.0*mass_bin_mid #* (mass_bin_mid - 9.4)

        f_g = (0.1/2.3) * (mass_bin_mid - 9.2) #1.5 * 0.1 + 0.0 * mass_bin_mid  #0.01*(mass_bin_mid - 9.0)#cenque.prev_fquenching
        f_g[np.isnan(f_g)] = 0.0
    else: 
        f_g = 0.0
    
    
    ngal_g_massbin = np.rint(f_g * ngal_massbin.astype(float)).astype(int)
    ngal_q_massbin = np.rint(qf_massbin * (ngal_massbin - ngal_g_massbin).astype(float)).astype(int)
    ngal_sf_massbin = ngal_massbin - ngal_g_massbin - ngal_q_massbin


        #i_sf_end = ngal_q_massbin[i_m] + ngal_sf_massbin[i_m]

        # artificial green valley 
        if ngal_g_massbin[i_m] > 0: 
            g_massbin = shuffled_massbin_index[i_sf_end:]
            i_g_massbin = (massbin_unassigned[i_m])[g_massbin]
            
            cenque.gal_type[i_g_massbin] = 'quiescent'

            mu_q_ssfr = util.get_q_ssfr_mean(cenque.mass[i_g_massbin]) 
            #mu_g_ssfr = 0.5 * (mu_q_ssfr + mu_sf_ssfr)
            #sigma_g_ssfr = (mu_sf_ssfr - mu_q_ssfr)/3.
        
            #cenque.ssfr[i_g_massbin] = sigma_g_ssfr * np.random.randn(ngal_g_massbin[i_m]) + mu_g_ssfr 
            #cenque.sfr[i_g_massbin]  = cenque.ssfr[i_g_massbin] + cenque.mass[i_g_massbin]
            
            cenque.ssfr[i_g_massbin] = (mu_sf_ssfr - mu_q_ssfr) * np.random.uniform(size=len(mu_q_ssfr)) + mu_q_ssfr
            cenque.sfr[i_g_massbin]  = cenque.ssfr[i_g_massbin] + cenque.mass[i_g_massbin]

            cenque.tau[i_g_massbin] = get_quenching_efold(cenque.mass[i_g_massbin], tau_param = tau_prop)
            cenque.q_ssfr[i_g_massbin] = 0.18 * np.random.randn(ngal_g_massbin[i_m]) + mu_q_ssfr 
""" 
