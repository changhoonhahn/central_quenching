'''

Star-Formation main sequence codes 

Author(s): ChangHoon Hahn

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import os 
import h5py
from scipy import signal

# --- Local ---
from util import mpfit
from util.mass_bins import simple_mass_bin
from group_catalog import group_catalog as cq_group
from group_catalog.group_catalog import sf_centrals 
from util import cenque_utility as util 
from defutility.fitstables import mrdfits

def get_sfr_mstar_z_groupcat(m_stars, Mrcut=18, clobber=False): 
    ''' SFR(M*, z) from SDSS Group Catalog SFMS. calculate 
    the average SFR, sigma SFR, and Ngal in the mass bin encapsulating 
    given M*. 

    ----------------------------------------------------------------
    Parameters 
    ----------------------------------------------------------------
    m_star: stellar mass 
    Mrcut : absolute magnitude cut that specifies the group catalog
    clobber : If True,____. If False 

    ----------------------------------------------------------------
    Return
    ----------------------------------------------------------------
    sfr_out: star formation rate

    ----------------------------------------------------------------
    Notes
    ----------------------------------------------------------------
    * Average is skewed due to the Star-Forming classification which includes many 'transitioning galaxies' 

    '''

    massbins = simple_mass_bin()

    sf_data = sf_centrals(Mrcut=Mrcut, clobber=clobber) 
    
    if isinstance(m_stars, float): 
        m_stars = [m_stars]
    
    avg_sfrs, sig_sfrs, n_gals = [], [], [] 

    for m_star in m_stars:

        mbin_index = np.where(
                (massbins.mass_low <= m_star) &
                (massbins.mass_high > m_star) 
                )

        if (len(mbin_index[0]) != 1):
            avg_sfrs.append(-999.)
            sig_sfrs.append(-999.)

            continue

        mass_low  = (massbins.mass_low[mbin_index])[0]
        mass_high = (massbins.mass_high[mbin_index])[0]

        in_massbin = np.where(
                (sf_data.mass >= mass_low) & 
                (sf_data.mass < mass_high)
                )

        sfr_massbin = sf_data.sfr[in_massbin]

        n_gal = len(in_massbin[0])

        if n_gal > 30: 
            
            avg_sfrs.append( 
                    np.median(sfr_massbin)
                    )
            sig_sfrs.append(
                    np.std(sfr_massbin)
                    )
            
            print np.average(sfr_massbin), np.median(sfr_massbin)
            # Remove outliers by recomputing summary statistics 
            # for +/- 3 sigma
            #for iter in xrange(3): 
            #    sfr_range = np.where(
            #            (sfr_massbin > med_sfr - 3.0 * var_sfr) & 
            #            (sfr_massbin < med_sfr + 3.0 * var_sfr)
            #            ) 
            #    
            #    if np.median(sfr_massbin[sfr_range]) > med_sfr: 
            #        med_sfr = np.median(sfr_massbin[sfr_range]) 
            #        avg_sfr = np.average(sfr_massbin[sfr_range]) 
            #        sig_sfr = np.std(sfr_massbin[sfr_range])

        else: 
            avg_sfrs.append(-999.)
            sig_sfrs.append(-999.)

    return [avg_sfrs, sig_sfrs, n_gals] 

def get_sfr_mstar_z_envcount(m_stars, z_ins): 
    ''' SFR(m_star, z_in) from the SDSS/PRIMUS envcount catalog. 
    In this calcluation, we use an isolation criteria of central 
    galaxies for envcount catalog (envcount = 0). 

    ----------------------------------------------------------------
    Parameters 
    ----------------------------------------------------------------
    m_star: stellar mass 
    z_in: redshift
    
    ----------------------------------------------------------------
    Return
    ----------------------------------------------------------------
    sfr_out: star formation rate
    '''
    if isinstance(m_stars, float): 
        m_stars = [m_stars] 

    if len(m_stars) != len(z_ins): 
        raise ValueError()

    file_dir = 'dat/wetzel_tree/envcount/'

    sdss_file = ''.join([file_dir, 
            'envcount_cylr2.5h35_thresh75_sdss_active_z0.05_0.12_primuszerr.fits']) 
    primus_file = ''.join([file_dir, 
            'envcount_cylr2.5h35_thresh75_active_z0.2_1.0_lit.fits']) 

    if np.max(z_ins) < 0.2:
        sdss_sf_data = mrdfits(sdss_file) 
    elif np.min(z_ins) > 0.2: 
        primus_sf_data = mrdfits(primus_file)
    else: 
        sdss_sf_data = mrdfits(sdss_file) 
        primus_sf_data = mrdfits(primus_file)
    
    massbins = simple_mass_bin()
    zbins_low = np.array([ 0.0, 0.2, 0.4, 0.6, 0.8 ])
    zbins_high = np.array([ 0.2, 0.4, 0.6, 0.8, 1.0 ])
    
    avg_sfrs, sig_sfrs, n_gals = [], [], [] 

    for i_obj in xrange(len(m_stars)): 

        mbin_index = np.where(
                (massbins.mass_low <= m_stars[i_obj]) &
                (massbins.mass_high > m_stars[i_obj]) 
                )
        
        zbin_index = np.where(
                (zbins_low <= z_ins[i_obj]) & 
                (zbins_high > z_ins[i_obj]) 
                )
        
        if (len(zbin_index[0]) != 1) or (len(mbin_index[0]) != 1):

            print m_stars[i_obj], z_ins[i_obj]

            avg_sfrs.append(-999.)
            sig_sfrs.append(-999.)

            continue

        mass_low  = (massbins.mass_low[mbin_index])[0]
        mass_high = (massbins.mass_high[mbin_index])[0]
        z_low = (zbins_low[zbin_index])[0]
        z_high = (zbins_high[zbin_index])[0]

        if z_ins[i_obj] < 0.2: 
            sf_data = sdss_sf_data
        else: 
            sf_data = primus_sf_data
        
        # slice data into bins of mass and redshift 
        mass_z_slice = np.where(
                (sf_data.mass >= mass_low) & (sf_data.mass < mass_high) & 
                (sf_data.redshift >= z_low) & (sf_data.redshift < z_high) & 
                (sf_data.envcount == 0.0) & (sf_data.mass > sf_data.masslimit) & 
                (sf_data.edgecut == 1) 
                )

        n_gal = len(mass_z_slice[0])

        n_gals.append(n_gal)

        if n_gal > 30: 

            avg_sfr = np.average(
                    sf_data.sfr[mass_z_slice], 
                    weights = sf_data.weight[mass_z_slice]
                    ) 
            avg_sfrs.append(
                    np.average(
                        sf_data.sfr[mass_z_slice], 
                        weights = sf_data.weight[mass_z_slice]
                        ) 
                    )

            sig_sfrs.append(
                    np.sqrt(np.average(
                        (sf_data.sfr[mass_z_slice] - avg_sfr)**2, 
                        weights = sf_data.weight[mass_z_slice]
                        ))
                    )
        else: 
            avg_sfrs.append(-999.)
            sig_sfrs.append(-999.)

    return [avg_sfrs, sig_sfrs, n_gals] 


"""
def sf_duty_test(): 
    ''' Test of Star-Formation duty cycle

    '''
    prettyplot()
    pretty_colors = prettycolors()

    n_snaps, z_snap, t_snap, t_wid = np.loadtxt('snapshot_table.dat', 
            unpack=True, usecols=[0, 2, 3, 4])

    t_snaps = t_snap[n_snaps < 13]
    nsnaps = n_snaps[n_snaps < 13]

    #sfr_gauss = 0.3 * np.random.randn(10000)   # SFR gaussian with 0.3 dex sigma 
    sfr_w = (2.0 * np.pi)/np.random.uniform(0.01, 0.1, 10000)  # frequency
    sfr_d = np.random.uniform(0.0, 1.0, 10000)
    sfr_A = 0.3 * np.random.randn(10000)
    sfr_t = lambda t: sfr_A * signal.square(sfr_w * (t - t_snap[n_snaps == 12] + sfr_d))

    #sfr_amp = sfr_gauss / np.sin( sfr_d ) 
    #sfr_t = lambda t: sfr_amp * np.sin( sfr_w * (t - t_snap[n_snaps == 12]) - sfr_d )

    #sfr_t = lambda t: sfr_A * signal.sawtooth(sfr_w * (t - t_snap[n_snaps == 12] + sfr_d), width=0.5)
    #sfr_t = lambda t: sfr_a * (t - t_snap[n_snaps == 12] + sfr_b) % sfr_A
    #sfr_t = lambda t: sfr_gauss + 0.3
    
    fig = plt.figure(1, figsize = (10,10))
    sub = fig.add_subplot(111)
    sfr_tsnap = np.zeros(10000) 
    for i_t, t in enumerate(t_snaps): 

        sfr_hist, sfr_bin_edges = np.histogram(sfr_tsnap + sfr_t(t), range=[-1.0, 1.0], bins=100)
        sfr_bin_low = sfr_bin_edges[:-1]
        sfr_bin_high = sfr_bin_edges[1:]
        sfr_bin_mid = [ 0.5*(sfr_bin_low[i] + sfr_bin_high[i]) 
                for i in range(len(sfr_bin_low)) ] 
    
        sub.plot(sfr_bin_mid, sfr_hist, color=pretty_colors[i_t+1], lw=2, label='Snapshot '+str(n_snaps[i_t]))
    
        sfr_tsnap += sfr_t(t)
        print np.float(len(sfr_tsnap[(sfr_tsnap < 0.1) & (sfr_tsnap > -0.1)]))/np.float(len(sfr_tsnap))
    
    sfr_hist, sfr_bin_edges = np.histogram(0.3 * np.random.randn(10000), 
            range=[-1.0, 1.0], bins=100)
    sub.plot(sfr_bin_mid, sfr_hist, 
            color='black', lw=4, ls='--', label='Gaussian')
    sub.set_xlim([-1.0, 1.0])
    sub.set_xlabel('SFR - average SFR') 
    sub.legend(loc='upper right') 

    fig2 = plt.figure(2, figsize=(10,10))
    sub = fig2.add_subplot(111)

    for i in range(1,11): 
        sub.plot(np.arange(t_snap[n_snaps==12], 15., 0.01), [ sfr_t(t)[i] for t in np.arange(t_snap[n_snaps==12], 15., 0.01) ], 
                color=pretty_colors[i+1]) 
    sub.set_xlim([6.8, 9.0])
    sub.set_ylabel('SFR - average SFR') 
    sub.set_xlabel('cosmic time (Gyr)') 
    plt.show() 

def get_ssfr_mstar_qgroupcat(m_star, Mrcut=18, clobber=False): 
    ''' Return SSFR(m_star) from the Quiescent SDSS Group Catalog

    Parameters 
    ----------
    m_star: stellar mass 
    Mrcut : absolute magnitude cut that specifies the group catalog
    clobber : If True,____. If False 
    
    Return
    ------
    ssfr_out: Star formation rate (SFR) at m_star

    Notes
    -----
    * Output sSFR is the median sSFR of the sSFR in the mass bin  

    '''
    q_data = q_centrals(Mrcut=Mrcut, clobber=clobber) 
    
    # determine mass bin and redshift bin to splice the data
    massbins = simple_mass_bin()
    mbin_index = (np.array(massbins.mass_low) <= m_star) & \
            (np.array(massbins.mass_high) > m_star) 
    mass_low = (np.array(massbins.mass_low)[mbin_index])[0]
    mass_high = (np.array(massbins.mass_high)[mbin_index])[0]

    # splice data
    bin_index = (q_data.mass >= mass_low) & \
            (q_data.mass < mass_high)

    n_gal = len(q_data.mass[bin_index])

    if len(q_data.sfr[bin_index]) > 0: 
        avg_ssfr = np.average(q_data.ssfr[bin_index]) 
        median_ssfr = np.median(q_data.ssfr[bin_index]) 
        var_ssfr = np.average( (q_data.ssfr[bin_index] - avg_ssfr)**2 ) 
    else: 
        avg_ssfr = -999.
        median_ssfr = -999.
        var_ssfr = -999.

    return [median_ssfr, var_ssfr, n_gal] 

"""
