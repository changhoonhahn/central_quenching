'''

Star-Formation main sequence codes 

Author(s): ChangHoon Hahn

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import os 
import mpfit
import cenque_utility as util 


def get_sfr_mstar_z(m_star, z_in, machine='harmattan'): 
    ''' Get SFR(m_star, z_in) from SDSS/PRIMUS envcount data 

    input
    -----
    m_star: stellar mass 
    z_in: redshift
    
    output
    ------
    sfr_out: star formation rate
    '''
    # PRIMUS Star-Forming environment count data  
    if machine == 'harmattan': 
        file_dir = '/data1/hahn/wetzel_tree/envcount/'
    else: 
        raise NotImplementedError('asdfasdf')
    if z_in < 0.2:
        file = ''.join([file_dir, 
            'envcount_cylr2.5h35_thresh75_sdss_active_z0.05_0.12_primuszerr.fits']) 
    else: 
        file = ''.join([file_dir, 
            'envcount_cylr2.5h35_thresh75_active_z0.2_1.0_lit.fits']) 

    sf_data = util.mrdfits(file) 
    
    # determine mass bin and redshift bin to splice the data
    massbins = util.simple_mass_bin()
    mbin_index = (np.array(massbins.mass_low) <= m_star) & \
            (np.array(massbins.mass_high) > m_star) 

    mass_low = (np.array(massbins.mass_low)[mbin_index])[0]
    mass_high = (np.array(massbins.mass_high)[mbin_index])[0]

    zbins_low = np.array([ 0.0, 0.2, 0.4, 0.6, 0.8 ])
    zbins_high = np.array([ 0.2, 0.4, 0.6, 0.8, 1.0 ])
    zbin_index = (zbins_low <= z_in) & (zbins_high > z_in) 

    z_low = (zbins_low[zbin_index])[0]
    z_high = (zbins_high[zbin_index])[0]
    
    # splice data
    bin_index = (sf_data.mass >= mass_low) & \
            (sf_data.mass < mass_high) & \
            (sf_data.redshift >= z_low) & \
            (sf_data.redshift < z_high) & \
            (sf_data.envcount == 0.0) & \
            (sf_data.mass > sf_data.masslimit) & \
            (sf_data.edgecut == 1) 

    n_gal = len(sf_data.weight[bin_index])

    if len(sf_data.sfr[bin_index]) > 0: 
        avg_sfr = np.average(sf_data.sfr[bin_index], weights=sf_data.weight[bin_index]) 
        var_sfr = np.average( (sf_data.sfr[bin_index] - avg_sfr)**2, 
                weights=sf_data.weight[bin_index]) 
    else: 
        avg_sfr = -999.
        var_sfr = -999.

    return [avg_sfr, var_sfr, n_gal] 

def sdss_groupcat_sfms_bestfit():
    ''' SF MS linear fit 
    Using sSFR from SDSS Group Catalog 
    '''
    mass = np.array([10.25, 10.75, 11.25])
    SFR = np.array([-0.0777, 0.248, 0.483])

    p0 = [0.5607, 0.0775917]
    fa = {'x': mass-10.5, 'y': SFR}
    bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, nprint=0) 

    return [bestfit.params[0], bestfit.params[1]]

def line_fixedslope(x, p): 
    ''' Line with slope of SF MS linear fit 
    '''
    slope, yint = sdss_groupcat_sfms_bestfit()
    return slope * x + p[0]

def mpfit_line_fixedslope(p, fjac=None, x=None, y=None, err=None): 
    model = line_fixedslope(x, p) 
    status = 0 
    return([status, (y-model)/err]) 

def build_sfr_mstar_z(lowz_slope, lowz_yint):
    ''' Given z~0.1 SF MS slope and yint
    fit fixed slope line to other redshift limits
    
    input 
    -----
    slope: slope of z ~ 0.1 SF-MS  
    yint: yint of z ~ 0.1 SF-MS

    '''

    mass_bin = util.simple_mass_bin() 
    zbins = [ (0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0) ] 
    
    zmids = [] 
    slopes = [] 
    yints = []
    for zbin in zbins: 
        
        mass = [] 
        avg_sfrs = []
        var_sfrs = [] 

        for i_mass in range(len(mass_bin.mass_low)): 

            avg_sfr, var_sfr, ngal = get_sfr_mstar_z(mass_bin.mass_mid[i_mass], 
                    0.5*(zbin[0] + zbin[1]))

            if ngal < 10: 
                continue 
            if mass_bin.mass_mid[i_mass] < 9.5: 
                continue

            mass.append(mass_bin.mass_mid[i_mass]) 
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)

        p0 = [0.1] 
        fa = {'x': np.array(mass)-10.5, 'y': np.array(avg_sfrs), 'err': np.array(var_sfrs)} 

        bestfit = mpfit.mpfit(mpfit_line_fixedslope, p0, functkw=fa, nprint=0)
        
        zmids.append(0.5 * (zbin[0] + zbin[1])) 
        slopes.append(lowz_slope) 
        if zbin[0] == 0.0: 
            yints.append(lowz_yint) 
        else: 
            yints.append(bestfit.params[0]) 

    return [zmids, slopes, yints]

if __name__=='__main__':
    # SDSS group catalog best fit 
    groupcat_slope, groupcat_yint = sdss_groupcat_sfms_bestfit()
    print build_sfr_mstar_z(groupcat_slope, groupcat_yint)
