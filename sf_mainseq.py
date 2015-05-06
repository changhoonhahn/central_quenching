'''

Star-Formation main sequence codes 

Author(s): ChangHoon Hahn

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import os 
import mpfit
import h5py

# --- Local ---
import cenque as cq
import cenque_groupcat as cq_group
import cenque_utility as util 


def get_sfr_mstar_z_envcount(m_star, z_in, machine='harmattan'): 
    ''' Return SFR(m_star, z_in) from SDSS/PRIMUS envcount data 

    Parameters 
    ----------
    m_star: stellar mass 
    z_in: redshift
    machine : set to harmattan 
    
    Return
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

def line_fixedslope(x, p, slope=0.56): 
    ''' Line with slope of SF MS linear fit 
    '''
    return slope * x + p[0]

def mpfit_line_fixedslope(p, slope=0.56, fjac=None, x=None, y=None, err=None): 
    model = line_fixedslope(x, p, slope=slope) 
    status = 0 
    return([status, (y-model)/err]) 

def build_sfmsfit_sfr(lowz_slope, lowz_yint):
    ''' Given z~0.1 SF MS slope and yint fit fixed slope line to higher redshift bins
    to get best fit y-int values 
    
    Parameters
    -----------
    slope: slope of z ~ 0.1 SF-MS  
    yint: yint of z ~ 0.1 SF-MS

    '''
    fid_mass = 10.5

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

            avg_sfr, var_sfr, ngal = get_sfr_mstar_z_envcount(mass_bin.mass_mid[i_mass], 
                    0.5*(zbin[0] + zbin[1]))

            if ngal < 10: 
                continue 
            if mass_bin.mass_mid[i_mass] < 9.5: 
                continue

            mass.append(mass_bin.mass_mid[i_mass]) 
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)

        p0 = [0.1] 
        fa = {'slope': lowz_slope, 
                'x': np.array(mass)-fid_mass, 'y': np.array(avg_sfrs), 'err': np.array(var_sfrs)}

        bestfit = mpfit.mpfit(mpfit_line_fixedslope, p0, functkw=fa, nprint=0)
        
        zmids.append(0.5 * (zbin[0] + zbin[1])) 
        slopes.append(lowz_slope) 
        if zbin[0] == 0.0: 
            yints.append(lowz_yint) 
        else: 
            yints.append(bestfit.params[0]) 
    
    output_file = get_sfmsfit_slope_yint_file(lowz_slope, lowz_yint) 
    f = h5py.File(output_file, 'w') 
    grp = f.create_group('slope_yint')
    
    grp.attrs['fid_mass'] = fid_mass
    grp.create_dataset('zmid', data=zmids) 
    grp.create_dataset('slope', data=slopes) 
    grp.create_dataset('yint', data=yints) 

    f.close()

def get_sfmsfit_sfr(slope, yint, clobber=False):
    ''' Wrapper to extra slope, yint fit values
    '''
    fit_file = get_sfmsfit_slope_yint_file(slope, yint)

    if (os.path.isfile(fit_file) == False) or (clobber == True):   # file doesn't exist 
        build_sfmsfit_sfr(slope, yint) 
    
    f = h5py.File(fit_file, 'r') 

    zmids = f['slope_yint/zmid'][:]
    slopes = f['slope_yint/slope'][:]
    yints = f['slope_yint/yint'][:]
    
    f.close() 
    return [zmids, slopes, yints]

def get_sfmsfit_slope_yint_file(slope, yint): 
    ''' Given low redshift slope and y-int, get fit files
    input
    -----
    slope, yint

    output
    ------
    file name 
    '''

    slope_str = "%.2f" % slope
    yint_str = "%.2f" % yint 
    output_file = ''.join(['/data1/hahn/central_quenching/sf_ms/', 
        'sf_ms_fit_slope', slope_str, '_yint', yint_str, '.hdf5']) 

    return output_file 

# Group catalog SF-MS ----------------
def build_groupcat_sf(Mrcut=18): 
    ''' Build SF population for the SDSS group catalog for group catalog with specified Mrcut 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 

    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'

    # import centrals  
    centrals = cq_group.central_catalog(Mrcut=Mrcut, clobber=True) 
    
    # classification motivated by Salim et al. 
    sf_gals = centrals.sfr > -1.30 + 0.65*(centrals.mass-10.0)

    centrals_sfms = cq.CenQue() 

    group_catalog_columns = ['mass', 'sfr', 'ssfr']
    for column in group_catalog_columns: 
        column_data = getattr(centrals, column)[sf_gals]
        setattr(centrals_sfms, column, column_data) 

    output_file = ''.join(['/data1/hahn/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.starforming.hdf5']) 
    centrals_sfms.writeout(columns=group_catalog_columns,
            input_file=output_file) 

def sf_centrals(Mrcut=18, clobber=False): 
    ''' Read SDSS star-forming central group catalog into CenQue class

    Parameters
    ----------
    Mrcut : Absolute mangitude cut that specifies the group catalog 
    clobber : If True, re-construct the catalog. If False, just read catalog 

    Notes
    -----
    SF determined by a variation of the Salim et al. equation from Moustakas et al. 2013

    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    catalog_file = ''.join(['/data1/hahn/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.starforming.hdf5']) 
    if (os.path.isfile(catalog_file) == False) or (clobber == True): 
        build_groupcat_sf(Mrcut=Mrcut, central=True)
    
    f = h5py.File(catalog_file, 'r') 
    grp = f['cenque_data']
    mass = grp['mass'][:]
    sfr = grp['sfr'][:]
    ssfr = grp['ssfr'][:]

    catalog = cq.CenQue() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    
    f.close()
    return catalog 

def get_sfr_mstar_z_groupcat(m_star, Mrcut=18, clobber=False): 
    ''' Return SFR(m_star, z_in) from SDSS Group Catalog SFMS 

    Parameters 
    ----------
    m_star: stellar mass 
    Mrcut : absolute magnitude cut that specifies the group catalog
    clobber : If True,____. If False 
    
    Return
    ------
    sfr_out: star formation rate
    '''
    sf_data = sf_centrals(Mrcut=Mrcut, clobber=clobber) 
    
    # determine mass bin and redshift bin to splice the data
    massbins = util.simple_mass_bin()
    mbin_index = (np.array(massbins.mass_low) <= m_star) & \
            (np.array(massbins.mass_high) > m_star) 
    mass_low = (np.array(massbins.mass_low)[mbin_index])[0]
    mass_high = (np.array(massbins.mass_high)[mbin_index])[0]

    # splice data
    bin_index = (sf_data.mass >= mass_low) & \
            (sf_data.mass < mass_high)

    n_gal = len(sf_data.mass[bin_index])

    if len(sf_data.sfr[bin_index]) > 0: 
        avg_sfr = np.average(sf_data.sfr[bin_index]) 
        var_sfr = np.average( (sf_data.sfr[bin_index] - avg_sfr)**2 ) 
    else: 
        avg_sfr = -999.
        var_sfr = -999.

    return [avg_sfr, var_sfr, n_gal] 

def get_bestfit_groupcat_sfms(Mrcut=18, clobber=False):
    ''' Returns parameters for the best-fit line of the Star-forming SDSS Group Catalog (SF-MS)

    Parameters 
    ----------
    Mrcut : absolute magntiude cut that specifies the group catalog
    clobber : Rewrite if True

    Notes
    -----
    * Uses the SF SDSS Group Catalog 
    * Bestfit values are accessed from file,  unless it doesn't exist or clobber == True 

    '''
    fid_mass = 10.5         # fiducial mass 

    save_file = ''.join(['/data1/hahn/central_quenching/sf_ms/'
        'sf_ms_fit_starforming_groupcat.hdf5']) 
    
    if (os.path.isfile(save_file) == False) or (clobber == True): 
        # if file doesn't exist save to hdf5 file 
        f = h5py.File(save_file, 'w') 
        grp = f.create_group('slope_yint')      # slope-yint group
    
        grp.attrs['fid_mass'] = fid_mass    # fid mass meta data 

        avg_sfrs, var_sfrs, masses = [], [], [] 
        for mass in np.arange(9.5, 11.5, 0.25): 
            avg_sfr, var_sfr, ngal = get_sfr_mstar_z_groupcat(mass, Mrcut=Mrcut)
            
            if ngal < 100: 
                continue 

            masses.append(mass)
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)

        p0 = [0.5607, 0.0775917]
        fa = {'x': np.array(masses)-10.5, 'y': np.array(avg_sfrs)}
        bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, nprint=0) 
        
        # save to file 
        grp.create_dataset('zmid', data=[0.1]) 
        grp.create_dataset('slope', data=[bestfit.params[0]]) 
        grp.create_dataset('yint', data=[bestfit.params[1]]) 
        
        return [bestfit.params[0], bestfit.params[1]]
    else: 
        # if file already exists then just access the numbers from file 
        f = h5py.File(save_file, 'r') 

        return [f['slope_yint/slope'][:], f['slope_yint/yint'][:]] 

def build_groupcat_q(Mrcut=18): 
    ''' Build Q population for the SDSS group catalog for group catalog with specified Mrcut 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 

    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'

    # import centrals  
    centrals = cq_group.central_catalog(Mrcut=Mrcut, clobber=True) 
    
    # classification motivated by Salim et al. 
    sf_gals = centrals.sfr < -1.30 + 0.65*(centrals.mass-10.0)

    centrals_sfms = cq.CenQue() 

    group_catalog_columns = ['mass', 'sfr', 'ssfr']
    for column in group_catalog_columns: 
        column_data = getattr(centrals, column)[sf_gals]
        setattr(centrals_sfms, column, column_data) 

    output_file = ''.join(['/data1/hahn/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.quiescent.hdf5']) 
    centrals_sfms.writeout(columns=group_catalog_columns,
            input_file=output_file) 

def q_centrals(Mrcut=18, clobber=False): 
    ''' Read SDSS quiescent central group catalog into CenQue class

    Parameters
    ----------
    Mrcut : Absolute mangitude cut that specifies the group catalog 
    clobber : If True, re-construct the catalog. If False, just read catalog 

    Notes
    -----
    Q determined by a variation of the Salim et al. equation from Moustakas et al. 2013

    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    catalog_file = ''.join(['/data1/hahn/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.quiescent.hdf5']) 

    if (os.path.isfile(catalog_file) == False) or (clobber == True): 
        # (re)construct quiescent group catalog 
        build_groupcat_q(Mrcut=Mrcut)
    
    f = h5py.File(catalog_file, 'r') 
    grp = f['cenque_data']
    mass = grp['mass'][:]
    sfr = grp['sfr'][:]
    ssfr = grp['ssfr'][:]

    catalog = cq.CenQue() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    
    f.close()
    return catalog 

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
    massbins = util.simple_mass_bin()
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

def get_bestfit_qgroupcat_ssfr(Mrcut=18, clobber=False):
    ''' Returns parameters for the best-fit line of the mass vs sSFR relation of the 
    quiescent SDSS Group Catalog

    Parameters 
    ----------
    Mrcut : absolute magntiude cut that specifies the group catalog
    clobber : Rewrite if True

    Notes
    -----
    * Uses the Q SDSS Group Catalog 
    * Bestfit values are accessed from file,  unless it doesn't exist or clobber == True 

    '''
    fid_mass = 10.5         # fiducial mass 

    save_file = ''.join(['/data1/hahn/central_quenching/sf_ms/'
        'ssfr_mass_fit_quiescent_groupcat.hdf5']) 
    
    if (os.path.isfile(save_file) == False) or (clobber == True): 
        # if file doesn't exist save to hdf5 file 
        f = h5py.File(save_file, 'w') 
        grp = f.create_group('slope_yint')      # slope-yint group
    
        grp.attrs['fid_mass'] = fid_mass    # fid mass meta data 

        med_ssfrs, var_ssfrs, masses = [], [], [] 
        for mass in np.arange(9.5, 11.5, 0.25): 
            med_ssfr, var_ssfr, ngal = get_ssfr_mstar_qgroupcat(mass, Mrcut=Mrcut)
            
            if ngal < 10: 
                continue 

            masses.append(mass)
            med_ssfrs.append(med_ssfr)
            var_ssfrs.append(var_ssfr)

        p0 = [0.5, -12.0]
        fa = {'x': np.array(masses)-10.5, 'y': np.array(med_ssfrs)}
        bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, nprint=0) 
        
        # save to file 
        grp.create_dataset('zmid', data=[0.1]) 
        grp.create_dataset('slope', data=[bestfit.params[0].item()]) 
        grp.create_dataset('yint', data=[bestfit.params[1].item()]) 
        
        return [bestfit.params[0], bestfit.params[1]]
    else: 
        # if file already exists then just access the numbers from file 
        f = h5py.File(save_file, 'r') 

        return [f['slope_yint/slope'][:], f['slope_yint/yint'][:]] 

# EnvCount SF-MS --------------------
def get_bestfit_envcount_sfms(clobber=False):
    ''' Returns parameters for the best-fit line of the Star-forming Envcount (SF-MS)

    Parameters 
    ----------
    clobber : If True, reconstruct file 

    Returns 
    -------
    [ slope, yint ] : slope and yint of SF-MS bestfit line of SDSS data 

    Notes
    -----
    * Uses Star-forming sample of QFENV project
    * Reads from file if available. Otherwise constructs file

    '''
    fid_mass = 10.5         # fiducial mass 

    save_file = ''.join(['/data1/hahn/central_quenching/sf_ms/'
        'sf_ms_fit_starforming_envcount.hdf5']) 

    if (os.path.isfile(save_file) == False) or (clobber == True): 
        # if file doesn't exist save to hdf5 file 
        f = h5py.File(save_file, 'w') 
        grp = f.create_group('slope_yint')      # slope-yint group
    
        grp.attrs['fid_mass'] = fid_mass    # fid mass meta data 

        avg_sfrs, var_sfrs, masses = [], [], [] 
        for mass in np.arange(9.0, 11.5, 0.25): 
            avg_sfr, var_sfr, ngal = get_sfr_mstar_z_envcount(mass, 0.1)
            
            if ngal < 100: 
                continue 

            masses.append(mass)
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)

        p0 = [0.5607, 0.0775917]
        fa = {'x': np.array(masses)-10.5, 'y': np.array(avg_sfrs)}
        bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, nprint=0) 
        
        # save to file 
        grp.create_dataset('zmid', data=[0.1]) 
        grp.create_dataset('slope', data=[bestfit.params[0]]) 
        grp.create_dataset('yint', data=[bestfit.params[1]]) 
        
        return [bestfit.params[0], bestfit.params[1]]
    else: 
        # if file already exists then just access the numbers from file 
        f = h5py.File(save_file, 'r') 

        return [f['slope_yint/slope'][:], f['slope_yint/yint'][:]] 

if __name__=='__main__':
    print get_bestfit_groupcat_sfms(Mrcut=18) 
    print get_bestfit_groupcat_sfms(Mrcut=19) 
    print get_bestfit_groupcat_sfms(Mrcut=20) 
    #build_groupcat_sf(Mrcut=18)
    #build_groupcat_sf(Mrcut=19)
    #build_groupcat_sf(Mrcut=20)

    # SDSS group catalog best fit 
    #groupcat_slope, groupcat_yint = sdss_groupcat_sfms_bestfit()
    #print get_sfmsfit_sfr(groupcat_slope, groupcat_yint)
