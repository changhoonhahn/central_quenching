'''

Utility functions for Central Quenching project 

Author(s): ChangHoon Hahn

'''
import subprocess
import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import os 

import mpfit
import h5py
from scipy import signal
from scipy import interpolate

# ---- Local -----

# direcotry stuff
def code_dir(): 
    return os.path.dirname(os.path.realpath(__file__)).split('util')[0]

def intersection_index(arr1, arr2):  
    """ 
    Find the indicies of the intersecting elements of arr1 and arr2.
    Takes approximately < 1 second
    """
    sort_arr1_indices = np.argsort(arr1)
    sort_arr2_indices = np.argsort(arr2)

    sorted_arr1 = arr1[sort_arr1_indices]
    sorted_arr2 = arr2[sort_arr2_indices]

    arr1_in1d = np.in1d(sorted_arr1, sorted_arr2)
    arr2_in1d = np.in1d(sorted_arr2, sorted_arr1)

    arr1_intersect_indices = sort_arr1_indices[arr1_in1d]
    arr2_intersect_indices = sort_arr2_indices[arr2_in1d]

    return arr1_intersect_indices, arr2_intersect_indices 

# SF property functions ------------------------------------------------------------------------
def get_sfr_mstar_z_flex(mstar, z_in, built_sfms_fit):
    ''' given built SFMS fits from sf_mainseq function build_sfr_mstar_z 
    output SFR of mstar and z_in 
    
    Parameters
    ----------
    mstar : 
    z_in : 
    built_sfms_fit: (z_mid, slope, yint) 

    '''
    fid_mass = 10.5 
    
    SFR_amp = np.interp(z_in, np.array(built_sfms_fit[0]), np.array(built_sfms_fit[2])) 

    zmids = built_sfms_fit[0]

    # closest redshift index
    closest_i = min(range(len(zmids)), key=lambda i: abs(zmids[i] - z_in))

    # assuming slope doesn't change calculate average SFR
    avg_sfr = (built_sfms_fit[1])[closest_i] * (mstar - fid_mass) + SFR_amp
        
    return [avg_sfr, 0.28]       # 0.3 dex scatter hard-coded

def get_sfr_mstar_z(mstar, z_in, deltamass=0.2, deltaz=0.2, lit='primusfit', machine='harmattan'):
    ''' Get SFR using SF main sequence as a function of mass and redshift
    outputs [average SFR, standard deviation SFR] for mass and redshift bin 
    '''
    sig_mass = 0.5 * deltamass
    sig_z = 0.5 * deltaz

    if lit == 'primusfit':     
        ''' SFR value determined from linear fits to the PRIMUS data SF main sequence

        '''
        # import SF-MainSequence best-fit 
        sf_ms_file = ''.join(['dat/wetzel_tree/envcount/', 
            'sfr_mstar_fit_param_fixedslope.fits']) 
        sf_ms = mrdfits(sf_ms_file) 
        
        fid_mass = 10.5         # fid mass (part of the fititng) 
        
        sf_ms.yint[sf_ms.z == 0.1] = 0.134      # totally unjustified hack
        SFR_amp = np.interp(z_in, sf_ms.z, sf_ms.yint)  # interpolate SFR amplitude by z input

        # closest redshift index
        closest_i = min(range(len(sf_ms.z)), key=lambda i: abs(sf_ms.z[i] - z_in))
    
        # assuming slope doesn't change calculate average SFR
        avg_sfr = (sf_ms.slope)[closest_i] * (mstar - fid_mass) + SFR_amp

        return [avg_sfr, 0.3]       # 0.3 dex scatter hard-coded

    elif lit == 'primus':       # using PRIMUS central galaxy data 
        
        if machine == 'harmattan': 
            data_dir = 'dat/wetzel_tree/envcount/'
        else: 
            data_dir = '/global/data/scr/chh327/primus/data/envcount/'
        # iSEDfit galaxy set 
        sf_galaxy_file = ''.join([data_dir, 
            'envcount_cylr2h25_thresh75_active_lit_EDP-primus-z0210-numden.fits'
            ]) 
        sf_galaxy = mrdfits(sf_galaxy_file) 

        mass_z_bin = (sf_galaxy.mass > mstar - sig_mass) & \
                (sf_galaxy.mass <= mstar + sig_mass) &\
                (sf_galaxy.z > z_in - sig_z) & \
                (sf_galaxy.z <= z_in + sig_z) 

        bin_sfr = sf_galaxy.sfr[mass_z_bin]     # SFRs of SF galaxies within the bin 
        
        avg_sfr = np.mean(bin_sfr) 
        std_sfr = np.std(bin_sfr) 

        return [avg_sfr, std_sfr]
    else: 
        raise NameError("Not yet coded") 

def get_q_ssfr_mean(masses, Mrcut=18): 
    ''' Return average/median sSFR of quiscent population of the SDSS Group Catalog for an
    array of mass 

    Parameters 
    ----------
    masses : array of masses 
    Mrcut : 

    Notes
    -----
    * A little bit of hardcoded tweaking but these should not matter because ultimately  

    '''
    if isinstance(masses, list): 
        masses = np.array(masses)

    fit_line_param = get_bestfit_qgroupcat_ssfr(Mrcut=Mrcut) 
    fit_slope = fit_line_param[0].item() 
    fit_yint = fit_line_param[1].item() 

    #q_ssfr = (fit_slope + 0.05) * (masses - 10.5) + fit_yint# - 0.1  
    #q_ssfr = (fit_slope + 0.05) * (masses - 10.4) + fit_yint  
    q_ssfr = (fit_slope + 0.15) * (masses - 10.4) + fit_yint - 0.1
    
    mass = np.arange(9.5, 11.5, 0.1)
    #print fit_slope, fit_yint
    #print (fit_slope + 0.05) * (mass - 10.4) + fit_yint
    #print (fit_slope + 0.125) * (mass - 10.4) + fit_yint - 0.075
    #print (fit_slope + 0.15) * (mass - 10.4) + fit_yint - 0.1

    #q_ssfr = np.array([ (-0.7 * mass) - 4.625 for mass in masses ])  
    return q_ssfr 


# integrated mass ---------------------------------------------
def sfr_avg_residual(mass, z_cosmic, **sfr_param): 
    ''' SFR = average SFR + sampled residual  

    Parameters
    ----------
    mass : stellar mass 
    z_cosmic : redshift 
    sfr_param : SFR parameters  (residual, 'resid') 

    Notes
    -----

    '''
    if 'resid' in sfr_param.keys(): 
        sfr_resid = sfr_param['resid']
    else: 
        raise NameError('asdflkjaskdf') 

    sfr, sig_sfr = get_sfr_mstar_z_bestfit( mass, z_cosmic, Mrcut=18)
    #print 'SF AVG RESIDUAL function ', sfr

    return sfr + sfr_resid 

def sfr_squarewave(mass, tcosmic, **sfr_param): 
    ''' SFR determined by square function 

    Parameters
    ----------
    mass : stellar mass 
    tcosmic : cosmic time 
    sfr_param : SFR parameters  (amplitude, phase, and frequency) 

    Notes
    -----
    * Arbitrary sfr_param will be fed to produce random SFRs
    * 

    '''
    sfr_A = sfr_param['amp']
    sfr_d = sfr_param['phase']
    sfr_w = sfr_param['freq']

    if not isinstance(sfr_w, float): 
        if len(sfr_w) != len(mass): 
            return None

    dSFR = sfr_A * signal.square(sfr_w * (tcosmic - 6.9048 + sfr_d))    # 
     
    z_cosmic = get_zsnap(tcosmic) 

    sfr, sig_sfr = get_sfr_mstar_z_bestfit( mass, z_cosmic, Mrcut=18)

    return sfr + dSFR

# z_snap <--> t_snap ---------------------------------------------

def snapshottable():
    '''
    '''
    return os.path.dirname(os.path.realpath(__file__)).split('CenQue')[0]+'dat/snapshot_table.dat'

def get_zsnap(tcosmic): 
    ''' Given cosmic time return redshift using spline interpolation of snapshot table 

    Parameters
    ----------
    tcosmic : cosmic time 

    Notes
    -----
    * Only worry is that it may take too long
    '''
    # read in snapshot table file 
    z = np.array([0.0000, 0.0502, 0.1028, 0.1581, 0.2162, 0.2771, 0.3412, 0.4085, 0.4793, 0.5533, 0.6313, 0.7132, 0.7989, 0.8893, 0.9841, 1.0833, 1.1882, 1.2978, 1.4131, 1.5342, 1.6610, 1.7949, 1.9343, 2.0817, 2.2362, 2.3990, 2.5689, 2.7481, 2.9370, 3.1339, 3.3403, 3.5579, 3.7870])
    t = np.array([13.8099, 13.1328, 12.4724, 11.8271, 11.1980, 10.5893, 9.9988, 9.4289, 8.8783, 8.3525, 7.8464, 7.3635, 6.9048, 6.4665, 6.0513, 5.6597, 5.2873, 4.9378, 4.6080, 4.2980, 4.0079, 3.7343, 3.4802, 3.2408, 3.0172, 2.8078, 2.6136, 2.4315, 2.2611, 2.1035, 1.9569, 1.8198, 1.6918])

    z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 

    return z_of_t(tcosmic) 

def get_tsnap(redshift): 
    ''' Given redshift, return cosmic time using spline interpolation of snapshot table

    Parameters
    ----------
    redshift : redshift 

    Notes
    -----
    * Only worry is that it may take too long

    '''
    # read in snapshot table file 
    z = np.array([0.0000, 0.0502, 0.1028, 0.1581, 0.2162, 0.2771, 0.3412, 0.4085, 0.4793, 0.5533, 0.6313, 0.7132, 0.7989, 0.8893, 0.9841, 1.0833, 1.1882, 1.2978, 1.4131, 1.5342, 1.6610, 1.7949, 1.9343, 2.0817, 2.2362, 2.3990, 2.5689, 2.7481, 2.9370, 3.1339, 3.3403, 3.5579, 3.7870])
    t = np.array([13.8099, 13.1328, 12.4724, 11.8271, 11.1980, 10.5893, 9.9988, 9.4289, 8.8783, 8.3525, 7.8464, 7.3635, 6.9048, 6.4665, 6.0513, 5.6597, 5.2873, 4.9378, 4.6080, 4.2980, 4.0079, 3.7343, 3.4802, 3.2408, 3.0172, 2.8078, 2.6136, 2.4315, 2.2611, 2.1035, 1.9569, 1.8198, 1.6918])

    t_of_z = interpolate.interp1d(z, t, kind='cubic') 

    return t_of_z(redshift) 

def get_z_nsnap(nsnap): 
    ''' Given snapshot, return redshift using snapshot table

    Parameters
    ----------
    redshift : redshift 

    Notes
    -----
    * Only worry is that it may take too long

    '''
    # read in snapshot table file 
    z = np.loadtxt(snapshottable(), unpack=True, usecols=[2]) 

    return z[nsnap]

def get_t_nsnap(nsnap): 
    ''' Given snapshot, return cosmic time using snapshot table

    Parameters
    ----------
    redshift : redshift 

    Notes
    -----
    * Only worry is that it may take too long

    '''
    # read in snapshot table file 
    t = np.loadtxt(snapshottable(), unpack=True, usecols=[3]) 

    return t[nsnap]

def get_nsnap_t(tcosmic): 
    ''' Given cosmic time, return closest snapshot number 
    '''
    t = np.array([13.8099, 13.1328, 12.4724, 11.8271, 11.1980, 10.5893, 9.9988, 9.4289, 8.8783, 8.3525, 7.8464, 7.3635, 6.9048, 6.4665, 6.0513, 5.6597, 5.2873, 4.9378, 4.6080, 4.2980, 4.0079, 3.7343, 3.4802, 3.2408, 3.0172, 2.8078, 2.6136, 2.4315, 2.2611, 2.1035, 1.9569, 1.8198, 1.6918])
    
    if isinstance(tcosmic, float): 
        snapshot = np.array([(np.abs(t - tcosmic)).argmin()])

    elif isinstance(tcosmic, np.ndarray): 
        snapshot = [] 
        for i in range(len(tcosmic)): 
            snapshot.append((np.abs(t - tcosmic[i])).argmin())

        snapshot = np.array(snapshot)
    return snapshot 

def line(x, p): 
    # just line function 
    return p[0]*x + p[1]

def line_fixedslope(x, p, slope=0.56): 
    # Line with specified fixed slope 
    return slope * x + p[0]

def mpfit_line(p, fjac=None, x=None, y=None, err=None): 
    model = line(x, p) 
    status = 0 
    if err is None: 
        err = np.repeat(1., len(model))
    return([status, (y-model)/err]) 

def mpfit_line_fixedslope(p, slope=0.56, fjac=None, x=None, y=None, err=None): 
    model = line_fixedslope(x, p, slope=slope) 
    status = 0 
    return([status, (y-model)/err]) 

def png2pdf(png_filename): 
    ''' Convert png file to pdf 
    '''
    pdf_filename = png_filename.replace('.png', '.pdf')

    convert_cmd = ' '.join(['convert', png_filename, pdf_filename])
    
    subprocess.call(convert_cmd.split())
    return None 



