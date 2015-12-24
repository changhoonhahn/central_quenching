'''

Utility functions for CenQue

Author(s): ChangHoon Hahn

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import os 
import mpfit
import h5py
from scipy import signal
from scipy import interpolate

# ---- Local -----

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

def get_bestfit_qgroupcat_ssfr(Mrcut=18, fid_mass=10.5, clobber=False):
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

    save_file = ''.join([
        'dat/central_quenching/sf_ms/'
        'ssfr_mass_fit_quiescent_groupcat.hdf5'
        ]) 
    
    if not os.path.isfile(save_file) or clobber: 
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

        p0 = [-0.5, -12.0]
        fa = {'x': np.array(masses)-10.5, 'y': np.array(med_ssfrs)}
        bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, nprint=0) 
        
        # save to file 
        grp.create_dataset('zmid', data=[0.1]) 
        grp.create_dataset('slope', data=[bestfit.params[0].item()]) 
        grp.create_dataset('yint', data=[bestfit.params[1].item()]) 
        
        return [bestfit.params[0], bestfit.params[1]]
    else: 

        f = h5py.File(save_file, 'r') 

        return [f['slope_yint/slope'][:], f['slope_yint/yint'][:]] 

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

def integrated_mass_rk4(sfr, mass, t0, tf, f_retain=0.6, **sfr_param): 
    ''' Integrated stellar mass using RK4 integration
    
    Parameters
    ----------
    sfr : SFR python function (for example sfr_squarewave)
    mass : initial stellar mass  
    t0 : initial cosmic time 
    tf : final cosmic time
    f_retain : fraction of stellar mass not lost from SNe and winds from Wetzel Paper
    sfr_param : parameters of the SFR function 

    '''

    delt = 0.025 # Gyr
    iter = int(np.round( (tf-t0)/delt )) 
    delt = (tf - t0)/np.float(iter) 
    
    t_n_1 = t0 
    SFR_n_1 = sfr(mass, t0, **sfr_param)
    M_n_1 = mass
    
    for i in range(iter): 
        t_n = t_n_1 + delt

        k1 = f_retain * (10.0 ** SFR_n_1)
        k2 = f_retain * (10.0 ** (
            sfr(np.log10(10.0**M_n_1 + 10**9 * delt/2.0 * k1), t_n_1 + delt/2.0, **sfr_param))) 
        k3 = f_retain * (10.0 ** (
            sfr(np.log10(10.0**M_n_1 + 10**9 * delt/2.0 * k2), t_n_1 + delt/2.0, **sfr_param))) 
        k4 = f_retain * (10.0 ** (
            sfr(np.log10(10.0**M_n_1 + 10**9 * delt * k3), t_n_1 + delt, **sfr_param))) 

        M_n = np.log10(10.0 ** M_n_1 + 1.0/6.0 * delt * 10**9 * (k1 + 2.0*k2 + 2.0*k3 + k4)) 
        
        if np.sum(np.isnan(M_n)) > 0: 
            raise NameError('asldkfjalkjsdflk;ajsdf') 

        SFR_n_1 = sfr(M_n, t_n, **sfr_param)

        M_n_1 = M_n
        t_n_1 = t_n
    
    return M_n

# z_snap <--> t_snap ---------------------------------------------

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
    z = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2]) 

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
    t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[3]) 

    return t[nsnap]

"""
def cenque_file( **kwargs ): 
    ''' Given kwargs get CenQue file name
    
    Parameters (try to keep this up-to-date!)
    ----------
    nsnap : snapshot number 
    file_type : "sf assign", "evol from"
    tau : "discrete", "linefit", or tau flag
    sfms_slope : slope of the SF Main Sequence
    sfms_yint : yint of the SF Main Sequence 

    Notes
    -----
    * MORE FILE TYPES WILL BE SPECIFIED

    '''

    if 'input_file' in kwargs.keys(): 
        cenque_filename = kwargs['input_file']
    else: 
        try: 
            kwargs['nsnap'] 
        except KeyError:
            raise KeyError("Specify input file or nsnap + file_type") 

        try: 
            kwargs['file_type']
        except KeyError:
            file_type_str = '' 
        else: 
            if min( (kwargs['file_type']).find('sf'), (kwargs['file_type']).find('ass') ) > -1: 
                # star-formation assign 
                file_type_str = '_sfpropassign'
                
                # Quenching Fraction specifier 
                if 'fqing_slope' in kwargs.keys(): 
                    fqing_slope_str = str("%.2f" % kwargs['fqing_slope'])
                else: 
                    fqing_slope_str = str("%.2f" % 0.63)

                if 'fqing_yint' in kwargs.keys(): 
                    fqing_yint_str = str("%.2f" % kwargs['fqing_yint'])
                else: 
                    fqing_yint_str = str("%.2f" % -6.04) 

                fqing_str = ''.join([fqing_slope_str, '_', fqing_yint_str, 'fqing']) 

                # combine specifiers
                file_type_str = ''.join(['_', fqing_str, file_type_str])
                
            elif min( (kwargs['file_type']).find('evo'), (kwargs['file_type']).find('from') ) > -1: 
                # evolved from nsnap
                original_nsnap = int(((kwargs['file_type']).split('from'))[-1]) 

                file_type_str = '_evol_from'+str(original_nsnap) 
               
                # Tau specifier
                if kwargs['tau'] == 'discrete': 
                    tau_str = '_'+'_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
                elif kwargs['tau'] == 'linefit':
                    tau_str = '_line'+'_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
                else: 
                    tau_str = '_'+kwargs['tau']+'tau'

                # Quenching Fraction specifier 
                if 'fqing_slope' in kwargs.keys(): 
                    fqing_slope_str = str(kwargs['fqing_slope'])
                else: 
                    fqing_slope_str = str(0.63)

                if 'fqing_yint' in kwargs.keys(): 
                    fqing_yint_str = str(kwargs['fqing_yint'])
                else: 
                    fqing_yint_str = str(-6.04) 

                fqing_str = '_'+fqing_slope_str+'_'+fqing_yint_str+'fqing'

                # combine specifiers
                file_type_str = ''.join([tau_str, fqing_str, file_type_str]) 

            else: 
                raise NameError("File not specified") 
    
        if 'sfms_slope' not in kwargs.keys(): 
            sfms_param_str = ''
        else: 
            slope_str = "%.2f" % kwargs['sfms_slope']
            yint_str = "%.2f" % kwargs['sfms_yint'] 
            sfms_param_str = '_sfms_slope'+slope_str+'_yint'+yint_str

        cenque_filename = ''.join(['dat/central_quenching/', 
            'cenque_centrals_snapshot', str(kwargs['nsnap']), file_type_str, sfms_param_str, 
            '.hdf5']) 

    return cenque_filename
    
"""
