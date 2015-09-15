'''

Fit SF-MS


'''
import h5py
# --- Local ---
import mpfit

from sf_mainseq import get_sfr_mstar_z_groupcat
from sf_mainseq import get_sfr_mstar_z_envcount

def get_sfr_mstar_z_bestfit(mstar, z_in, Mrcut=18, clobber=False):
    ''' Calculate average SFR and standard deviation of the SF main sequence as a 
    function of mass and redshift using linear bestfits of SFMS from the SDSS 
    group catalog and qf_env catalogs. See notes for detail.

    If 
    SFR_groupcat (M*, z ~ 0.1) = alpha * M* + beta 
    and 
    SFR_envcount (M*, z) = alpha' * M* + beta'(z)

    SFR(M*,z) = alpha * M* + beta + (beta'(z) - beta'(z ~ 0.1))

    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    mstar : Stellar mass of galaxy
    z_in : Redshift 
    Mrcut : Absolute magnitude cut that specified the group catalog 

    Returns
    -------
    [average SFR, standard deviation SFR]

    ----------------------------------------------------------------
    Notes
    ----------------------------------------------------------------
    * Best-fit SFMS y-int offsets for redshift bins determined from EnvCount project 
    * Slope and SDSS y-int determined from SF Main Sequence of Group Catalog 
    * Fiducial Mass = 10.5
    * Assumptions: 
        * The overall shifts in SFR observed in the iSEDfit sample is equivalent to that of the group catalog  

    '''

    fid_mass = 10.5

    # Best-fit slope and y-int of SF SDSS Group Catalog 
    groupcat_fit_param = build_bestfit_groupcat_sfms(Mrcut=Mrcut, clobber=clobber)

    # Best-fit slope and y-int of SF EnvCount
    envcount_fit_param = sfms.get_bestfit_envcount_sfms()
    zmids, slopes, yints = sfms.get_sfmsfit_sfr(
            (envcount_fit_param[0]).item(), (envcount_fit_param[1]).item(), 
            clobber=clobber)
        
    d_yints = np.interp(z_in, zmids, yints) - yints[0] 
    SFR_amp = groupcat_fit_param[1] + d_yints

    avg_SFR = groupcat_fit_param[0] * (mstar - fid_mass) + SFR_amp
        
    return [avg_SFR, 0.3] 

def get_bestfit_sfms_param(catalog, clobber=False):
    ''' Calculate the linear bestfit parameters for the StarForming 
    Main Sequence of the SDSS Group Catalog specified by Mrcut or 
    SDSS+PRIMUS envcount catalog from qf_env project. Fitting is done
    using MPFit
    
    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    Mrcut : absolute magntiude cut that specifies the group catalog
    clobber : Rewrite if True
    ----------------------------------------------------------------
    Notes
    ----------------------------------------------------------------
    * Bestfit values are accessed from file,  unless it doesn't exist or clobber == True 

    '''
    fid_mass = 10.5         # fiducial mass 
    
    if catalog == 'groupcat': 
        bestfit_file = ''.join([
            'dat/central_quenching/sf_ms/'
            'sf_ms_fit_starforming_groupcat_', str(Mrcut), '.hdf5'
            ]) 
    elif catalog == 'envcount': 
        bestfit_file = ''.join([
            'dat/central_quenching/sf_ms/'
            'sf_ms_fit_starforming_envcount.hdf5'
            ]) 

    if not os.path.isfile(bestfit_file): 

        avg_sfrs, var_sfrs, masses = [], [], [] 

        for mass in np.arange(9.0, 11.5, 0.25): 
            
            if catalog == 'groupcat': 
                avg_sfr, var_sfr, ngal = get_sfr_mstar_z_groupcat(mass, Mrcut=Mrcut)
            elif catalog == 'envcount': 
                avg_sfr, var_sfr, ngal = get_sfr_mstar_z_envcount(mass, 0.1)
            
            if ngal < 100: 
                continue 

            masses.append(mass)
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)
        
        p0 = [0.5607, 0.0775917]    # guess 
        fa = {
                'x': np.array(masses) - fid_mass, 
                'y': np.array(avg_sfrs), 
                'err': np.array(var_sfrs)
                }
        bestfit = mpfit.mpfit(
                mpfit_line, 
                p0, 
                functkw = fa
                ) 
        
        # save data to h5py file 
        f = h5py.File(bestfit_file, 'w') 
        grp = f.create_group('slope_yint')
        grp.attrs['fid_mass'] = fid_mass 

        grp.create_dataset('zmid', data=[0.1]) 
        grp.create_dataset('slope', data=[bestfit.params[0]]) 
        grp.create_dataset('yint', data=[bestfit.params[1]]) 
       
        f.close() 

        return None 

def line(x, p): 
    return p[0]*x + p[1]

def mpfit_line(p, fjac=None, x=None, y=None, err=None): 
    model = line(x, p) 
    status = 0 

    return([status, (y-model)/err]) 

if __name__=="__main__":
    build_bestfit_sfms_groupcat(Mrcut=18, clobber=True)
    build_bestfit_sfms_envcount(clobber=True)
