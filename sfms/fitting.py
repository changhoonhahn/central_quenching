'''

Fit SF-MS


'''
import os 
import h5py
import numpy as np
# --- Local ---
import mpfit

from sf_mainseq import get_sfr_mstar_z_groupcat
from sf_mainseq import get_sfr_mstar_z_envcount

def get_bestfit_sfr_mstar_z(mstar, z_in, Mrcut=18, fid_mass=10.5, clobber=False):
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

    ----------------------------------------------------------------
    Returns
    ----------------------------------------------------------------
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

    # Best-fit slope and y-int of SF SDSS Group Catalog and EnvCount catalog 
    gc_zmid, gc_slope, gc_yint = get_bestfit_sfms_groupcat(Mrcut=Mrcut, fid_mass=fid_mass, clobber=clobber)
    ec_zmid, ec_slope, ec_yint = get_bestfit_sfms_envcount(fid_mass=fid_mass, clobber=clobber)
    print gc_zmid, gc_slope, gc_yint
    print ec_zmid, ec_slope, ec_yint

    delta_yint = np.interp(z_in, ec_zmid, ec_yint) - ec_yint[0] 

    avg_SFR = gc_slope * (mstar - fid_mass) + gc_yint + delta_yint 
        
    return [avg_SFR, 0.3] 

def get_bestfit_sfms_groupcat(Mrcut=18, fid_mass=10.5, clobber=False):
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
    
    bestfit_file = ''.join([
        'dat/central_quenching/sf_ms/'
        'sf_ms_fit_starforming_groupcat_', str(Mrcut), '.hdf5'
        ]) 

    if not os.path.isfile(bestfit_file) or clobber: 

        print 'Writing '
        print bestfit_file

        avg_sfrs, var_sfrs, ngal = get_sfr_mstar_z_groupcat(
                np.arange(9.0, 11.5, 0.25), 
                Mrcut=Mrcut
                )
            
        enough_gal = np.where(ngal > 100)
            
        masses = np.arange(9.0, 11.5, 0.25 )[enough_gal]
        avg_sfrs = np.array(avg_sfrs)[enough_gal]
        var_sfrs = np.array(var_sfrs)[enough_gal]
        
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
        grp.create_dataset('slope', data=[0.65])
        grp.create_dataset('yint', data=[0.0]) 
        #grp.create_dataset('slope', data=[bestfit.params[0]]) 
        #grp.create_dataset('yint', data=[bestfit.params[1]]) 
       
        f.close() 

        #return [0.1, bestfit.params[0], bestfit.params[1]]
        return [0.1, 0.65, 0.0]

    else: 
        f = h5py.File(bestfit_file, 'r') 

        zmids = f['slope_yint/zmid'][:]
        slopes = f['slope_yint/slope'][:]
        yints = f['slope_yint/yint'][:]
        
        f.close() 

        return [zmids, slopes, yints]

def get_bestfit_sfms_envcount(fid_mass = 10.5, clobber = False):
    ''' Calculate linear bestfit parameters for SF-MS fits for redshift 
    bins with delta z = 0.2. The slope and y-int are fit for z ~ 0.1. 
    While only y-int is fit for z > 0.2 
    
    '''

    bestfit_file = ''.join([
        'dat/central_quenching/sf_ms/'
        'sf_ms_fit_starforming_envcount.hdf5'
        ]) 
    
    if not os.path.isfile(bestfit_file) or clobber:
        print 'Writing '
        print bestfit_file

        # first calculate average SFR, sigma_sfr as a function of mass 
        # for SDSS redshift bin of the envcount catalog. Afterwards 
        # fit slope and yint to the SF-MS. Notes that mass range is higher
        # than later in the function. 
        avg_sfrs, sig_sfrs, ngal = get_sfr_mstar_z_envcount(
                np.arange(9.0, 11.5, 0.25), 
                [0.1 for i in xrange(len(np.arange(9.0, 11.5, 0.25)))]
                )
            
        enough_gal = np.where(np.array(avg_sfrs) != -10.)
            
        masses = np.arange(9.0, 11.5, 0.25 )[enough_gal]
        avg_sfrs = np.array(avg_sfrs)[enough_gal]
        sig_sfrs = np.array(sig_sfrs)[enough_gal]
        
        p0 = [0.0775917]    # guess 
        fa = {
                'x': np.array(masses) - fid_mass, 
                'y': np.array(avg_sfrs), 
                'err': np.array(sig_sfrs)
                }
        bestfit = mpfit.mpfit(
                mpfit_line_fixedslope, 
                p0, 
                functkw = fa
                ) 

        sdsszbin_slope = 0.65 #bestfit.params[0]
        sdsszbin_yint = 0.2 #bestfit.params[0]
        
        # use bestfit slope of low z SDSS bin to fit fixed slope 
        # lines to the rest of the redshift bins
        zmids, slopes, yints  = [], [], [] 
        zbins = [0.3, 0.5, 0.7, 0.9] 

        for zbin in zbins: 
            
            avg_sfrs, sig_sfrs, ngal = get_sfr_mstar_z_envcount(
                    np.arange(9.0, 11.5, 0.25), 
                    [zbin for i in xrange(len(np.arange(9.0, 11.5, 0.25)))]
                    )
                
            enough_gal = np.where(np.array(avg_sfrs) > -10.)
                
            masses = np.arange(9.0, 11.5, 0.25 )[enough_gal]
            avg_sfrs = np.array(avg_sfrs)[enough_gal]
            sig_sfrs = np.array(sig_sfrs)[enough_gal]

            p0 = [0.5]  # bad guess
            fa = {
                    'slope': sdsszbin_slope, 
                    'x': np.array(masses) - fid_mass, 
                    'y': np.array(avg_sfrs), 
                    'err': np.array(sig_sfrs)
                    }

            bestfit = mpfit.mpfit(
                    mpfit_line_fixedslope, 
                    p0, 
                    functkw = fa, 
                    quiet=1 
                    )
            
            zmids.append(zbin) 
            slopes.append(sdsszbin_slope) 
            yints.append(bestfit.params[0]) 

        zmids.insert(0, 0.1)
        slopes.insert(0, sdsszbin_slope)
        yints.insert(0, sdsszbin_yint)
        
        f = h5py.File(bestfit_file, 'w') 
        grp = f.create_group('slope_yint')
        
        grp.attrs['fid_mass'] = fid_mass
        grp.create_dataset('zmid', data=zmids) 
        grp.create_dataset('slope', data=slopes) 
        grp.create_dataset('yint', data=yints) 

        f.close()
        return [zmids, slopes, yints]

    else:
        f = h5py.File(bestfit_file, 'r') 

        zmids = f['slope_yint/zmid'][:]
        slopes = f['slope_yint/slope'][:]
        yints = f['slope_yint/yint'][:]
        
        f.close() 

        return [zmids, slopes, yints]

# ---- MPfit stuff ----
def line(x, p): 
    # just line function 
    return p[0]*x + p[1]

def line_fixedslope(x, p, slope=0.56): 
    # Line with specified fixed slope 
    return slope * x + p[0]

def mpfit_line(p, fjac=None, x=None, y=None, err=None): 
    model = line(x, p) 
    status = 0 
    return([status, (y-model)/err]) 

def mpfit_line_fixedslope(p, slope=0.56, fjac=None, x=None, y=None, err=None): 
    model = line_fixedslope(x, p, slope=slope) 
    status = 0 
    return([status, (y-model)/err]) 

if __name__=="__main__":
    print get_bestfit_sfr_mstar_z(10.5, 0.7, clobber=True)
