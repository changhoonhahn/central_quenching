'''

Fit SF-MS


'''
import h5py
# --- Local ---
import mpfit
from sf_mainseq import get_sfr_mstar_z_groupcat
from sf_mainseq import get_sfr_mstar_z_envcount

def build_bestfit_sfms_groupcat(Mrcut=18, clobber=False):
    ''' Calculate the linear bestfit for the StarForming Main Sequence of the 
    SDSS Group Catalog specified by Mrcut

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

    bestfit_file = ''.join(['dat/central_quenching/sf_ms/'
        'sf_ms_fit_starforming_groupcat_', str(Mrcut), '.hdf5']) 
    
    avg_sfrs, var_sfrs, masses = [], [], [] 

    for mass in np.arange(9.5, 11.5, 0.25): 
        avg_sfr, var_sfr, ngal = get_sfr_mstar_z_groupcat(mass, Mrcut=Mrcut)
        
        if ngal < 500: 
            continue 

        masses.append(mass)
        avg_sfrs.append(avg_sfr)
        var_sfrs.append(var_sfr)
    
    p0 = [0.5607, 0.0775917]    # guess 
    fa = {
            'x': np.array(masses)-10.5, 
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

def build_bestfit_sfms_envcount(clobber=False):
    ''' Calculate the linear bestfit for the StarForming Main Sequence of the 
    QfEnv project
    
    ----------------------------------------------------------------
    Notes
    ----------------------------------------------------------------
    * Uses Star-forming sample of QFENV project
    * Reads from file if available. Otherwise constructs file

    '''

    fid_mass = 10.5         # fiducial mass 
    
    # average SFR, sigma SFR for array of masses 
    avg_sfrs, var_sfrs, masses = [], [], [] 

    for mass in np.arange(9.0, 11.5, 0.25): 

        avg_sfr, var_sfr, ngal = get_sfr_mstar_z_envcount(mass, 0.1)
        
        if ngal < 100: 
            continue 

        masses.append(mass)
        avg_sfrs.append(avg_sfr)
        var_sfrs.append(var_sfr)

    p0 = [0.5607, 0.0775917]
    fa = {
            'x': np.array(masses) - fidmass, 
            'y': np.array(avg_sfrs), 
            'err': np.array(var_sfrs)
            }
    bestfit = mpfit.mpfit(
            mpfit_line, 
            p0, 
            functkw=fa
            ) 
    
    # save to hdf5 file using h5py
    bestfit_file = ''.join([
        'dat/central_quenching/sf_ms/'
        'sf_ms_fit_starforming_envcount.hdf5'
        ]) 
    f = h5py.File(save_file, 'w') 
    grp = f.create_group('slope_yint')      # slope-yint group
    grp.attrs['fid_mass'] = fid_mass    # fid mass meta data 
    grp.create_dataset('zmid', data=[0.1]) 
    grp.create_dataset('slope', data=[bestfit.params[0]]) 
    grp.create_dataset('yint', data=[bestfit.params[1]]) 
        
def line(x, p): 
    return p[0]*x + p[1]

def mpfit_line(p, fjac=None, x=None, y=None, err=None): 
    model = line(x, p) 
    status = 0 

    return([status, (y-model)/err]) 

if __name__=="__main__":
    build_bestfit_sfms_groupcat(Mrcut=18, clobber=True)
    build_bestfit_sfms_envcount(clobber=True)
