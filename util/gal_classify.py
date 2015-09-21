'''

Code for classifying galaxies into 
star-forming or quiescent based on 
their SF properties and stellar mass 

'''
import numpy as np

def sfq_classify(mstar, sfr, z_in, Mrcut=18, clobber=False):
    ''' Classify galaxies' SFRs and stellar masses into
    star-forming or quiescent based on evolve sSFR cut.

    Parameters
    ----------
    mstar : Stellar mass of galaxy (array)
    sfr : Star-formation rate of galaxies (array)
    z_in : Redshift 
    --Mrcut : Absolute magnitude cut that specified the group catalog --

    Returns
    -------
    [average SFR, standard deviation SFR]

    Notes
    -----
    * Best-fit SFMS y-int offsets for redshift bins determined from EnvCount project 
    * Slope and SDSS y-int determined from SF Main Sequence of Group Catalog 
    * Fiducial Mass = 10.5
    * Assumptions: 
        * The overall shifts in SFR observed in the iSEDfit sample is equivalent to that of the group catalog  

    '''
    
    #fid_mass = 10.5

    ## Best-fit slope and y-int of SF SDSS Group Catalog 
    #groupcat_fit_param = sfms.get_bestfit_groupcat_sfms(Mrcut=Mrcut, clobber=clobber)

    ## Best-fit slope and y-int of SF EnvCount
    #envcount_fit_param = sfms.get_bestfit_envcount_sfms()
    #zmids, slopes, yints = sfms.get_sfmsfit_sfr(
    #        (envcount_fit_param[0]).item(), (envcount_fit_param[1]).item(), 
    #        clobber=clobber)
    #    
    #d_yints = np.interp(z_in, zmids, yints) - yints[0] 
    #SFR_amp = groupcat_fit_param[1] + d_yints

    #SFR_cut = groupcat_fit_param[0] * (mstar - fid_mass) + SFR_amp - 1.0
    #ssfr_cut = -11.35 + 0.76*(z_in-0.05) - 0.49*(mstar-10.5)
    #ssfr_cut = -11.35 + 0.76*(z_in-0.05) - 0.35*(mstar-10.5)
    sfr_class = sfr_cut(mstar, z_in)

    sfq = np.empty(len(mstar), dtype=(str,16))
    sf_index = np.where(sfr > sfr_class)
    sfq[sf_index] = 'star-forming'
    q_index = np.where(sfr <= sfr_class)
    sfq[q_index] = 'quiescent'

    return sfq 

def sfr_cut(mstar, zin):
    """ Specific SFR cut off used to classify SF or Quiescent 
    galaxies 
    """
    return -0.75 + 0.76*(zin-0.05) + 0.65*(mstar-10.5)
