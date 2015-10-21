"""

Code to analyze the quiecsent fraction 


"""
import numpy as np

def get_fq(Mstar, z_in, lit='cosmosinterp'): 
    ''' Calculate the quiescent fraction as a funcnction of 
    stellar mass and redshift. Different methods of calculating 
    quiescent fraction is implemented. 

    f_Q ( M_star, z) 

    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    Mstar : stellar mass
    z_in : redshift 
    lit : 'cosmosinterp', 

    '''

    if lit == 'cosmosinterp': 
        zbins = [0.36, 0.66, 0.88] 

        fq_z = [] 
        for zbin in zbins: 
            fq_file = ''.join([ 'dat/wetzel_tree/', 
                'qf_z', str(zbin), 'cen.dat' ]) 
           
            # read in mass and quiescent fraction
            mass, fq = np.loadtxt(fq_file, unpack=True, usecols=[0,1])  
    
            fq_z.append( np.interp(Mstar, mass, fq) )   # interpolate to get fq(Mstar)
        
        return np.interp(z_in, zbins, fq_z) 

    elif lit == 'cosmosfit': 
        zbins = [0.36, 0.66, 0.88] 
        exp_sigma = [1.1972271, 1.05830526, 0.9182575] 
        exp_sig = np.interp(z_in, zbins, exp_sigma) 
        output = np.exp( ( Mstar - 12.0 )/exp_sig)
        if Mstar > 12.0: 
            output = 1.0

        return output

    elif lit == 'wetzel':       # Wetzel et al. 2013
        qf_z0 = -6.04 + 0.63*Mstar

        if Mstar < 9.5: 
            alpha = -2.3
        elif (Mstar >= 9.5) & (Mstar < 10.0): 
            alpha = -2.1
        elif (Mstar >= 10.0) & (Mstar < 10.5): 
            alpha = -2.2
        elif (Mstar >= 10.5) & (Mstar < 11.0): 
            alpha = -2.0
        elif (Mstar >= 11.0) & (Mstar <= 11.5): 
            alpha = -1.3
        else: 
            raise NameError('Mstar is out of range')

        output = qf_z0 * ( 1.0 + z_in )**alpha 
        if output < 0.0: 
            output = 0.0
        elif output > 1.0: 
            output = 1.0 

        return output 
    
    elif lit == 'wetzelsmooth': 

        #qf_z0 = -6.04 + 0.63*Mstar
        qf_z0 = -6.04 + 0.64*Mstar
        alpha = -1.75

        output = qf_z0 * ( 1.0 + z_in )**alpha 
        if min(output) < 0.0: 
            output[np.where(output < 0.0)] = 0.0
        if max(output) > 1.0: 
            output[np.where(output > 1.0)] = 1.0

        return output 

    else: 
        raise NameError('Not yet coded') 

def get_fq_nsnap(Mstar, nsnap, **kwargs): 
    ''' Calculate the quiescent fraction as a funcnction of 
    stellar mass and redshift. Different methods of calculating 
    quiescent fraction is implemented. 

    f_Q ( M_star, z) 

    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    Mstar : stellar mass
    nsnap : redshift 
    lit : 'cosmosinterp', 

    '''
    z_snap = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])

    return get_fq(Mstar, z_snap[nsnap], **kwargs) 

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
    #return -0.75 + 0.76*(zin-0.05) + 0.5*(mstar-10.5)
    return -0.75 + 0.76*(zin-0.04) + 0.5*(mstar-10.5)

def cq_fq(cenque): 
    """
    Calculate the quiescent fraction of CenQue object. Note
    that the CenQue object has the following attributes: 
    mass, SFR, zsnap (all the necessary properties for f_Q 
    classification). 
    """

    # check that all necessary attributes are present
    if cenque.zsnap is None: 
        raise ValueError
    if cenque.mass is None: 
        raise ValueError
    if cenque.sfr is None: 
        raise ValueError

    # Star-forming or Quiescent    
    sf_q = sfq_classify(cenque.mass, cenque.sfr, cenque.zsnap, Mrcut=18)

    mass_bins = cenque.mass_bins
    
    masses, f_q = [], [] 
    for i_bin in xrange(mass_bins.nbins): 

        in_massbin = np.where( 
                (cenque.mass > mass_bins.mass_low[i_bin]) &
                (cenque.mass <= mass_bins.mass_high[i_bin])
                )

        n_massbin = len(in_massbin[0])

        if n_massbin == 0: 
            continue

        q = len(np.where(sf_q[in_massbin] == 'quiescent')[0])
        
        masses.append( mass_bins.mass_mid[i_bin] ) 
        f_q.append( np.float(q)/np.float(n_massbin) )

    return np.array(masses), np.array(f_q)

"""
def get_fq_alpha(Mstar, z_in, alpha): 
    ''' Quiescent fraction evolved from z = 0.88 by (1+z)^alpha where alpha is a free parameter
    '''

    fq_file = ''.join([ 'dat/wetzel_tree/', 
        'qf_z0.88cen.dat' ]) 
           
    # read in mass and quiescent fraction at z = 0.88 
    mass, fq = np.loadtxt(fq_file, unpack=True, usecols=[0,1])  
     
    fq_mstar_highz = np.interp(Mstar, mass, fq)  # interpolate to get fq(Mstar)
        
    output = fq_mstar_highz * np.abs(z_in + 0.16)**alpha 

    if output < 0.0: 
        output = 0.0
    elif output > 1.0: 
        output = 1.0 

    return output 

# 'quenching' fraction 
def get_fquenching(Mstar, z_in, **kwargs): 
    ''' Return the *quenching* fraction for given stellar mass and redshift 
    
    Parameters
    ----------
    Mstar : Stellar mass 
    z_in : Redshift

    Returns
    -------
    fquenching : quenching fraction 

    Notes
    -----
    * *Quenching* fraction is *not* quiescent fraction
    * Based on Wetzel et al. Quiescent Fraction parameterization 
    * Redshift evolution is the same as Wetzel     
    * As a first test use wetzel slope while varying yint 

    '''
    if 'slope' in kwargs.keys(): 
        slope = kwargs['slope']
    else: 
        slope = 0.63                # Wetzel slope

    if 'yint' in kwargs.keys(): 
        yint = kwargs['yint'] 
    else: 
        yint = -6.04

    qf_z0 = yint + slope*Mstar
    
    if Mstar < 9.5: 
        alpha = -2.3
    elif (Mstar >= 9.5) & (Mstar < 10.0): 
        alpha = -2.1
    elif (Mstar >= 10.0) & (Mstar < 10.5): 
        alpha = -2.2
    elif (Mstar >= 10.5) & (Mstar < 11.0): 
        alpha = -2.0
    elif (Mstar >= 11.0) & (Mstar < 11.5): 
        alpha = -1.3
    else: 
        raise NameError('Mstar is out of range')

    output = qf_z0 * ( 1.0 + z_in )**alpha 
    if output < 0.0: 
        output = 0.0
    elif output > 1.0: 
        output = 1.0 

    return output 
"""
