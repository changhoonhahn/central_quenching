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
        if output.min() < 0.0: 
            output[np.where(output < 0.0)] = 0.0
        elif output.max() > 1.0: 
            output[np.where(output > 1.0)] = 1.0

        return output 

    else: 
        raise NameError('Not yet coded') 

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
