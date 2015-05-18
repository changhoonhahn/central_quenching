'''

Utility functions for CenQue

Author(s): ChangHoon Hahn

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import os 
import mpfit

# ---- Local -----
import cenque as cq 
import cenque_groupcat as cq_group
import sf_mainseq as sfms 


# quiescent fraction ------------------------------------------------------------
class fq: 
    ''' quiescent fraction 
    '''
    def __init__(self): 
        self.z = None       # redshift
        self.mass = None
        self.fq = None      # quiescent fraction 

def get_fq(Mstar, z_in, lit='cosmosinterp'): 
    ''' Return quiescent fraction from selected literature

    Parameters
    ----------
    Mstar : stellar mass
    z_in : redshift 
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

# SSFR distribution ----------------------------------------------------------------

class CQssfr: 
    ''' SSFR distribution for CenQue class
    '''
    def __init__(self): 
        self.ssfr = None 
        self.bin_low = None     # SSFR bin low 
        self.bin_high = None    # SSFR bin high 
        self.bin_mid = None     # SSFR bin mid 
        self.ssfr_hist = None   # histogram 

    def histogram(self): 
        ''' Calculate SSFR histogram 
        '''
        if self.ssfr == None: 
            raise NameError('Not SSFR specified') 

    def readin(self, **kwargs): 
        ''' Read in ssfr histogram 
        '''

    def writeout(self, **kwargs): 
        ''' Writes out ssfr histogram 
        '''

# mass bins ------------------------------------------------------------------------ 
class mass_bin:         
    ''' Mass bin class  
    '''
    def __init__(self): 
        self.mass_mid = None 
        self.mass_low = None 
        self.mass_high = None 
        self.nbins = None

    def build_delta_mass(delta_mass, mass_min=9.5, mass_max=11.5): 
        ''' build mass bins given min, max mass and delta 
        '''
        pass

def simple_mass_bin(): 
    ''' Simple mass bins 

    output 
    ------
    mass_bin class that contains mass bin information
    '''
    simple_mass_binsize = 0.2
    simple_mass_bin = mass_bin()
    simple_mass_bin.mass_low = [ 9.0 + np.float(i)*simple_mass_binsize for i in range(13) ]
    simple_mass_bin.mass_high = [ simple_mass_bin.mass_low[i] + simple_mass_binsize
            for i in range(len(simple_mass_bin.mass_low)) ]
    simple_mass_bin.mass_mid = [
            0.5 * (simple_mass_bin.mass_low[i] + simple_mass_bin.mass_high[i]) 
            for i in range(len(simple_mass_bin.mass_low))]

    simple_mass_bin.mass_wid = [simple_mass_bin.mass_high[i] - simple_mass_bin.mass_low[i] 
            for i in range(len(simple_mass_bin.mass_high))]

    simple_mass_bin.nbins = len(simple_mass_bin.mass_low)

    return simple_mass_bin 

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

def get_sfr_mstar_z_bestfit(mstar, z_in, Mrcut=18, clobber=False):
    ''' Return SFR of SF main sequence as a function of mass and redshift

    Parameters
    ----------
    mstar : Stellar mass of galaxy
    z_in : Redshift 
    Mrcut : Absolute magnitude cut that specified the group catalog 

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
    fid_mass = 10.5

    # Best-fit slope and y-int of SF SDSS Group Catalog 
    groupcat_fit_param = sfms.get_bestfit_groupcat_sfms(Mrcut=Mrcut, clobber=clobber)

    # Best-fit slope and y-int of SF EnvCount
    envcount_fit_param = sfms.get_bestfit_envcount_sfms()
    zmids, slopes, yints = sfms.get_sfmsfit_sfr(
            (envcount_fit_param[0]).item(), (envcount_fit_param[1]).item(), 
            clobber=clobber)
        
    d_yints = np.interp(z_in, zmids, yints) - yints[0] 
    SFR_amp = groupcat_fit_param[1] + d_yints

    avg_SFR = groupcat_fit_param[0] * (mstar - fid_mass) + SFR_amp
        
    return [avg_SFR, 0.3] 

def sfq_classify(mstar, sfr, z_in, Mrcut=18, clobber=False):
    ''' Return SF or Q classification given M* and SFR 

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
    ssfr_cut = -11.7 + 0.98*z_in - 0.25*(mstar-10.25)
    sfr_cut = ssfr_cut + mstar 

    sfq = np.empty(len(mstar), dtype=(str,16))
    sf_index = sfr > sfr_cut 
    sfq[sf_index] = 'star-forming'
    q_index = sfr <= sfr_cut
    sfq[q_index] = 'quiescent'
    return sfq 

def line(x, p): 
    return p[0]*x + p[1]

def mpfit_line(p, fjac=None, x=None, y=None): 
    #model = line_fixedslope(x, p) 
    model = line(x, p) 
    status = 0 
    return([status, (y-model)]) 

'''
def sdss_sf_ms_fit(): 

    p0 = [0.02]
    fa = {'x': np.array([10.25, 10.75, 11.25]), 'y': np.array([-0.0777, 0.248, 0.483])}
    bestfit_pars = mpfit.mpfit(mpfit_line, p0, functkw=fa, nprint=0)
    print bestfit_pars.params[0]
    print line_fixedslope(10.25, bestfit_pars.params) 
    print line_fixedslope(10.75, bestfit_pars.params) 
    print line_fixedslope(11.0, bestfit_pars.params)
'''

def get_quenching_efold(mstar, type='constant', param=None): 
    ''' get quenching efold based on stellar mass of galaxy 
    '''

    if type == 'constant':      # constant tau 

        n_arr = len(mstar) 
        tau = np.array([0.5 for i in range(n_arr)]) 

    elif type == 'linear':      # lienar tau(mass) 

        tau = -(0.9 / 1.5) * ( mstar - 9.5) + 1.0
        if np.min(tau) < 0.1:
            tau[ tau < 0.1 ] = 0.1
         
    elif type == 'instant':     # instant quenching 

        n_arr = len(mstar) 
        tau = np.array([0.001 for i in range(n_arr)]) 

    elif type == 'discrete': 
        # param will give 4 discrete tau at the center of mass bins 
        masses = np.array([9.75, 10.25, 10.75, 11.25]) 

        if param is None: 
            raise ValueError('asdfasdfa') 

        tau = np.interp(mstar, masses, param) 
        tau[ tau < 0.05 ] = 0.05

    elif type == 'linefit': 
        # param will give slope and yint of pivoted tau line 
        
        tau = param[0] * (mstar - 10.5) + param[1]
        tau[ tau < 0.05 ] = 0.05

    else: 
        raise NotImplementedError('asdf')

    return tau 

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
    fit_line_param = sfms.get_bestfit_qgroupcat_ssfr(Mrcut=Mrcut) 
    fit_slope = fit_line_param[0].item() 
    fit_yint = fit_line_param[1].item() 

    q_ssfr = np.array([ fit_slope * (mass - 10.5) + fit_yint - 0.1 for mass in masses ])  

    #q_ssfr = np.array([ (-0.7 * mass) - 4.625 for mass in masses ])  
    return q_ssfr 

# fits file treatment -----------------------------------------------------------------------
class FitsTable:
    def __init__(self): 
        pass
    def columns(self): 
        return self.__dict__.keys()

def mrdfits(filename): 
    output = FitsTable()
    fitsdata = fits.open(filename)[1].data
    for name in fitsdata.names: 
        setattr(output, name.lower(), fitsdata.field(name))
    return output 

def mwrfits(fitstable, filename, columns=[], clobber=1): 
    typedict = dict([(d, type(np.zeros(1,d).tolist()[0])) for d in (np.float32,np.float64,np.uint32, np.int16)])
   
    if len(columns) == 0: 
        # if columns aren't specified then all type of columns will be inputed
        fitscolumn = fitstable.__dict__.keys() 
    else: 
        if all(isinstance(x,str) for x in columns): 
            fitscolumn = columns
        else: 
            raise TypeError('columns need to be strings') 

    table_arrays = [] 
    for column in fitscolumn: 
        column_array = getattr(fitstable, column)

        if type(column_array[0]) in typedict.keys(): 
            column_type = typedict[type(column_array[0])]
        else: 
            column_type = type(column_array[0])
        
        if column_type == str: 
            max_str_len = np.max(np.array([len(column_array[i]) for i in range(len(column_array))]))
            format_str = str(max_str_len)+'A'
            column_array = np.array([column_array])
        elif column_type == int: 
            format_str = 'K'  
        elif column_type == float: 
            format_str = 'E'
        else:                               # in the case that it's an array of arrays
            array_len = len(column_array[0])    
            if type((column_array[0])[0]) == int: 
                format_str = str(array_len)+'K'  
            elif type((column_array[0])[0]) == float: 
                format_str = str(array_len)+'E'

        column_array = fits.Column(name=column, format=format_str, array=column_array) 
        table_arrays.append(column_array)

    table_arrays_combined = fits.ColDefs(table_arrays)
    real_fitstable = fits.new_table(table_arrays_combined)
    real_fitstable.writeto(filename, clobber=clobber) 

# CenQue file treatment ----------------------------------------------------------------------- 
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
