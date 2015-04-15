'''

Utility functions for CenQue

Author(s): ChangHoon Hahn

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import os 

# quiescent fraction ------------------------------------------------------------
class fq: 
    ''' quiescent fraction 
    '''
    def __init__(self): 
        self.z = None       # redshift
        self.mass = None
        self.fq = None      # quiescent fraction 

def get_fq(Mstar, z_in, lit='cosmosinterp'): 
    ''' get quiescent fraction from literature
    '''
    if lit == 'cosmosinterp': 
        zbins = [0.36, 0.66, 0.88] 

        fq_z = [] 
        for zbin in zbins: 
            fq_file = ''.join([ '/data1/hahn/wetzel_tree/', 
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
    elif lit == 'wetzel': 
        qf_z0 = -6.04 + 0.63*Mstar

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

        return output 
    else: 
        raise NameError('Not yet coded') 

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
    '''
    simple_mass_binsize = 0.1
    simple_mass_bin = mass_bin()
    simple_mass_bin.mass_low = [ 9.0 + np.float(i)*simple_mass_binsize for i in range(25) ]
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
        sf_ms_file = ''.join(['/data1/hahn/wetzel_tree/envcount/', 
            'sfr_mstar_fit_param_fixedslope.fits']) 
        sf_ms = mrdfits(sf_ms_file) 
        
        fid_mass = 10.5         # fid mass (part of the fititng) 
            
        SFR_amp = np.interp(z_in, sf_ms.z, sf_ms.yint)  # interpolate SFR amplitude by z input

        # closest redshift index
        closest_i = min(range(len(sf_ms.z)), key=lambda i: abs(sf_ms.z[i] - z_in))
    
        # assuming slope doesn't change calculate average SFR
        avg_sfr = (sf_ms.slope)[closest_i] * (mstar - fid_mass) + SFR_amp

        return [avg_sfr, 0.3]       # 0.3 dex scatter hard-coded

    elif lit == 'primus':       # using PRIMUS central galaxy data 
        
        if machine == 'harmattan': 
            data_dir = '/data1/hahn/wetzel_tree/envcount/'
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

def get_quenching_efold(mstar, type='constant'): 
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
    else: 
        raise NotImplementedError('asdf')

    return tau 

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

    MORE FILE TYPES WILL BE SPECIFIED
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
            
                file_type_str = ''.join(['_', kwargs['fq'], 'fq', file_type_str])
                
            elif min( (kwargs['file_type']).find('evo'), (kwargs['file_type']).find('from') ) > -1: 
                # evolved from nsnap
                original_nsnap = int(((kwargs['file_type']).split('from'))[-1]) 

                file_type_str = '_evol_from'+str(original_nsnap) 
                
                file_type_str = ''.join(['_', kwargs['tau'], 'tau_', 
                    kwargs['fq'], 'fq', file_type_str])
            else: 
                raise NameError("File not specified") 

        cenque_filename = ''.join(['/data1/hahn/central_quenching/', 
            'cenque_centrals_snapshot', str(kwargs['nsnap']), file_type_str, 
            '.dat']) 

    return cenque_filename
