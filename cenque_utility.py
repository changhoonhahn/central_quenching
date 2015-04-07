'''

Utility functions for CenQue

Author(s): ChangHoon Hahn

'''

import numpy as np
import pyfits as fits
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
        
        return np.interp(z_in, zbins,  fq_z) 
    else: 
        raise NameError('Not yet coded') 
        

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
    simple_mass_bin = mass_bin()
    simple_mass_bin.mass_low = [ 9.0, 9.5, 10.0, 10.5, 11.0, 11.5]
    simple_mass_bin.mass_high = [ 9.5, 10.0, 10.5, 11.0, 11.5, 12.0]
    simple_mass_bin.mass_mid = [
            0.5 * (simple_mass_bin.mass_low[i] + simple_mass_bin.mass_high[i]) 
            for i in range(len(simple_mass_bin.mass_low))]

    simple_mass_bin.mass_wid = [simple_mass_bin.mass_high[i] - simple_mass_bin.mass_low[i] 
            for i in range(len(simple_mass_bin.mass_high))]

    simple_mass_bin.nbins = len(simple_mass_bin.mass_low)

    return simple_mass_bin 

# SFR-Mstar ---------------------------------------------------------------------------
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
