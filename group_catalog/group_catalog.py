'''

Package to analyze Jeremy's SDSS group catalog.
Works clunk-ily with CenQue Class object

Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import os 
import matplotlib.pyplot as plt
import h5py

#---- Local ----

class GroupCat: 
    def __init__(self): 
        pass

def read_group_catalog_galdata(Mrcut=18): 
    ''' wrapper for reading in group catalog galaxy data 
    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'

    h = 0.7
    
    file = ''.join(['dat/group_catalog/', 
        'clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.galdata_corr.fits']) 

    gal_data = mrdfits(file) 
    for column in gal_data.__dict__.keys(): 
        column_data = getattr(gal_data, column)
        if column == 'stellmass': 
            # stellmass is in units of Msol/h^2
            column_data = column_data / h**2
            setattr(gal_data, 'mass', np.log10(column_data))    # convert to log Mass
        elif column == 'ssfr': 
            column_data = column_data + np.log10(h**2) 
            setattr(gal_data, 'ssfr', column_data)    # convert to log Mass

        elif column == 'cz': 
            setattr(gal_data, 'z', column_data/299792.458)          # convert to z else: 
            pass

    setattr(gal_data, 'sfr', gal_data.mass + gal_data.ssfr)     # get sfr values

    return gal_data 

def read_group_catalog_prob(Mrcut=18): 
    ''' wrapper for reading in probability galaxy data 
    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    file = ''.join(['dat/group_catalog/', 
        'clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.prob.fits']) 

    prob_data = mrdfits(file)            # import probability file 

    return prob_data

def build_group_catalog(Mrcut=18, central=True):
    ''' only keep satellite or central galaxies based on p_sat
    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'

    # import gal_data and prob_data
    gal_data = read_group_catalog_galdata(Mrcut=Mrcut)
    prob_data = read_group_catalog_prob(Mrcut=Mrcut)

    # central or satellite 
    if central == True: 
        prob_index = prob_data.p_sat <= 0.5 
    else: 
        prob_index = prob_data.p_sat > 0.5
    
    central_catalog = GroupCat()       # save as CenQue class

    group_catalog_columns = ['mass', 'sfr', 'ssfr', 'z']
    for column in group_catalog_columns: 
        column_data = getattr(gal_data, column)[prob_index]
        setattr(central_catalog, column, column_data)

    output_file = ''.join(['dat/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.hdf5']) 
    central_catalog.writeout(columns=group_catalog_columns,
            input_file=output_file) 

def central_catalog(Mrcut=18, clobber=False): 
    ''' read in central group catalog into CenQue class
    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    catalog_file = ''.join(['dat/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.hdf5']) 

    if (os.path.isfile(catalog_file) == False) or (clobber == True): 
        build_group_catalog(Mrcut=Mrcut, central=True)
    
    f = h5py.File(catalog_file, 'r') 
    grp = f['cenque_data']
    mass = grp['mass'][:]
    sfr = grp['sfr'][:]
    ssfr = grp['ssfr'][:]
    z = grp['z'][:]
    #mass, sfr, ssfr = np.loadtxt(catalog_file, skiprows=2, delimiter=',', 
    #        unpack=True, usecols=[0,1,2]) 

    catalog = GroupCat() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    setattr(catalog, 'z', z) 
    
    f.close()
    return catalog 

def build_groupcat_sf(Mrcut=18): 
    ''' Build SF population for the SDSS group catalog for group catalog with specified Mrcut 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 

    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'

    # import centrals  
    centrals = central_catalog(Mrcut=Mrcut, clobber=True) 
    
    # classification motivated by Salim et al. 
    sf_gals = centrals.sfr > -1.30 + 0.65*(centrals.mass-10.0)

    centrals_sfms = GroupCat() 

    group_catalog_columns = ['mass', 'sfr', 'ssfr']
    for column in group_catalog_columns: 
        column_data = getattr(centrals, column)[sf_gals]
        setattr(centrals_sfms, column, column_data) 

    output_file = ''.join(['dat/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.starforming.hdf5']) 
    centrals_sfms.writeout(columns=group_catalog_columns,
            input_file=output_file) 

def sf_centrals(Mrcut=18, clobber=False): 
    ''' Read SDSS star-forming central group catalog into CenQue class

    Parameters
    ----------
    Mrcut : Absolute mangitude cut that specifies the group catalog 
    clobber : If True, re-construct the catalog. If False, just read catalog 

    Notes
    -----
    SF determined by a variation of the Salim et al. equation from Moustakas et al. 2013

    '''

    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    catalog_file = ''.join([
        'dat/group_catalog/', 
        'massSFR_clf_groups_M', 
        str(Mrcut), '_', str(masscut), 
        '_D360.central.starforming.hdf5'
        ]) 

    if not os.path.isfile(catalog_file) or clobber: 

        build_groupcat_sf(Mrcut=Mrcut, central=True)
    
    f = h5py.File(catalog_file, 'r') 
    grp = f['cenque_data']
    mass = grp['mass'][:]
    sfr = grp['sfr'][:]
    ssfr = grp['ssfr'][:]

    catalog = GroupCat() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    
    f.close()

    return catalog 

def build_groupcat_q(Mrcut=18): 
    ''' Build Q population for the SDSS group catalog for group catalog with specified Mrcut 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 

    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'

    # import centrals  
    centrals = cq_group.central_catalog(Mrcut=Mrcut, clobber=True) 
    
    # classification motivated by Salim et al. 
    sf_gals = centrals.sfr < -1.30 + 0.65*(centrals.mass-10.0)

    centrals_sfms = GroupCat() 

    group_catalog_columns = ['mass', 'sfr', 'ssfr']
    for column in group_catalog_columns: 
        column_data = getattr(centrals, column)[sf_gals]
        setattr(centrals_sfms, column, column_data) 

    output_file = ''.join(['dat/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.quiescent.hdf5']) 
    centrals_sfms.writeout(columns=group_catalog_columns,
            input_file=output_file) 

def q_centrals(Mrcut=18, clobber=False): 
    ''' Read SDSS quiescent central group catalog into CenQue class

    Parameters
    ----------
    Mrcut : Absolute mangitude cut that specifies the group catalog 
    clobber : If True, re-construct the catalog. If False, just read catalog 

    Notes
    -----
    Q determined by a variation of the Salim et al. equation from Moustakas et al. 2013

    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    catalog_file = ''.join(['dat/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.quiescent.hdf5']) 

    if (os.path.isfile(catalog_file) == False) or (clobber == True): 
        # (re)construct quiescent group catalog 
        build_groupcat_q(Mrcut=Mrcut)
    
    f = h5py.File(catalog_file, 'r') 
    grp = f['cenque_data']
    mass = grp['mass'][:]
    sfr = grp['sfr'][:]
    ssfr = grp['ssfr'][:]

    catalog = GroupCat() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    
    f.close()
    return catalog 

if __name__=='__main__': 
    central_ssfr = central_catalog(Mrcut=18, clobber=True) 
    central_ssfr = central_catalog(Mrcut=19, clobber=True) 
    central_ssfr = central_catalog(Mrcut=20, clobber=True) 
    central_catalog_match2isedfit(Mrcut=18, clobber=True) 
    central_catalog_match2isedfit(Mrcut=19, clobber=True) 
    central_catalog_match2isedfit(Mrcut=20, clobber=True) 
    #print double_gaussian_fit(central_ssfr) 
