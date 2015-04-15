'''

Central_Quenchign project
codes for SDSS group catalog


Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import os 

#---- Local ----
import cenque_utility as util
import cenque as cq 

def read_group_catalog_galdata(Mrcut=18): 
    ''' wrapper for reading in group catalog galaxy data 
    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    file = ''.join(['/data1/hahn/group_catalog/', 
        'clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.galdata_corr.fits']) 

    gal_data = mrdfits(file) 
    for column in gal_data.__dict__.keys(): 
        column_data = getattr(gal_data, column)
        if column == 'stellmass': 
            setattr(gal_data, 'mass', np.log10(column_data))    # convert to log Mass
        elif column == 'cz': 
            setattr(gal_data, 'z', column_data/299792.458)          # convert to z 
        else: 
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
    
    file = ''.join(['/data1/hahn/group_catalog/', 
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
    
    central_catalog = cq.CenQue()       # save as CenQue class

    group_catalog_columns = ['mass', 'sfr', 'ssfr']
    for column in group_catalog_columns: 
        column_data = getattr(gal_data, column)[prob_index]
        setattr(central_catalog, column, column_data)

    output_file = ''.join(['/data1/hahn/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.dat']) 
    central_catalog.writeout(columns=group_catalog_columns,
            input_file=output_file) 

def central_catalog(Mrcut=18): 
    ''' read in central group catalog into CenQue class
    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    catalog_file = ''.join(['/data1/hahn/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.central.dat']) 

    if os.path.isfile(catalog_file) == False: 
        build_group_catalog(Mrcut=Mrcut, central=True)

    mass, sfr, ssfr = np.loadtxt(catalog_file, skiprows=2, delimiter=',', 
            unpack=True, usecols=[0,1,2]) 

    catalog = cq.CenQue() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    
    return catalog 
