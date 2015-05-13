'''

Central_Quenchign project
codes for SDSS group catalog


Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import os 
import matplotlib.pyplot as plt
import h5py

#---- Local ----
import cenque_utility as util
import cenque as cq 
import mpfit 
import pyspherematch as pysph

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
    
    file = ''.join(['/dat/group_catalog/', 
        'clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.galdata_corr.fits']) 

    gal_data = util.mrdfits(file) 
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
    
    file = ''.join(['/dat/group_catalog/', 
        'clf_groups_M', str(Mrcut), '_', str(masscut), '_D360.prob.fits']) 

    prob_data = util.mrdfits(file)            # import probability file 

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

    catalog = cq.CenQue() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    setattr(catalog, 'z', z) 
    
    f.close()
    return catalog 

def double_gaussian(x, p): 
    '''
    double gaussian function
    '''
    # p includes the parameters
    # p[0] = fred 
    # p[1] = mean SSFR for quiescent
    # p[2] = mean SSFR for star-forming  

    fred = p[0]
    ssfr_q = p[1]
    ssfr_sf = p[2]
    ssfr_q_sig = p[3]
    ssfr_sf_sig = p[4]

    phi = fred/(np.sqrt(2*np.pi)*ssfr_q_sig) * np.exp(-1.0*( x - ssfr_q )**2 / (2*ssfr_q_sig**2)) + \
            (1 - fred)/(np.sqrt(2*np.pi)*ssfr_sf_sig) * np.exp(-1.0*( x - ssfr_sf )**2 / (2*ssfr_sf_sig**2))

    return phi

def mpfit_doubleG(p, fjac=None, x=None, y=None): 
    model = double_gaussian(x, p) 
    status = 0 
    return([status, (y-model)]) 

def double_gaussian_fit(cen_cat): 
    ''' Fits a double guassian fit to the SSFR distribution of the central catalog 
    '''
    
    # hardcoded mass bins
    panel_mass_bins = [
            [10.0, 10.5], [10.5, 11.0], [11.0, 11.5]
            ]
    
    output = []
    for i_mass, panel_mass in enumerate(panel_mass_bins):       # loop through each panel 

        mass_limit = (cen_cat.mass >= panel_mass[0]) & (cen_cat.mass < panel_mass[1])  
        
        # calculate SSFR histogram for mass bin 
        ssfr_hist, ssfr_bin_edges = np.histogram(cen_cat.ssfr[mass_limit], 
                range=[-13.0, -7], bins=40, normed=True)

        ssfr_bin_low = ssfr_bin_edges[:-1]
        ssfr_bin_high = ssfr_bin_edges[1:]
        ssfr_bin_mid = [ 0.5*(ssfr_bin_low[i] + ssfr_bin_high[i]) 
                for i in range(len(ssfr_bin_low)) ] 
    
        # guesses for parameters
        fred_guess = util.get_fq(np.mean(panel_mass), 0.05, lit='wetzel')   # wetzel's Fq

        ssfr_q_guess = -12.0            # hardcoded -12.0

        [sf_avg_sfr, sf_sig_sfr] = util.get_sfr_mstar_z(np.mean(panel_mass), 0.05) 

        ssfr_sf_guess = sf_avg_sfr - np.mean(panel_mass)             # from SF MS 

        p0 = [fred_guess, ssfr_q_guess, ssfr_sf_guess, 0.2, 0.2]            # initial guess

        print 'Guess ', p0
    
        fa = {'x': ssfr_bin_mid, 'y': ssfr_hist}

        bestfit_pars = mpfit.mpfit(mpfit_doubleG, p0, functkw=fa, nprint=0)
        
        print 'Bestfit ', bestfit_pars.params

        output.append([ panel_mass[0], panel_mass[1], bestfit_pars.params]) 

    return output 

def match2_primus_mfdata(Mrcut=18): 
    ''' Experimental data to match Jeremy's group catalog to John's iSEDfit SDSS Dr7 catalog
    '''
    # import group catalog 

    gal_data = read_group_catalog_galdata(Mrcut=Mrcut)
    prob_data = read_group_catalog_prob(Mrcut=Mrcut)
    
    gal_ra = gal_data.ra * 57.2957795 
    gal_dec = gal_data.dec * 57.2957795 

    # import sdss mf data 
    mfdata = util.mrdfits('dat/group_catalog/mfdata_all_supergrid01_sdss.fits.gz') 
    mf_ra = mfdata.ra
    mf_dec = mfdata.dec

    mf_index, gal_index, d = pysph.spherematch(mf_ra, mf_dec, gal_ra, gal_dec, tol=0.0001)
    print len(mf_index), ' spherematches'
    
    fig = plt.figure(figsize=[8,8]) 
    sub = fig.add_subplot(111)

    sub.scatter(gal_data.mass[gal_index], mfdata.mass[mf_index], c='k') 
    
    sub.set_xlabel('k-correct Mass') 
    sub.set_ylabel('iSEDfit Mass') 
    sub.set_ylim([9.0, 12.0])
    sub.set_xlim([9.0, 12.0])
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'primus_mfdata_', str(Mrcut), '_galdatacorr_mass_comp.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 

def build_group_catalog_match2isedfit(Mrcut=18, central=True):
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

    gal_ra = gal_data.ra * 57.2957795 
    gal_dec = gal_data.dec * 57.2957795 
    
    # import sdss mf data 
    mfdata = util.mrdfits('dat/group_catalog/mfdata_all_supergrid01_sdss.fits.gz') 

    mf_ra = mfdata.ra
    mf_dec = mfdata.dec

    mf_index, gal_index, d = pysph.spherematch(mf_ra, mf_dec, gal_ra, gal_dec, tol=0.0001)

    index = np.arange(len(gal_ra))
    match_index = index[gal_index]
    # central or satellite 
    if central == True: 
        prob_index = (prob_data.p_sat[gal_index] <= 0.5)
    else: 
        prob_index = (prob_data.p_sat[gal_index] > 0.5)
    
    central_catalog = cq.CenQue()       # save as CenQue class

    group_catalog_columns = ['mass', 'sfr', 'ssfr', 'kcorrect_mass', 'mpajhu_sfr']
    for column in group_catalog_columns: 
        #column_data = getattr(gal_data, column)[prob_index]
        #print len(column_data)
        if column == 'mass': 
            setattr(central_catalog, column, mfdata.mass[mf_index[prob_index]])
        elif column == 'sfr': 
            setattr(central_catalog, column, mfdata.sfr[mf_index[prob_index]])
        elif column == 'ssfr': 
            setattr(central_catalog, column, mfdata.sfr[mf_index[prob_index]] - mfdata.mass[mf_index[prob_index]])

        elif column == 'kcorrect_mass': 
            setattr(central_catalog, column, gal_data.mass[gal_index[prob_index]])

        elif column == 'mpajhu_sfr': 
            setattr(central_catalog, column, gal_data.sfr[gal_index[prob_index]])

        else: 
            raise NameError('asdfasdfasdf') 

    output_file = ''.join(['dat/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360_match2isedfit.central.hdf5']) 
    central_catalog.writeout(columns=group_catalog_columns,
            input_file=output_file) 

def central_catalog_match2isedfit(Mrcut=18, clobber=False): 
    ''' read in central group catalog into CenQue class
    '''
    if Mrcut == 18: 
        masscut='9.4'
    elif Mrcut == 19: 
        masscut='9.8'
    elif Mrcut == 20: 
        masscut='10.2'
    
    catalog_file = ''.join(['dat/group_catalog/', 
        'massSFR_clf_groups_M', str(Mrcut), '_', str(masscut), '_D360_match2isedfit.central.dat']) 

    if (os.path.isfile(catalog_file) == False) or (clobber == True): 
        build_group_catalog_match2isedfit(Mrcut=Mrcut, central=True)

    f = h5py.File(catalog_file) 
    mass = f['cenque_data/mass'][:]
    sfr = f['cenque_data/sfr'][:]
    ssfr = f['cenque_data/ssfr'][:]
    kcorrect_mass = f['cenque_data/kcorrect_mass'][:]
    mpajhu_sfr = f['cenque_data/mpajhu_sfr'][:]

    catalog = cq.CenQue() 
    setattr(catalog, 'mass', mass) 
    setattr(catalog, 'sfr', sfr) 
    setattr(catalog, 'ssfr', ssfr) 
    setattr(catalog, 'kcorrect_mass', kcorrect_mass) 
    setattr(catalog, 'mpajhu_sfr', mpajhu_sfr) 
    return catalog 

if __name__=='__main__': 
    central_ssfr = central_catalog(Mrcut=18, clobber=True) 
    central_ssfr = central_catalog(Mrcut=19, clobber=True) 
    central_ssfr = central_catalog(Mrcut=20, clobber=True) 
    central_catalog_match2isedfit(Mrcut=18, clobber=True) 
    central_catalog_match2isedfit(Mrcut=19, clobber=True) 
    central_catalog_match2isedfit(Mrcut=20, clobber=True) 
    #print double_gaussian_fit(central_ssfr) 
