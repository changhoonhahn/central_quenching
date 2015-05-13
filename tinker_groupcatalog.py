import numpy as np
import pyfits as fits

class groupcatalog: 
    # import tinker et al. group catalog 
    def __init__(self): 
        pass    
    def readfile(self, Mrcut=18): 
        if Mrcut == 18: 
            masscut='9.4'
        elif Mrcut == 19: 
            masscut='9.8'
        elif Mrcut == 20: 
            masscut='10.2'
        catalogDir = 'dat/group_catalog/'
        filePre = 'clf_groups_M'
        fileExt = '.galdata_corr.fits'
        filename = ''.join([catalogDir, filePre, str(Mrcut), '_', str(masscut), '_D360', fileExt])
        self.file = filename 
        groupdata = mrdfits(filename) 
        columns = groupdata.__dict__.keys() 
        for column in columns:                  # import data from group catalog
            if column == 'stellmass':  
                setattr(self, 'mass', np.log10(getattr(groupdata, column))) 
            elif column == 'cz': 
                setattr(self, 'z', getattr(groupdata, column)/299792.458)
            else: 
                setattr(self, column, getattr(groupdata, column)) 
        #self.ssfr = np.log10(np.abs(self.sfr/(10.0**self.mass)))         # Originaly for calculating SSFR, but SSFR data is contained within group catalog corrdata
    def readProbFile(self): 
        '''
        read in .prob file for corresponding group catalog and extra p_sat data
        '''
        probfile = (self.file).rstrip('.galdata_corr.fits')
        probfile = ''.join([probfile, '.prob.fits'])
        self.probfile = probfile
        probdata = mrdfits(probfile)            # import probability file 
        probdata_column = probdata.__dict__.keys()
        for column in probdata_column:  
            setattr(self, column, getattr(probdata, column))
    def onlykeep(self, type='central'): 
        '''
        only keep satellite or central galaxies based on p_sat
        '''
        if type == 'central': 
            prob_index = self.p_sat <= 0.5 
        elif type == 'satellite': 
            prob_index = self.p_sat > 0.5
        self_columns = self.__dict__.keys() 
        for column in self_columns: 
            if not ((column == 'file') or (column  == 'probfile')): 
                setattr(self, column, getattr(self, column)[prob_index])

class treepm_evomock: 
    # import mocks with evolved galaxy properties 
    def __init__(self): 
        pass
    def readfile(self, nsnap, tau='constant', mass_bin=0.01): 
        '''
        Read mock file with evolved list of galaxy properties that are hardcoded in the code 
        '''
        evomock_file = ''.join(['dat/wetzel_tree/', 
            'subhalo_sham_centrals_snapshot', str(nsnap), '_ssfr_mbin', str(mass_bin), '_evolfrom_snapshot13_', tau, 'tau.fits'])
        evomock_data = mrdfits(evomock_file) 
        gal_props = ['sfq', 'sfr', 'mass', 'ssfr']
        for gal_prop in gal_props: 
            setattr(self, gal_prop, getattr(evomock_data, gal_prop)) 

class ssfr_hist: 
    # SSFR histogram object 
    def __init__(self, ssfr_array): 
        '''
        Given array of SSFR values, calculates the sSFR histogram 
        '''
        self.ssfr = ssfr_array 
        # hardcoded min and max SSFR range values 
        ssfr_min = -13.0   
        ssfr_max = -7.0 
        print np.min(ssfr_array), np.max(ssfr_array)
        ssfr_binsize = 0.1
        ssfr_nbin = (ssfr_max - ssfr_min)/ssfr_binsize
        
        hist_ssfr = np.histogram(ssfr_array, bins=ssfr_nbin, range=[ssfr_min, ssfr_max])
        self.hist_ssfr = hist_ssfr[0]
        self.ssfr_low = (hist_ssfr[1])[0:-1]
        self.ssfr_high = (hist_ssfr[1])[1:]  
        self.ssfr_mid =  np.array([0.5*(self.ssfr_low[i]+self.ssfr_high[i]) for i in range(len(self.ssfr_low))])

def ssfr_chi2(ssfr1, ssfr2): 
    '''
    Given two SSFR histogram objects calculate the chi-squared value
    '''
    chi2 = np.sum((ssfr1.hist_ssfr - ssfr2.hist_ssfr)**2) 
    return chi2

def group_evomock_bestfit(): 
    '''
    Find the best fit 
    '''

def writeCentralSatellite(Mrcut=18): 
    '''
    read group catalog and write central and satellite galaxies 
    '''
    for galtype in ['central', 'satellite']: 
        data = groupcatalog() 
        data.readfile(Mrcut=Mrcut)
        fname = data.file
        data.readProbFile()
        data.onlykeep(type=galtype)
        filedir = fname.rsplit("/",1)[0]
        filename = (fname.rsplit("/",1)[1]).rsplit(".",1)[0]
        fname = ''.join([filedir, '/massSFR_', filename, '.', galtype, '.fits'])
        mwrfits(data, fname) 

if __name__=="__main__": 
    #writeCentralSatellite(Mrcut=18)
    for Mr in [18, 19, 20]: 
        writeCentralSatellite(Mrcut=Mr) 

        #tinkergroup = groupcatalog() 
        #tinkergroup.readfile(Mrcut=Mr) 
        #tinkergroup.readProbFile()
        #tinkergroup.onlykeep(type='central')
        #tinkergroup_ssfr = ssfr_hist(tinkergroup.ssfr) 

        #evomock = treepm_evomock() 
        #evomock.readfile(1)
        #evomock_ssfr = ssfr_hist(evomock.ssfr) 
        #print ssfr_chi2(tinkergroup_ssfr, evomock_ssfr) 
