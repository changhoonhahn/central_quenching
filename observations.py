'''


Module for observed galaxy catalogs with galactic properties. 
Jeremy Tinker's Group Catalog along with John Moustakas's 
iSEDfit SDSS + PRIMUS catalogs with my environment values. 


'''
import os 
import h5py
import time 
import numpy as np
from pydl.pydlutils.spheregroup import spherematch

from util import mpfit 
from util import util 

from gal_prop import Fq
from gal_prop import Ssfr

from defutility.fitstables import mrdfits 


class GroupCat(object): 
    def __init__(self, Mrcut=18, position='central', **kwargs): 
        ''' Class object that describes SDSS Group Catalog galaxy data
        '''
        self.kwargs = kwargs.copy()
        # central/satellite/all
        self.position = position
        
        self.Mrcut = Mrcut  # Mr cut
        # mass cuts based on Mr cuts
        if self.Mrcut == 18: 
            self.masscut='9.4'
        elif Mrcut == 19: 
            self.masscut='9.8'
        elif Mrcut == 20: 
            self.masscut='10.2'
        
        self.columns = ['ra', 'dec', 'mass', 'sfr', 'ssfr', 'z']

    def Ssfr(self): 
        ''' Calculate SSFR distribution of group catalog 
        '''
        try: 
            self.mass
        except AttributeError:
            self.Read()

        ssfr_obj = Ssfr()
        return ssfr_obj.Calculate(self.mass, self.ssfr)

    def File(self): 
        ''' GroupCat File 
        '''
        output_file = ''.join([
            'dat/group_catalog/', 
            'GroupCat.', 
            'Mr', str(self.Mrcut), 
            '.Mass', str(self.masscut), 
            '.D360.', 
            self.position, 
            '.hdf5']) 
        return output_file 

    def Write(self): 
        ''' Write compiled GroupCat to hdf5 file. 
        '''
        try: 
            self.mass
        except AttributeError:
            raise ValueError
        
        f = h5py.File(self.File(), 'w')
        grp = f.create_group('data')
        
        for column in self.columns: 
            column_attr = getattr(self, column)
            grp.create_dataset(column, data=column_attr)

        f.close()
        return None

    def Read(self): 
        ''' Read GroupCat from hdf5 file 
        '''

        f = h5py.File(self.File(), 'r') 
        grp = f['data']

        for column in grp.keys(): 
            setattr(self, column, grp[column][:])
        
        f.close()
        return None 
    
    def CompileGroupCat(self): 
        ''' Build Central/all Group Catalog 
        '''
        # import gal_data and prob_data
        gal_data = self._ReadGroupCat_GalData()
        prob_data = self._ReadGroupCat_Prob()

        # central or satellite 
        if self.position == 'central': 
            prob_index = np.where(prob_data.p_sat <= 0.5)
        elif self.position == 'satellite': 
            prob_index = np.where(prob_data.p_sat > 0.5)
        elif self.position == 'all': 
            prob_index = np.range(len(prob_data.p_sat))
    
        for column in self.columns: 
            column_data = getattr(gal_data, column)[prob_index]
            if column in ['ra', 'dec']: 
                setattr(self, column, column_data * 57.2957795)
            else: 
                setattr(self, column, column_data)
        return None

    def _ReadGroupCat_GalData(self, Mrcut=18): 
        '''
        '''
        h = 0.7
        
        galdata_file = ''.join([
            'dat/group_catalog/', 
            'clf_groups_M', str(self.Mrcut), '_', str(self.masscut), '_D360.', 
            'galdata_corr.fits']) 
        gal_data = mrdfits(galdata_file) 
        #print gal_data.__dict__.keys()
        for column in gal_data.__dict__.keys(): 
            column_data = getattr(gal_data, column)
            if column == 'stellmass': 
                # stellmass is in units of Msol/h^2
                column_data = column_data / h**2
                # convert to log Mass
                setattr(gal_data, 'mass', np.log10(column_data))
            elif column == 'ssfr': 
                column_data = column_data + np.log10(h**2)
                # convert to log Mass 
                setattr(gal_data, 'ssfr', column_data)    

            elif column == 'cz': 
                # convert to z else: 
                setattr(gal_data, 'z', column_data/299792.458)
            else: 
                pass

        setattr(gal_data, 'sfr', gal_data.mass + gal_data.ssfr)     # get sfr values
        
        return gal_data

    def _ReadGroupCat_Prob(self, Mrcut=18): 
        ''' wrapper for reading in probability galaxy data 
        '''
        file = ''.join([
            'dat/group_catalog/', 
            'clf_groups_M', str(self.Mrcut), '_', str(self.masscut), '_D360.', 
            'prob.fits']) 

        prob_data = mrdfits(file)            # import probability file 
        return prob_data

    def _iSEDfitMatch(self): 
        ''' Match the GroupCat galaxies with iSEDfit galaxy properties from 
        John Moustakas's MFData objects. The matching is done using PyDL's 
        spherematch
        '''
        # import SDSS MFdata catalog
        mfdata = mrdfits('dat/group_catalog/mfdata_all_supergrid01_sdss.fits.gz') 
        spherematch_time = time.time()
        match1, match2, d_match = spherematch(
                self.ra, self.dec, mfdata.ra, mfdata.dec, 0.001)
        print 'spherematch took ', time.time() - spherematch_time

        iSEDfit_mass = np.repeat(-999., len(self.ra))
        iSEDfit_SFR = np.repeat(-999., len(self.ra))
        iSEDfit_SSFR = np.repeat(-999., len(self.ra))
        
        if np.max(self.z[match1] - mfdata.z[match2]) > 0.1: 
            raise ValueError
        #wrong = np.argmax(self.z[match1] - mfdata.z[match2])
        #print self.ra[match1[wrong]], self.dec[match1[wrong]], self.z[match1[wrong]]
        #print mfdata.ra[match2[wrong]], mfdata.dec[match2[wrong]], mfdata.z[match2[wrong]]
        iSEDfit_mass[match1] = mfdata.mass[match2]
        iSEDfit_SFR[match1] = mfdata.sfr[match2]
        iSEDfit_SSFR[match1] = iSEDfit_SFR[match1]-iSEDfit_mass[match1]
        
        setattr(self, 'iSEDfit_mass', iSEDfit_mass) 
        setattr(self, 'iSEDfit_sfr', iSEDfit_SFR) 
        setattr(self, 'iSEDfit_ssfr', iSEDfit_SSFR) 
        return None


class PrimusSDSS(object): 
    def __init__(self, redshift, environment='no', **kwargs):
        ''' Class that describes SDSS and PRIMUS galaxy data used in the 
        galaxy environment paper. 
        '''
        self.kwargs = kwargs
        self.redshift = redshift
        self.environment = environment 

    def Read(self): 
        '''
        '''
        file_dir = 'dat/wetzel_tree/envcount/'
        sdss_file = ''.join([file_dir, 
                'envcount_cylr2.5h35_thresh75_sdss_active_z0.05_0.12_primuszerr.fits']) 
        primus_file = ''.join([file_dir, 
                'envcount_cylr2.5h35_thresh75_active_z0.2_1.0_lit.fits']) 
        # redshift determines SDSS or PRIMUS sample
        if self.redshift < 0.2: 
            galdata = mrdfits(sdss_file) 
        else: 
            galdata = mrdfits(primus_file)
        # environment cuts 
        if self.environment == 'no': 
            envcuts = (galdata.envcount == 0.)
        else: 
            raise NotImplementedError
        # redshift cuts  
        zlow = [0., 0.2, 0.4, 0.6, 0.8]
        zhigh = [0.2, 0.4, 0.6, 0.8, 1.0]
        i_z = int(np.floor(self.redshift/0.2))
        zcuts = (galdata.redshift >= zlow[i_z]) & (galdata.redshift < zhigh[i_z])

        all_cuts = np.where(
                envcuts & zcuts & 
                (galdata.mass > galdata.masslimit) & 
                (galdata.edgecut == 1) 
                ) 
        for key in galdata.__dict__.keys(): 
            setattr(self, key, getattr(galdata, key)[all_cuts])
        
        sfr_class = np.chararray(len(all_cuts[0]), itemsize=16)
        sfr_class[:] = 'star-forming'
        setattr(self, 'sfr_class', sfr_class)
        return None
        

def ObservedSFMS(observable, **kwargs): 
    ''' Calculate the observed SF-MS function relation for either
    Jeremy's group catalog or the SDSS+PRIMUS sample. The SF-MS 
    relation is calculated by the median SFR in mass bins. 
    '''
    if observable == 'groupcat': 
        if 'Mrcut' not in kwargs.keys(): 
            raise ValueError('Mrcut kwarg needs to be specified') 
        if 'position' not in kwargs.keys(): 
            raise ValueError('position kwarg needs to be specified') 
        # import Group Catalog
        galdata = GroupCat(**kwargs)
        galdata.Read()
        # classify galaxies into Starforming or quiescent based on  
        # SFRcut in gal_prop module 
        if 'isedfit' in kwargs.keys() and kwargs['isedfit']: 
            # use iSEDfit M*, SFR and SSFR instead 
            galdata._iSEDfitMatch()
            galdata.mass = galdata.iSEDfit_mass
            galdata.sfr = galdata.iSEDfit_sfr
        fq_obj = Fq()
        setattr(galdata, 
                'sfr_class', 
                fq_obj.Classify(galdata.mass, galdata.sfr, np.mean(galdata.z))
                )

    elif observable == 'sdssprimus': 
        if 'redshift' not in kwargs.keys():
            raise ValueError('Specify redshift bin') 
        if 'environment' not in kwargs.keys(): 
            raise ValueError('envrionment kwarg needs to be specifieid') 
        # import PRIMUS/SDSS galaxy catalog 
        galdata = PrimusSDSS(**kwargs)
        galdata.Read()

    else: 
        raise ValueError
    onlySF = (galdata.sfr_class == 'star-forming') & (galdata.sfr != -999.)
        
    # mass bins 
    mstar = np.arange(6.0, 12.0, 0.2) 
    m_low = mstar[:-1]
    m_high = mstar[1:]
    m_mid = [] 

    muSFR, sigSFR, ngal = [], [], [] 
    # loop thorugh mass bins 
    for i_m in range(len(mstar)-1): 
        # star forming galaxies within the mass bin 
        bin = np.where(
                (galdata.mass > m_low[i_m]) & 
                (galdata.mass <= m_high[i_m]) &
                onlySF
                )
        ngal_bin = len(bin[0]) 
        if ngal_bin < 10: 
            continue 
        
        ngal.append(ngal_bin)
        m_mid.append(mstar[i_m])
        muSFR.append(np.median(galdata.sfr[bin]) ) 
        sigSFR.append(np.std(galdata.sfr[bin]))

    return [np.array(m_mid), np.array(muSFR), np.array(sigSFR), np.array(ngal)]


def FitObservedSFMS(observable, Mfid=10.5, slope=None, mass_range=None, **kwargs): 
    ''' Fit the observed SFMS relation with a linear function.   

    Parameters 
    ----------
    observable : str
        String that specifies the observables. SDSS/PRIMUS or SDSS Group Catalog 
    Mfid : float
        Fiducial mass where the SFMS is pivoted around. A( M* - M_fid) + B. 
    '''
    # read in the observed SFMS
    mass, muSFR, sigSFR, ngal = ObservedSFMS(observable, **kwargs) 
    if mass_range is not None: 
        m_min, m_max = mass_range
        mlim = np.where((mass > m_min) & (mass < m_max))
        mass = mass[mlim]
        muSFR = muSFR[mlim]
        sigSFR = sigSFR[mlim]
        ngal = ngal[mlim]

    # amplitude of SFMS depends on z so the initial guess
    # depends on z 
    if observable == 'groupcat': 
        guess_z = 0.05
    elif observable == 'sdssprimus': 
        guess_z = kwargs['redshift'] 
    guess_z *= 0.7

    if slope is None: 
        p0 = [0.6, guess_z]
        fa = {'x': mass - Mfid, 'y': muSFR, 'err': sigSFR}
        bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, quiet=True) 
    else: 
        p0 = [guess_z]
        fa = {'x': mass - Mfid, 'y': muSFR, 'err': sigSFR, 'slope': slope}
        bestfit = mpfit.mpfit(util.mpfit_line_fixedslope, p0, functkw=fa, quiet=True) 
    return bestfit.params
