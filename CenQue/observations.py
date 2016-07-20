'''


Module for observed galaxy catalogs with galactic properties. 
Jeremy Tinker's Group Catalog along with John Moustakas's 
iSEDfit SDSS + PRIMUS catalogs with my environment values. 


'''
import os 
import h5py
import time 
import numpy as np
from scipy import interpolate
from pydl.pydlutils.spheregroup import spherematch

from util import mpfit 
from util import util 
from util.util import code_dir

from gal_prop import Fq
from gal_prop import Ssfr

import matplotlib.pyplot as plt
from ChangTools.fitstables import mrdfits 


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
            code_dir(), 'dat/observations/', 
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

    def _ReadGroupCat_GalData(self): 
        '''
        '''
        h = 0.7
        
        galdata_file = ''.join([
            code_dir(), 'dat/observations/', 
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

    def _ReadGroupCat_Prob(self): 
        ''' wrapper for reading in probability galaxy data 
        '''
        file = ''.join([
            code_dir(), 'dat/observations/', 
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
        mfdata = mrdfits(code_dir(), 'dat/observations/mfdata_all_supergrid01_sdss.fits.gz') 
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
        self.redshift_bin = redshift
        self.environment = environment 

    def Read(self): 
        '''
        '''
        sfr_class = [] 
        for sfq in ['star-forming', 'quiescent']:
            file_dir = code_dir()+'dat/observations/envcount/'
            if sfq == 'star-forming': 
                sfq_str = 'active'
            else: 
                sfq_str = sfq 

            sdss_file = ''.join([file_dir, 
                    'envcount_cylr2.5h35_thresh75_sdss_', sfq_str, '_z0.05_0.12_primuszerr.fits']) 
            primus_file = ''.join([file_dir, 
                    'envcount_cylr2.5h35_thresh75_', sfq_str, '_z0.2_1.0_lit.fits']) 
            # redshift determines SDSS or PRIMUS sample
            if self.redshift_bin < 0.2: 
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
            i_z = int(np.floor(self.redshift_bin/0.2))
            zcuts = (galdata.redshift >= zlow[i_z]) & (galdata.redshift < zhigh[i_z])

            all_cuts = np.where(
                    envcuts & zcuts & 
                    (galdata.mass > galdata.masslimit) & 
                    (galdata.edgecut == 1) 
                    ) 
            for key in galdata.__dict__.keys(): 
                try: 
                    setattr(self, key, 
                            np.concatenate([getattr(self,key), getattr(galdata, key)[all_cuts]]))
                except AttributeError: 
                    setattr(self, key, getattr(galdata, key)[all_cuts])
            
            sfr_class_tmp = np.chararray(len(all_cuts[0]), itemsize=16)
            sfr_class_tmp[:] = sfq

            sfr_class.append(sfr_class_tmp)

        setattr(self, 'sfr_class', np.concatenate(sfr_class))

        return None
        
def ObservedSFMS(observable, sfms_prop=None, **kwargs): 
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
                fq_obj.Classify(galdata.mass, galdata.sfr, np.mean(galdata.z), 
                    sfms_prop=sfms_prop) 
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

def ObservedSSFR(observable, **kwargs): 
    ''' Calculate the observed SSFR function relation for either
    Jeremy's group catalog or the SDSS+PRIMUS sample. 
    '''
    if 'groupcat' in observable: 
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

    mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

    ssfr_range = [-13.0, -7.0]
    ssfr_nbin = 40

    ssfr = galdata.sfr - galdata.mass

    ssfr_dist, ssfr_bin_mid, ssfr_bin_edges = [], [], [] 
    # loop through the mass bins
    for i_m, mass_bin in enumerate(mass_bins): 
        mass_lim = np.where(
                (galdata.mass >= mass_bin[0]) & 
                (galdata.mass < mass_bin[1])
                )
        n_bin = len(mass_lim[0])

        # calculate SSFR distribution  
        dist, bin_edges = np.histogram(
                ssfr[mass_lim], 
                range=ssfr_range, 
                bins=ssfr_nbin, 
                normed=True)

        ssfr_dist.append(dist)
        ssfr_bin_mid.append(0.5 * (bin_edges[:-1] + bin_edges[1:]))
        ssfr_bin_edges.append(bin_edges)
    
    return [ssfr_bin_mid, ssfr_dist]

def ObservedSSFR_Peaks(observable, sfq='star-forming', **kwargs): 
    ''' Estimate the positions of the SF and Q peaks of the SSFR distribution. 
    '''
    # calculate the SSFR distribution 
    ssfr_bin_mid, ssfr_dist = ObservedSSFR(observable, **kwargs) 
    # mass bins of the SSFR distribution 
    mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

    Mass = np.zeros(len(mass_bins))
    SSFR = np.zeros(len(mass_bins))
    SFR = np.zeros(len(mass_bins))
    for i_mass, mbin in enumerate(mass_bins): 
        if sfq == 'star-forming': 
            # estimate of the SF peak 
            SFQpeak = np.where(ssfr_bin_mid[i_mass] > -11.5) 
        elif sfq == 'quiescent': 
            SFQpeak = np.where(ssfr_bin_mid[i_mass] < -11.5) 
        SFQdist = ssfr_dist[i_mass][SFQpeak]
        Mass[i_mass] = np.mean(mbin)
        SSFR[i_mass] = (ssfr_bin_mid[i_mass][SFQpeak])[SFQdist.argmax()]  
        SFR[i_mass] = (ssfr_bin_mid[i_mass][SFQpeak])[SFQdist.argmax()] + np.mean(mbin)

    output_dict = {
            'sfr_class': sfq,
            'mass': Mass, 
            'ssfr': SSFR, 
            'sfr': SFR
            }
    return output_dict 

def FitObservedSSFR_Peaks(observable='groupcat', sfq='star-forming', Mrcut=18, position='central', Mfid=10.5): 
    ''' Fit the (M*, SFR) positions of the SSFR distribution SF peak to a
    SFR(M*) linear parameterization. This is only for the SDSS group catalog.  
    '''
    if 'groupcat' not in observable: 
        raise ValueError

    SSFRpeak = ObservedSSFR_Peaks(observable, sfq=sfq, Mrcut=Mrcut, position=position)

    SFRs = SSFRpeak['sfr']
    Mass = SSFRpeak['mass']
    
    # remove outliers
    if sfq  == 'star-forming': 
        # for the SDSS group catalog, the last mass bin has nearly no SF galaxies 
        # so the estimate is inaccurate 
        SFRs = SFRs[:-1]
        Mass = Mass[:-1]
        SFRs[0] = -0.4
        SFRs[1] = -0.275
        p0 = [0.6, 0.]  # guess
    elif sfq == 'quiescent': 
        SFRs[0] = -1.925
        SFRs[-1] = -1.475
        p0 = [-0.4, 0.] 
    print SFRs

    fa = {'x': Mass - Mfid, 'y': SFRs}
    bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, quiet=True) 

    return bestfit.params


def Lee2015_SFMS_zslope(): 
    ''' Calculate the sloep of the redshift dependence of the Lee et al. (2015) SFMS parameterizations 

    Lee et al. (2015) parameterized SFMS: 

    log SFR(M*, z) = S0(z) - log( 1 + (M*/M0)^-gamma)
    
    S0(z) quantifies the redshift dependence of the SFMS. Lee et al. (2015) fits S0 for each redshift bin. 
    I fit:

    S0(z) = A_z * ( z - 0.0502) + C

    returns [A_z, C]
    '''
    z_mid = np.array([0.36, 0.55, 0.70, 0.85, 0.99, 1.19])
    S0 = np.array([0.80, 0.99, 1.23, 1.35, 1.53, 1.72])
    S0_err = np.array([0.019, 0.015, 0.016, 0.014, 0.017, 0.024])

    p0 = [1.5, 0.0]
    fa = {'x': z_mid-0.05, 'y': S0, 'err': S0_err}
    bestfit = mpfit.mpfit(util.mpfit_line, p0, functkw=fa, quiet=True) 

    return bestfit.params


def FqCen_bestfit(clobber=False): 
    ''' We parameterize the central galaxy quiescent fraction as 
    fQ(M*, z) = fQ(M*, z = 0.) * (1+z)^(alpha(M*)). 

    Here we fit for alpha

    log fQ(M*,z) - log fQ(M*, z=0) = alpha(M*) * log(1+z) 
    '''
    # mass binnning we impose 
    m_low = np.array([9.5, 10., 10.5, 11., 11.5]) 
    m_high = np.array([10., 10.5, 11., 11.5, 12.0])
    m_mid = 0.5 * (m_low + m_high) 

    # SDSS 
    fq_file = ''.join(['dat/observations/cosmos_fq/', 'fcen_red_sdss_scatter.dat']) 
    m_sdss, fqcen_sdss, N_sdss = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2])
    
    fqcen_sdss_rebin = [] 
    for im, m_mid_i in enumerate(m_mid): 
        sdss_mbin = np.where(
                (m_sdss >= m_low[im]) & 
                (m_sdss < m_high[im])) 
    
        fqcen_sdss_rebin.append(
                np.sum(fqcen_sdss[sdss_mbin] * N_sdss[sdss_mbin].astype('float'))/np.sum(N_sdss[sdss_mbin].astype('float'))
                )
    fqcen_sdss_rebin = np.array(fqcen_sdss_rebin)
    
    fqcen_cosmos_rebin = [] 
    fqcen_low_cosmos_rebin = [] 
    fqcen_high_cosmos_rebin = [] 

    for iz, z in enumerate([0.36, 0.66, 0.88]): 
        fq_file = ''.join(['dat/observations/cosmos_fq/', 
            'stats_z', str(iz+1), '.fq_cen']) 
        
        m_cosmos, fqcen_cosmos, fqcen_cosmos_low, fqcen_cosmos_high = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2,3])
        m_cosmos = np.log10(m_cosmos)

        fqcen_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos) 
        fqcen_low_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos_low) 
        fqcen_high_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos_high) 
    
        fqcen_cosmos_rebin.append(fqcen_interp(m_mid))
        fqcen_low_cosmos_rebin.append(fqcen_low_interp(m_mid))
        fqcen_high_cosmos_rebin.append(fqcen_high_interp(m_mid))

    norm_fqcen_cosmos_rebin = np.array([
            fqcen_cosmos_rebin[ii]/fqcen_sdss_rebin for ii in range(len(fqcen_cosmos_rebin))
            ])
    norm_fqcen_low_cosmos_rebin = np.array([
            fqcen_low_cosmos_rebin[ii]/fqcen_sdss_rebin for ii in range(len(fqcen_cosmos_rebin))
            ])
    norm_fqcen_high_cosmos_rebin = np.array([
            fqcen_high_cosmos_rebin[ii]/fqcen_sdss_rebin for ii in range(len(fqcen_cosmos_rebin))
            ])

    z_arr = np.array([0.0, 0.36, 0.66, 0.88])
    alpha_m = [] 
    for im, mm in enumerate(m_mid): 
        # fit the redshift evolution with a power law for each mass bin 
        # using mpfit with asymmetric errors from Tinker et al. (2013)'s 
        # errors
        p0 = [-2.] 
        fa = {
                'x': z_arr, 
                'y': np.log10(np.array([1.] + list(norm_fqcen_cosmos_rebin[:,im]))), 
                'y_low': np.log10(np.array([.999] + list(norm_fqcen_low_cosmos_rebin[:,im]))), 
                'y_high': np.log10(np.array([1.001] + list(norm_fqcen_high_cosmos_rebin[:,im])))
                }
        bestfit = mpfit.mpfit(mpfit_z_powerlaw, p0, functkw=fa, quiet=True) 
        alpha_m.append(bestfit.params[0])

    output_file = ''.join([os.path.dirname(os.path.realpath(__file__)).split('CenQue')[0], 
        'dat/fqcen_alphaM.dat']) 
    if not os.path.isfile(output_file) or clobber:
        np.savetxt(output_file, np.array([m_mid, alpha_m]).T, fmt=['%10.5f', '%10.5f']) 
    
    return [m_mid, alpha_m]


def z_powerlaw(z, p): 
    # just line function 
    return p[0] * np.log10(1. + z) 

def mpfit_z_powerlaw(p, fjac=None, x=None, y=None, y_low=None, y_high=None): 
    model = z_powerlaw(x, p) 
    status = 0 
    
    err_low = y - y_low  
    err_high = y_high - y
    
    resid = y-model

    err = err_low 
    positive = np.where(resid >= 0) 
    err[positive] = err_high 
    return([status, (y-model)/err]) 
