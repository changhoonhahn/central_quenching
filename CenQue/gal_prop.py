'''

Module for calculating properties of galaxy populations.
Properties include, Specific Star Formation Rate distribution, 
Quiescent Fraction as a funciton of mass, 


'''
import numpy as np 
import util.util as Util
from scipy.interpolate import interp1d

from sfr_evol import AverageLogSFR_sfms

# TreePM
from sham_hack import LFClass 
from sham_hack import SMFClass 


class Ssfr(object): 
    def __init__(self, **kwargs): 
        ''' Class object that describes the sSFR distribution of a 
        galaxy population
        '''
        self.kwargs = kwargs.copy()

        # mass bins 
        self.mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

        self.ssfr_range = [-13.0, -7.0]
        self.ssfr_nbin = 40
        
        self.ssfr_dist = None
        self.ssfr_bin_edges = None 
        self.ssfr_bin_mid = None 

    def Calculate(self, mass, ssfr): 
        ''' Calculate the SSFR distribution for the four hardcoded mass bins 
        from the mass and ssfr values 
        '''
        if len(mass) != len(ssfr): 
            raise ValueError

        self.ssfr_dist = [] 
        self.ssfr_bin_mid = [] 
        self.ssfr_bin_edges = [] 

        # loop through the mass bins
        for i_m, mass_bin in enumerate(self.mass_bins): 
            mass_lim = np.where(
                    (mass >= mass_bin[0]) & 
                    (mass < mass_bin[1])
                    )
            n_bin = len(mass_lim[0])

            # calculate SSFR distribution  
            dist, bin_edges = np.histogram(
                    ssfr[mass_lim], 
                    range=self.ssfr_range, 
                    bins=self.ssfr_nbin, 
                    normed=True)

            self.ssfr_dist.append(dist)
            self.ssfr_bin_mid.append(0.5 * (bin_edges[:-1] + bin_edges[1:]))
            self.ssfr_bin_edges.append(bin_edges)
        
        return [self.ssfr_bin_mid, self.ssfr_dist]

class Fq(object): 
    def __init__(self, **kwargs): 
        ''' Class that describes the quiescent fraction 
        '''
        self.kwargs = kwargs
        # mass bin 
        mb = np.arange(9.0, 12.0, 0.2)

        self.mass_low  = mb[:-1]
        self.mass_high = mb[1:]
        self.mass_mid  = 0.5 * (self.mass_low + self.mass_high) 
    
    def Calculate(self, mass=None, sfr=None, z=None, sfr_class=None, sfms_prop=None):
        ''' Calculate the quiescent fraction 
        '''
        if sfms_prop is None: 
            raise ValueError
        if sfr_class is None: 
            if sfr is None or z is None: 
                raise ValueError
            sfq = self.Classify(mass, sfr, z, sfms_prop=sfms_prop)
        else: 
            sfq = sfr_class  
            
        masses, f_q = [], [] 
        for i_m in xrange(len(self.mass_mid)):

            masslim = np.where(
                    (mass > self.mass_low[i_m]) & 
                    (mass <= self.mass_high[i_m]) 
                    )
            ngal_mass = len(masslim[0])
            if ngal_mass == 0:  # no galaxy in mass bin 
                continue 

            ngal_q = len(np.where(sfq[masslim] == 'quiescent')[0])
            masses.append( self.mass_mid[i_m] ) 
            f_q.append( np.float(ngal_q)/np.float(ngal_mass) )

        return [np.array(masses), np.array(f_q)]
    
    def Classify(self, mstar, sfr, z_in, sfms_prop=None):
        ''' Classify galaxies based on M*, SFR, and redshift inputs.
        Returns an array of classifications
        '''
        sfr_class = self.SFRcut(mstar, z_in, sfms_prop=sfms_prop)

        sfq = np.empty(len(mstar), dtype=(str,16))

        sf_index = np.where(sfr > sfr_class)
        sfq[sf_index] = 'star-forming'
        q_index = np.where(sfr <= sfr_class)
        sfq[q_index] = 'quiescent'

        return sfq 

    def SFRcut(self, mstar, zin, sfms_prop=None):
        ''' Specific SFR cut off used to classify SF or Quiescent 
        galaxies 
        ''' 
        #lowmass = np.where(mstar < 9.5)
        #factor = np.repeat(0.8, len(mstar))
        #factor[lowmass] = 1.0 
        #return -0.75 + 0.76*(zin-0.05) + 0.5*(mstar-10.5)
        #return -0.75 + 0.76*(zin-0.04) + factor*(mstar-9.5) - 0.8
        mu_sfr = AverageLogSFR_sfms(mstar, zin, sfms_prop=sfms_prop)
        #offset = -0.75
        offset = -0.9
        return mu_sfr + offset

    def model(self, Mstar, z_in, lit='cosmosinterp'):
        ''' Model quiescent fraction as a funcnction of 
        stellar mass and redshift from literature. Different methods 

        f_Q ( M_star, z) 

        Parameters
        ----------
        Mstar : array
            Galaxy stellar mass

        z_in : array
            Galaxy redshifts  

        lit : string
            String that specifies the model from literature 'cosmosinterp'
        '''

        if lit == 'cosmos_tinker': 
            qf_z0 = -6.04 + 0.63*Mstar
            
            try: 
                alpha = np.repeat(-2.57, len(Mstar))

                w2 = np.where((Mstar >= 10.) & (Mstar < 10.5))
                alpha[w2] = -2.52
                w3 = np.where((Mstar >= 10.5) & (Mstar < 11.))
                alpha[w3] = -1.47
                w4 = np.where((Mstar >= 11.) & (Mstar < 11.5))
                alpha[w4] = -0.55
                w5 = np.where(Mstar > 11.5)
                alpha[w5] = -0.12
            except TypeError: 
                if Mstar < 10.0: 
                    alpha = -2.57
                elif (Mstar >= 10.0) & (Mstar < 10.5): 
                    alpha = -2.52
                elif (Mstar >= 10.5) & (Mstar < 11.0): 
                    alpha = -1.47
                elif (Mstar >= 11.0) & (Mstar <= 11.5): 
                    alpha = -0.55
                elif (Mstar >= 11.5):
                    alpha = -0.12
            #else: 
            #    raise NameError('Mstar is out of range')

            output = qf_z0 * ( 1.0 + z_in )**alpha 
            try: 
                if output.min() < 0.0: 
                    output[np.where(output < 0.0)] = 0.0
                if output.max() > 1.0: 
                    output[np.where(output > 1.0)] = 1.0
            except TypeError:  
                if output < 0.0: 
                    output = 0.0
                elif output > 1.0: 
                    output = 1.0 

            return output 

        elif lit == 'cosmosinterp': 
            zbins = [0.36, 0.66, 0.88] 

            fq_z = [] 
            for zbin in zbins: 
                fq_file = ''.join([code_dir(), 'dat/wetzel_tree/', 
                    'qf_z', str(zbin), 'cen.dat' ]) 
               
                # read in mass and quiescent fraction
                mass, fq = np.loadtxt(fq_file, unpack=True, usecols=[0,1])  
                fq_z.append( np.interp(Mstar, mass, fq)[0] )   # interpolate to get fq(Mstar)
            interp_fq_z = interp1d(zbins, fq_z)#, kind='linear')#, fill_value='extrapolate') 
            if not isinstance(z_in, np.ndarray): 
                return interp_fq_z(np.array([z_in]))
            else: 
                return interp_fq_z(z_in)

        elif lit == 'cosmosfit': 
            zbins = [0.36, 0.66, 0.88] 
            exp_sigma = [1.1972271, 1.05830526, 0.9182575] 
            exp_sig = np.interp(z_in, zbins, exp_sigma) 
            output = np.exp( ( Mstar - 12.0 )/exp_sig)
            massive = np.where(Mstar > 12.0) 
            output[massive] = 1.0
            return output

        elif lit == 'wetzel':       # Wetzel et al. 2013
            qf_z0 = -6.04 + 0.63*Mstar
            
            try: 
                alpha = np.repeat(-2.3, len(Mstar))

                w1 = np.where((Mstar >= 9.5) & (Mstar < 10.0))
                alpha[w1] = -2.1
                w2 = np.where((Mstar >= 10.) & (Mstar < 10.5))
                alpha[w2] = -2.2
                w3 = np.where((Mstar >= 10.5) & (Mstar < 11.))
                alpha[w3] = -2.0
                w4 = np.where(Mstar >= 11.)
                alpha[w4] = -1.3
            except TypeError: 
                if Mstar < 9.5: 
                    alpha = -2.3
                elif (Mstar >= 9.5) & (Mstar < 10.0): 
                    alpha = -2.1
                elif (Mstar >= 10.0) & (Mstar < 10.5): 
                    alpha = -2.2
                elif (Mstar >= 10.5) & (Mstar < 11.0): 
                    alpha = -2.0
                elif (Mstar >= 11.0): # & (Mstar <= 11.5): 
                    alpha = -1.3
            #else: 
            #    raise NameError('Mstar is out of range')

            output = qf_z0 * ( 1.0 + z_in )**alpha 
            try: 
                if output.min() < 0.0: 
                    output[np.where(output < 0.0)] = 0.0
                if output.max() > 1.0: 
                    output[np.where(output > 1.0)] = 1.0
            except TypeError:  
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
            try: 
                if output.min() < 0.0: 
                    output[np.where(output < 0.0)] = 0.0
                if output.max() > 1.0: 
                    output[np.where(output > 1.0)] = 1.0
            except TypeError: 
                if output < 0.0: 
                    output = 0.0
                if output > 1.0: 
                    output = 1.0

            return output 

        elif lit == 'wetzel_alternate': 
            fqall = lambda A, alpha, z: A*(1.+z)**alpha
            fsat = lambda B0, B1, z: B0 + B1 * z
            fqsat = lambda C0, C1, M: C0 + C1 * M

            M_arr = np.array([9.75, 10.25, 10.75, 11.25]) 
            
            A_arr = np.repeat(0.227, len(M_arr)) 
            alpha_arr = np.repeat(-2.1, len(M_arr))
            B0_arr = np.repeat(0.33, len(M_arr))
            B1_arr = np.repeat(-0.055, len(M_arr))
            C0_arr = np.repeat(-3.26, len(M_arr))
            C1_arr = np.repeat(0.38, len(M_arr))
            
            # mass ranges
            w1 = np.where(M_arr < 10.0)
            w2 = np.where((M_arr >= 10.) & (M_arr < 10.5))
            w3 = np.where((M_arr >= 10.5) & (M_arr < 11.))
            w4 = np.where(M_arr >= 11.)

            A_arr[w1] = 0.227
            A_arr[w2] = 0.471
            A_arr[w3] = 0.775
            A_arr[w4] = 0.957
            alpha_arr[w1] = -2.1
            alpha_arr[w2] = -2.2
            alpha_arr[w3] = -2.0
            alpha_arr[w4] = -1.3
            B0_arr[w1] = 0.33
            B0_arr[w2] = 0.30
            B0_arr[w3] = 0.25
            B0_arr[w4] = 0.17
            B1_arr[w1] = -0.055 
            B1_arr[w2] = -0.073 
            B1_arr[w3] = -0.11 
            B1_arr[w4] = 0.1 

            
            out_arr = (fqall(A_arr, alpha_arr, z_in) - fqsat(C0_arr, C1_arr, M_arr) * fsat(B0_arr, B1_arr, z_in)) / (1 - fsat(B0_arr, B1_arr, z_in))
            out_arr = np.array([0.0]+list(out_arr))
            M_arr = np.array([9.0]+list(M_arr))

            interp_out = interp1d(M_arr, out_arr)
            within = np.where((Mstar > M_arr.min()) & (Mstar < M_arr.max())) 
            output = np.repeat(1.0, len(Mstar))
            output[within] = interp_out(Mstar[within])
            geha_lim = np.where(Mstar < 9.5) 
            output[geha_lim] = 0.

            massive_lim = np.where(Mstar > M_arr.max())
            extrap = lambda mm: (out_arr[-1] - out_arr[-2])/(M_arr[-1] - M_arr[-2]) * (mm - M_arr[-1]) + out_arr[-1]
            output[massive_lim] = extrap(Mstar[massive_lim])

            return output 
        else: 
            raise NameError('Not yet coded') 

def dFqdt(mstar, t_cosmic, lit='wetzelsmooth', dt=0.01): 
    ''' Calculate the derivative of the quiescent fraction, dfQ/dt, 
    as a function of cosmic time (t_cosmic) for some given parameteric 
    model of quiescent fraction evolution 

    Parameters
    ----------
    '''
    if isinstance(mstar, float): 
        mstar = np.array([mstar]) 

    qf = Fq() 
    fq_model = qf.model 
    try:  
        outofbound = np.where(t_cosmic > 13.1)
        t_cosmic[outofbound] = 13.1
    except TypeError: 
        if t_cosmic > 13.1: 
            t_cosmic = 13.1

    dfdt = (fq_model(mstar, Util.get_zsnap(t_cosmic+dt), lit=lit) - fq_model(mstar, Util.get_zsnap(t_cosmic), lit=lit))/dt

    neg = np.where(dfdt < 0.)
    dfdt[neg] = 0.
    return dfdt

class SMF(object): 
    def __init__(self, **kwargs): 
        self.kwargs = kwargs
        self.mass = None
        self.phi = None

    def Obj(self, obj, dlogm=None, box=None, m_arr=None, h=0.7, LF=False):
        ''' Calculate the SMF for any Class Object with attributes
        mass 

        Parameters
        ----------
        obj : 
            Any object with the attribute 'mass'
        dlogm : 
            log M bin size (default is 0.1)
        box : 
            Box length (default is 250 Mpc/h)
        '''
        if not LF: 
            return self._smf(obj.mass, dlogm=dlogm, box=box, m_arr=m_arr, h=h)
        else: 
            return self._smf(obj.mag_r, dlogm=dlogm, box=box, m_arr=m_arr, h=h)

    def analytic(self, redshift, dlogm='', source='li-drory-march'): 
        ''' Analytic SMF for a given redshift. 

        Return
        ------
        [masses, phi] : 
            array of masses, array of phi (smf number density) 
        '''
        if not dlogm: 
            dlogm = 0.1

        if redshift < 0.1:
            redshift = 0.1
    
        if source in ('cool_ages', 'blanton'): 
            m_arr = np.arange(-24., -16., dlogm)
            MF = LFClass(source=source, redshift=redshift) 
        else: 
            m_arr = np.arange(6.0, 12.1, dlogm)
            MF = SMFClass(source=source, redshift=redshift)
        
        mass, phi = [], [] 
        for mi in m_arr: 
            if source in ('cool_ages', 'blanton'): 
                mass.append(mi - 0.5 * dlogm)
                phi.append(MF.numden(-mi, -mi+dlogm)/dlogm) 
            else: 
                mass.append(mi + 0.5 * dlogm)
                phi.append(MF.numden(mi, mi+dlogm)/dlogm) 
        #print 'Analytic ', np.sum(np.array(phi))
        return np.array(mass), np.array(phi)

    def _smf(self, masses, dlogm=None, box=None, m_arr=None, h=0.7): 
        '''
        Calculate SMF given masses
        '''
        if not dlogm: 
            dlogm = 0.1 # log M bin size
        if not box: 
            box = 250   # 250 Mpc/h Box length
        
        if m_arr is None: 
            m_arr = np.arange(0.0, 12.1, dlogm)

        vol = box ** 3  # box volume
        mass, phi = [], [] 
        for mi in m_arr: 
            mass_bin = np.where(
                    (masses > mi) & 
                    (masses <= mi+dlogm)
                    )
            ngal_bin = np.float(len(masses[mass_bin]))
            
            mass.append(mi + 0.5 * dlogm)
            phi.append(np.float(ngal_bin)/vol/dlogm * h**3) 

        return np.array(mass), np.array(phi)
