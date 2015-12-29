'''

Central_Quenchign project


Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import h5py
import time
import os 
import json 
import matplotlib.pyplot as plt

import time
import random 
import numpy as np
import warnings

# --- Local ----
from smf import SMF
from ssfr import Ssfr
#from assign_sfr import assign_sfr 
from quiescent_fraction import get_fq
from util import cenque_utility as util
from util.mass_bins import simple_mass_bin
from central_subhalo import Subhalos
from central_subhalo import CentralSubhalos
from quiescent_fraction import sfq_classify
from sfms.fitting import get_param_sfr_mstar_z
from sfms.fitting import get_quiescent_mean_ssfr 

# -- plotting --
from plotting import plots 

class CenQue: 
    def __init__(self, **kwargs): 
        ''' Class that descirbes the central quenching (CenQue) galaxy 
        catalog
        '''
        self.kwargs = kwargs
    
        # TreePM properties
        self.mass = None
        self.sfr = None
        self.ssfr = None 
        self.parent = None
        self.child = None
        self.ilk = None
        self.snap_index = None
        self.pos = None     # position 
        self.halo_mass = None 
        self.sham_mass = None

        # star-formation properties
        self.gal_type = None    # quiescent/star-forming 
    
        # meta data
        if 'n_snap' in self.kwargs.keys(): 
            self.nsnap = self.kwargs['n_snap']
        else: 
            self.nsnap = None       # n_snapshot 

        if 'subhalo_prop' in self.kwargs.keys(): 
            self.subhalo_prop = self.kwargs['subhalo_prop']
        else: 
            self.subhalo_prop = None

        if 'sfr_prop' in self.kwargs.keys(): 
            self.sfr_prop = self.kwargs['sfr_prop']
        else: 
            self.sfr_prop = None

        self.zsnap = None           # z_snapshot
        self.t_cosmic = None    # t_cosmic for snapshot
        self.t_step = None      # t_cosmic step 

        self.mass_bins = self.set_mass_bins()
        self.data_columns = None 

        if 'cenque_type' in self.kwargs.keys(): 
            self.cenque_type = self.kwargs['cenque_type']
        else: 
            self.cenque_type = None

        self.ssfr_dist = None

    def import_treepm(self, nsnap, subhalo_prop=None): 
        ''' 
        Import Subhalo Snapshots with SHAM M* by interfacing with CentralSubhalos class 
        object. kwargs specify the parameters of the calculated SHAM subhalos. 
        
        Parameters
        ----------
        nsnap : (int)
            Snapshot number
        scatter : (float)
            Scatter between SMF and HMF in SHAM 
        source : (str)
            Specifies which analytic SMF to use for SHAM

        Notes 
        -----
        Also imports snapshot, redshift and t_cosmic metadata. 
        * SHAM stellar mass 
        * parent index 
        * child index 
        * halo mass 
        '''
        #start_time = time.time()
        if self.subhalo_prop is None: 
            self.subhalo_prop = subhalo_prop
        else: 
            if subhalo_prop:
                if self.subhalo_prop != subhalo_prop: 
                    raise ValueError
        
        centsub = CentralSubhalos()
        #centsub = Subhalos()
        centsub.read(nsnap, scatter=self.subhalo_prop['scatter'], source=self.subhalo_prop['source'])
        print centsub.file_name

        self.data_columns = [] 
        for col in centsub.data_columns: 
            if col == 'halo.m.max': 
                newcol = 'halo_mass'
            elif col == 'index': 
                newcol = 'snap_index'
            elif col == 'halo.m': 
                continue 
            else: 
                newcol = col
            setattr(self, newcol, getattr(centsub, col))
            self.data_columns.append(newcol)
        
        # Meta data 
        self.metadata = [ 'nsnap', 'zsnap', 't_cosmic', 't_step', 'cenque_type', 'subhalo_prop']
        # get snapshot redshift/cosmic time data using Andrew's table
        n_snaps, z_snap, t_snap, t_wid = np.loadtxt(
                'snapshot_table.dat', 
                unpack=True, 
                usecols=[0, 2, 3, 4]
                )
        self.nsnap = nsnap 
        self.zsnap = z_snap[(n_snaps.tolist()).index(nsnap)]        # redshift of snapshot
        self.t_cosmic = t_snap[(n_snaps.tolist()).index(nsnap)] # t_cosmic of snapshot 
        self.t_step = t_wid[(n_snaps.tolist()).index(nsnap)]    # Gyrs until next snapshot
        self.cenque_type = 'treepm_import'
        if self.subhalo_prop is None: 
            self.subhalo_prop = subhalo_prop
        else: 
            if subhalo_prop: 
                if self.subhalo_prop != subhalo_prop: 
                    raise ValueError

        #print "import_treepm takes ", time.time() - start_time
        return None

    def set_mass_bins(self): 
        """ Mass bin of 
        """
        massbins = simple_mass_bin() 
        return massbins

    def readin(self): 
        ''' Read in cenque data written by CenQue

        Get snapshot and column info from header then import column 
        into CenQue data 

        overly convoluted way to read in data -.-
        **UPDATED TO hdf5 FILE FORMAT**
        '''

        cq_file = self.file() 
        if os.path.isfile(cq_file): 
            pass
        else: 
            raise ValueError(cq_file+' does not exist') 

        f = h5py.File(cq_file, 'r') 
        #print 'Reading ', cq_file 

        grp = f['cenque_data']

        # read in meta data first
        for i_meta, metadatum in enumerate(grp.attrs.keys()): 
            if isinstance((grp.attrs.values())[i_meta], str):
                try: 
                    json_dict = json.loads((grp.attrs.values())[i_meta])
                    setattr(self, metadatum, json_dict) 
                except ValueError: 
                    setattr(self, metadatum, (grp.attrs.values())[i_meta]) 
            else:  
                setattr(self, metadatum, (grp.attrs.values())[i_meta]) 
        
        self.metadata = [str(key) for key in grp.attrs.keys()]

        for i_col, column in enumerate(grp.keys()): 
            setattr(self, column, grp[column][:])

        if self.cenque_type == 'treepm_import': 
            self.data_columns = [
                    'mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index', 'pos'
                    ]
        elif self.cenque_type == 'sf_assigned': 
            self.data_columns = [
                    'mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index', 'pos', 
                    'gal_type', 'sfr', 'ssfr'
                    ]
        else: 
            raise NotImplementedError
        
        #if 'sham_mass' not in grp.keys(): 
        #    setattr(self, 'sham_mass', grp['mass'][:])

        f.close() 

        return None 

    def writeout(self): 
        """ Write out self.data_columns to h5py file
        """
        
        output_file = self.file()  
        print 'Writing ', output_file 

        f = h5py.File(output_file, 'w')    
        grp = f.create_group('cenque_data')
    
        n_cols = len(self.data_columns)       # number of columns 
        col_fmt = []                # column format 
    
        for column in self.data_columns:      
            column_attr = getattr(self, column)

            # save to h5py data group 
            grp.create_dataset( column, data=column_attr ) 
        
        # save metadata 
        for metadatum in self.metadata: 
            if isinstance(getattr(self, metadatum), dict): 
                grp.attrs[metadatum] = json.dumps(getattr(self, metadatum))
            else: 
                grp.attrs[metadatum] = getattr(self, metadatum) 
        
        f.close()  
        return None 
    
    def file(self): 
        ''' 
        File name of CenQue object. Assigned based on the properties of the 
        CenQue object
        ''' 
        if self.nsnap is None: 
            raise ValueError()
        if self.cenque_type is None: 
            raise ValueError()
        if self.subhalo_prop is None: 
            raise ValueError()
     
        # Cenque Object with assigned starformation properties
        if self.cenque_type == 'treepm_import':
            file_spec_str = self._file_spec(subhalo_prop = self.subhalo_prop)
        elif self.cenque_type == 'sf_assigned':
            file_spec_str = self._file_spec(
                    subhalo_prop = self.subhalo_prop, 
                    sfr_prop = self.sfr_prop
                    )
        else: 
            raise NameError() 

        cenque_filename = ''.join([
            'dat/CenQue/', 
            'cenque_centrals.snapshot', 
            str(self.nsnap), file_spec_str, '.hdf5'
            ]) 

        return cenque_filename

    def _file_spec(self, subhalo_prop = None, sfr_prop  = None): 
        
        if subhalo_prop is None: 
            raise ValueError
    
        subhalo_str = ''.join([
            '_subhalo.scatter', str(subhalo_prop['scatter']), 
            '.source-', subhalo_prop['source']
            ])

        if sfr_prop is not None: 
            sfr_str = ''.join([
                '_sfrassign.fq-', sfr_prop['fq']['name'], 
                '.sfr-', sfr_prop['sfr']['name']
                ])
        else: 
            sfr_str = ''

        return ''.join([subhalo_str, sfr_str])

    def sample_trim(self, npwhere, quiet = False): 
        ''' Given numpy.where condition, apply numpy where condition
        to object data columns. 
        '''

        for column in self.data_columns:         
            obj_attr = getattr(self, column) 

            new_attr = obj_attr[npwhere]     

            setattr(self, column, new_attr)

        return None 

    # ----  SSFR -----
    def Ssfr(self, **ssfr_kwargs):          
        '''
        Calculate SSFR distribution 

        Returns 
        -------
        ssfr_bin_mid : 
            midpoints of SSFR bins
        ssfr_hist : 
            SSFR histogram 

        '''
        ssfr_dist = Ssfr()
        ssfr_bin_mid, ssfr_hist = ssfr_dist.cenque(self)
        
        #self.ssfr_dist = ssfr_dist  # save to class 

        return ssfr_bin_mid, ssfr_hist

    def plotSsfr(self, **pltkwargs):
        '''
        Plot SSFR of CenQue objectusing PlotSSFR plot object

        Parameters
        ----------
        pltkwargs : 
            - quenching : (True/False) If true plot stacked SSFR distribution 
                for CenQue data that differentiates the quiescent, quenching, 
                and star-fomring populations.
        '''
        plt.close() # in case there's another plot

        ssfr_plot = plots.PlotSSFR()
        ssfr_plot.cenque(self, **pltkwargs)

        if 'groupcat' in pltkwargs.keys(): 
            if pltkwargs['groupcat']: 
                ssfr_plot.groupcat(Mrcut=18)    # Mr cut is hardcoded
            else: 
                raise ValueError('Are you sure you wanted to specify groupcat=False?')

        ssfr_plot.set_axes()

        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                ssfr_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')
            return None
        else: 
            return ssfr_plot

    # --- SMF ---
    def SMF(self, **smf_kwargs):
        ''' 
        Calculate the Stellar Mass Function of Cenque Object
        '''
        if self.mass is None: 
            raise ValueError
        smf = SMF()
        return smf.cenque(self, **smf_kwargs)

    def plotSMF(self, **pltkwargs): 
        '''
        Plot SMF of CenQue object using PlotSMF object 
        '''
        plt.close() # in case there's another plot

        smf_plot = plots.PlotSMF()
        smf_plot.cenque(self)
        if self.zsnap < 0.1: 
            redshift = 0.1
        else: 
            redshift = self.zsnap
        smf_plot.analytic(redshift, source=self.subhalo_prop['source'])
        smf_plot.central_subhalo(self.nsnap, subhalo_prop=self.subhalo_prop)
        smf_plot.set_axes()

        if 'show' in pltkwargs.keys(): 
            smf_plot.show()
        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                smf_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')
            return None
        else: 
            return smf_plot

    # --- Quiescent Fraction ---
    def Fq(self, **fq_kwargs):
        '''
        Calculate quiescent fraction of CenQue class object
        '''
        if self.zsnap is None: 
            raise ValueError
        if self.mass is None: 
            raise ValueError
        if self.sfr is None: 
            raise ValueError

        # Star-forming or Quiescent    
        sf_q = sfq_classify(self.mass, self.sfr, self.zsnap, Mrcut=18)

        masses, f_q = [], [] 
        for i_bin in xrange(self.mass_bins.nbins): 

            in_massbin = np.where( 
                    (self.mass > self.mass_bins.mass_low[i_bin]) &
                    (self.mass <= self.mass_bins.mass_high[i_bin])
                    )

            n_massbin = len(in_massbin[0])

            if n_massbin == 0: 
                continue

            q = len(np.where(sf_q[in_massbin] == 'quiescent')[0])
            
            masses.append( self.mass_bins.mass_mid[i_bin] ) 
            f_q.append( np.float(q)/np.float(n_massbin) )
        
        return np.array(masses), np.array(f_q)
    
    def plotFq(self, **pltkwargs): 
        '''
        Plot f_q of CenQue object using PlotFq plot object
        '''
        plt.close() # in case there's another plot

        fq_plot = plots.PlotFq() 
        fq_plot.cenque(self, **pltkwargs)

        if 'param' in pltkwargs.keys():
            # hardcoded for convenience 
            if pltkwargs['param']: 
                fq_plot.param_fq(line_color='k', label=None)

        fq_plot.set_axes()

        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                fq_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return fq_plot

    def plotTau(self, tau_prop, **pltkwargs): 
        '''
        Plot quenching timescale 
        '''
        plt.close() 
        tau_dict = tau_prop.copy()
        tau_plot = plots.PlotTau(tau_dict)

        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                tau_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return tau_plot 

    def plotMstarMhalo(self, **pltkwargs): 
        '''
        Plot M* versus M_halo relationship

        Notes
        -----
        Uses bovy_plot.scatter_plot so things are a big clunky 
        '''

        plt.close()

        mm_plot = plots.PlotMstarMhalo()
        mm_plot.cenque(self, **pltkwargs)
        
        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                mm_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return mm_plot 

    def plotSFMS(self, **pltkwargs): 
        '''
        Plot the Star Forming Main Sequence of the CenQue object
        
        Notes
        -----
        Uses bovy_plot.scatter_plot so things are a big clunky 
        '''
        plt.close() 

        sfms_plot = plots.PlotSFMS()
        sfms_plot.cenque(self, **pltkwargs)
        
        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                sfms_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return sfms_plot 
    
    def _assign_sfr(self, sfr_prop = None, quiet=True, **kwargs):
        ''' 
        ****No longer supported****
        Assign Star Formation Rates to CenQue class object based on 
        properties listed in keyword sfr_prop, which specifies the 
        analytic quiescsent fraction prescription and the SFR assignment
        method. Namely, 'gal_type', 'sfr', and 'ssfr' are assigned to 
        M* and M_halo values.

        Goes through mass bins and then classifies unassigned 
        galaxies into quiescent/star-forming based on quiescent fraction 
        function. 

        Then SF properties are assigned to the galaxies. For quiescent 
        galaxies the sSFR is drawn from a log normal distribution about 
        some predeterimined mu_sSFR. For star-forming galaxies, SFR is
        sampled from a designated SF-MS model. 

        Parameters
        ----------
        sfr_prop : (dict)
            Dictionary that specifies the SFR assignment properties and the 
            quiescent fraction properties used in the method  

        '''
        # time the code 
        start_time = time.time()

        if self.mass == None: 
            raise ValueError()
        if self.zsnap == None: 
            raise ValueError()
        
        if self.sfr_prop is None: 
            self.sfr_prop = sfr_prop
            self.metadata.append('sfr_prop')
        else: 
            if sfr_prop: 
                if self.sfr_prop != sfr_prop: 
                    raise ValueError('SFR properties do not match')
        fq_dict = self.sfr_prop['fq']       # dictionary that specifies fq properties
        sfr_dict = self.sfr_prop['sfr']     # dictionary that specifies sfr properties

        mass_bin_low  = self.mass_bins.mass_low
        mass_bin_mid  = self.mass_bins.mass_mid
        mass_bin_high = self.mass_bins.mass_high

        #keep = np.where(
        #        (self.mass > mass_bin_low.min()) & 
        #        (self.mass <= mass_bin_high.max()) & 
        #        (self.child >= 0) 
        #        ) 
        #ngal_keep = len(keep[0])
        #if (np.min(self.mass) < mass_bin_low.min()) or (np.max(self.mass) > mass_bin_high.max()): 
        #    self.sample_trim(keep) 
        ngal_keep = len(self.mass)
        
        for attrib in ['gal_type', 'sfr', 'ssfr']: 
            if attrib not in self.data_columns: 
                self.data_columns.append(attrib)

        if self.gal_type is None:         
            self.gal_type = np.repeat('', ngal_keep).astype('|S16') 
            self.sfr      = np.repeat(-999., ngal_keep) 
            self.ssfr     = np.repeat(-999., ngal_keep)

        # analytic mu_SFR(M*, z) with randomly sampled normal delta_SFR. 
        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

        if 'delta_sfr' not in self.__dict__.keys(): 
            self.delta_sfr = np.repeat(-999., ngal_keep)
        if 'avg_sfr' not in self.__dict__.keys(): 
            self.avg_sfr = np.repeat(-999., ngal_keep)
        for attrib in ['avg_sfr', 'delta_sfr']: 
            if attrib not in self.data_columns: 
                self.data_columns.append(attrib)

        # f_Q(M_*mid, z_snapshot) 
        qf_massbin = get_fq(mass_bin_mid, self.zsnap, lit = fq_dict['name']) 

        unassigned = [ 
                np.where(
                    (self.mass > mass_bin_low[i_m]) &
                    (self.mass <= mass_bin_high[i_m]) &
                    (self.gal_type == ''))[0]
                for i_m in xrange(self.mass_bins.nbins)
                ]
        # Ngal(M_mid), Ngal,Q(M_mid) and Ngal,SF(M_mid)
        ngal_massbin = np.array([x.size for x in unassigned])
        ngal_q_massbin = np.rint(qf_massbin * ngal_massbin.astype(float)).astype(int)
        ngal_sf_massbin = ngal_massbin - ngal_q_massbin

        # fail-safe for ngal_q_massbin
        if len(np.where(ngal_q_massbin > ngal_massbin)[0]) > 0: 
            ngal_q_massbin[np.where(ngal_q_massbin > ngal_massbin)] = \
                    ngal_massbin[np.where(ngal_q_massbin > ngal_massbin)]

        for i_m in xrange(self.mass_bins.nbins):     # loop through mass bins
            if ngal_massbin[i_m] == 0: 
                continue
            #begin_loop_time = time.time()
            if not quiet: 
                print mass_bin_low[i_m], ' < M < ', mass_bin_high[i_m]
                print 'fQ = ', qf_massbin[i_m], ' Ngal = ', ngal_massbin[i_m]
                print 'Ngal,Q = ', ngal_q_massbin[i_m], ' Ngal,SF = ', ngal_sf_massbin[i_m]

            shuffled = np.arange(ngal_massbin[i_m])
            np.random.seed()
            np.random.shuffle(shuffled)
            i_q_end = ngal_q_massbin[i_m]
    
            # Randomly select ngal_q_massbin quiescent galaxies from the 
            # massbin. Assign gal_type = 'quiescent', sSFR and SFR 
            # based on log-normal distribution of SSFR with 0.18 dex scatter
            if i_q_end > 0: 
                i_q_massbin = unassigned[i_m][shuffled[:i_q_end]]
                self.gal_type[i_q_massbin] = 'quiescent'   

                mu_q_ssfr = util.get_q_ssfr_mean(self.mass[i_q_massbin]) 
                if len(i_q_massbin) != ngal_q_massbin[i_m]:
                    raise ValueError
                self.ssfr[i_q_massbin] = 0.18 * np.random.randn(ngal_q_massbin[i_m]) + mu_q_ssfr 
                self.sfr[i_q_massbin]  = self.ssfr[i_q_massbin] + self.mass[i_q_massbin]
    
            # ngal_sf_massbin starforming galaxies from the massbin. Assign 
            # them 'star-forming' gal_type and sSFR and SFR in some manner
            if ngal_sf_massbin[i_m] > 0: 
                i_sf_massbin = (unassigned[i_m])[shuffled[i_q_end:]]
                self.gal_type[i_sf_massbin] = 'star-forming'
                
                mu_sf_sfr = sfr_mstar_z(self.mass[i_sf_massbin], self.zsnap)
                if  sfr_dict['name'] == 'average': 
                    sigma_sf_sfr = sig_sfr_mstar_z(self.mass[i_sf_massbin], self.zsnap)
                elif sfr_dict['name'] == 'average_noscatter': 
                    sigma_sf_sfr = 0.0
                else:
                    raise NotImplementedError
                    
                self.avg_sfr[i_sf_massbin] = mu_sf_sfr
                self.delta_sfr[i_sf_massbin] = sigma_sf_sfr * np.random.randn(ngal_sf_massbin[i_m])
                self.sfr[i_sf_massbin] = mu_sf_sfr + self.delta_sfr[i_sf_massbin]

                if not quiet:
                    print 'Starforming galaxies: '
                    print 'Average(SFR) = ',np.mean(mu_sf_sfr),' sigma(SFR) = ',sigma_sf_sfr 

            self.ssfr[i_sf_massbin] = self.sfr[i_sf_massbin] - self.mass[i_sf_massbin]

        # double check that SF assign didn't fail anywhere
        #assign_fail = np.where((self.sfr == -999.0) | (self.ssfr == -999.0))
        #if len(assign_fail[0]) > 0: 
        #    raise ValueError('Function failed!')

        print 'Assign SFR function takes', (time.time()-start_time), ' seconds'

        if 'evol_from' in self.cenque_type: 
            pass
        else: 
            self.cenque_type = 'sf_assigned'

        return None

def AssignSFR(mass, redshift, sfr_prop=None):
    ''' Assign star-forming properties based on the input mass and redshift. 
    Return gal_type, SFR, and sSFR given a set of masses and redshifts.

    Parameters
    ----------
    mass : 
        Stellar mass
    redshift : 
        Corresponding redshift
    sfr_prop : 
        Dictionary that specifies the star-forming properties.

    Notes
    -----
    * Takes ~1.5 seconds for assign_sfr_ancestor in lineage object
    '''
    # time the code 
    start_time = time.time()
    if sfr_prop is None: 
        raise ValueError('Specify SFR Properties dictionary')
    fq_dict = sfr_prop['fq']       # dictionary that specifies fq properties
    sfr_dict = sfr_prop['sfr']     # dictionary that specifies sfr properties

    ngal = len(mass)
    gal_type = np.repeat('', ngal).astype('|S16') 
    sfr      = np.repeat(-999., ngal) 
    ssfr     = np.repeat(-999., ngal)
    if sfr_dict['name'] == 'average': 
        # analytic mu_SFR(M*, z) with randomly sampled normal delta_SFR. 
        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()
        delta_sfr   = np.repeat(-999., ngal)
        avg_sfr     = np.repeat(-999., ngal)
    else: 
        raise NotImplementedError

    # massive galaxies
    massive = np.where(mass > 0.0)[0]
    ngal_massive = len(massive)

    # f_Q(M*, z) for each massive galaxy
    qf = get_fq(mass[massive], redshift[massive], lit=fq_dict['name']) 

    # determine quiescent/star-forming based on random's 
    # comparison to quiescent fraction
    np.random.seed()
    rand = np.random.uniform(0., 1., ngal_massive)
    quiescent = np.where(rand <= qf)
    starforming = np.where(rand > qf)
    ngal_q = len(quiescent[0])
    ngal_sf = len(starforming[0])

    # quiescent 
    gal_type[massive[quiescent]] = 'quiescent'
    #mu_q_ssfr = util.get_q_ssfr_mean(mass[massive[quiescent]]) 
    mu_q_ssfr = get_quiescent_mean_ssfr(mass[massive[quiescent]])
    ssfr[massive[quiescent]] = 0.18 * np.random.randn(ngal_q) + mu_q_ssfr 
    sfr[massive[quiescent]]  = ssfr[massive[quiescent]] + mass[massive[quiescent]]
    
    # star-forming 
    gal_type[massive[starforming]] = 'star-forming'
    mu_sf_sfr = sfr_mstar_z(mass[massive[starforming]], redshift[massive[starforming]])
    if sfr_dict['name'] == 'average': 
        sigma_sf_sfr = sig_sfr_mstar_z(mass[massive[starforming]], redshift[massive[starforming]])
    else:
        raise NotImplementedError
    avg_sfr[massive[starforming]] = mu_sf_sfr
    delta_sfr[massive[starforming]] = sigma_sf_sfr * np.random.randn(ngal_sf)
    sfr[massive[starforming]] = mu_sf_sfr + delta_sfr[massive[starforming]]
    ssfr[massive[starforming]] = sfr[massive[starforming]] - mass[massive[starforming]]
    
    print 'Assign SFR took ', time.time() - start_time  
    if sfr_dict['name'] == 'average': 
        return [gal_type, sfr, ssfr, delta_sfr, avg_sfr] 
    else: 
        return [gal_type, sfr, ssfr]

"""
    def build_cenque_importsnap(snapshots, scatter = 0.0): 
        ''' Build CenQue snapshots with TreePM data imported
        
        Notes
        -----
        * pretty much hardcoded
        '''
        for i_snap in snapshots: 
            snap = CenQue() 
            snap.import_treepm(i_snap, scatter=scatter)
            snap.writeout()

    def build_cenque_original(snapshot, **kwargs):
        snap = CenQue() 
        snap.AssignSFR(i_snap, **kwargs) 
        snap.writeout(nsnap=i_snap, file_type='sf assign', **kwargs)

    if __name__=='__main__': 
        #build_cenque_importsnap(range(1,21), scatter = 0.2)

        blah = CenQue()
        blah.import_treepm(20, scatter=0.2)
        blah = assign_sfr(blah)
        blah.writeout()

    #EvolveCenQue(13, 1, fqing_yint=-5.84, tau='instant')  
    #tau='linefit', tau_param=[-0.5, 0.4]) 
    #EvolveCenQue(13, 1, fqing_yint=-5.84, tau='linefit', tau_param=[-0.4, 0.2])
    #build_cenque_importsnap() 
    #build_cenque_original(sfr='sfr_func') 
    #EvolveCenQue(13, 1, tau='linefit', tau_param=[-0.7, 0.4], 
    #        sfr='sfr_func', stellmass='sham') 
    #EvolveCenQue(13, 1, tau='linefit', tau_param=[-0.7, 0.4], 
    #        sfr='sfr_func', stellmass='integrated') 
"""
