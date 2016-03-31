'''


Module for galaxy populations 

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
#from smf import SMF
#from ssfr import Ssfr
from gal_prop import Fq
from gal_prop import SMF
from gal_prop import Ssfr
#from assign_sfr import assign_sfr 
from util import cenque_utility as util

from sfr_evol import AverageLogSFR_sfms
from sfr_evol import ScatterLogSFR_sfms
from sfr_evol import AverageLogSSFR_q_peak
from sfr_evol import ScatterLogSSFR_q_peak

from central_subhalo import Subhalos
from central_subhalo import CentralSubhalos
#from sfms.fitting import get_quiescent_mean_ssfr 

# -- plotting --
from plotting import plots 

class CGPop(object): 
    def __init__(self, **kwargs): 
        ''' Class object that descirbes the central galaxy population derived 
        from CentralSubhalos class, which is a wrapper for Andrew Wetzel's TreePM 
        Subhalo Catalog. 
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

        self.data_columns = None    # default has no data

        self.ssfr_dist = None

    def ImportSubhalo(self, nsnap, subhalo_prop=None): 
        ''' Import Central Subhalo Snapshots with SHAM M* through the CentralSubhalos class 
        object. The CentralSubhalos class is a wrapper for Andrew Wetzel's TreePM code.
        
        Parameters
        ----------
        nsnap : int
            Subhalo Snapshot number. Each snapshot represents a different redshift. 
            Higher number are at higher redshifts. 
        
        subhalo_prop : dict (optional) 
            Optional dictionary that specifies the subhalo snapshot properties. 
            The properties include, 
            - scatter : float
                Scatter between SMF and HMF in SHAM 
            - source : str
                Specifies which analytic SMF to use for SHAM

        Notes 
        -----
        Also imports snapshot, redshift and t_cosmic metadata. 
        * SHAM stellar mass 
        * parent index 
        * child index 
        * halo mass 
        '''
        if self.subhalo_prop is None: 
            self.subhalo_prop = subhalo_prop
        else: 
            if subhalo_prop and self.subhalo_prop != subhalo_prop: 
                # subhalo_prop values do not match 
                raise ValueError
        
        centsub = CentralSubhalos()
        centsub.Read(
                nsnap,  # snapshot 
                scatter=self.subhalo_prop['scatter'],   # scatter
                source=self.subhalo_prop['source'],      # SMF source
                nsnap_ancestor=self.subhalo_prop['nsnap_ancestor']
                )
        #print centsub.file_name
    
        # pass CentralSubhalos class attributes to 
        self.data_columns = self.get_datacol() 
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
            if newcol not in self.data_columns: 
                self.data_columns.append(newcol)

        # Meta data 
        self.metadata = [ 'nsnap', 'zsnap', 't_cosmic', 't_step', 'subhalo_prop']
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
        #print "import_treepm takes ", time.time() - start_time
        return None

    def Read(self): 
        ''' Read in data written specifically by the/for this class. That is in 
        hdf5 file format. The metadata for the data is in grp attrs.
        '''
        cgpop_file = self.File()       # file name of CGPop object
        if not os.path.isfile(cgpop_file): 
            raise ValueError(cgpop_file+' does not exist') 

        f = h5py.File(cgpop_file, 'r') 
        grp = f['cenque_data']
        # first read in meta data, which are saved as attrs 
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

        self.data_columns = self.get_datacol() 
        for i_col, column in enumerate(grp.keys()): 
            setattr(self, column, grp[column][:])
            if not isinstance(column, str):
                column = str(column).encode('ascii')
            if column not in self.data_columns: 
                self.data_columns.append(column)
        #print 'CGPop data columns', self.data_columns
        f.close() 
        return None 

    def Write(self): 
        ''' Write out data_columns to h5py file along with the appropriate
        metadata that describe the object.
        '''
        print 'Writing ', self.File()
        f = h5py.File(self.File(), 'w')    
        grp = f.create_group('cenque_data')
    
        if self.data_columns is None: 
            raise ValueError
        #self.data_columns = self.get_datacol()
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
    
    def File(self): 
        ''' File name of CGPop object. The file name acts as a tool to 
        label the properties of the object.  
        ''' 
        if self.nsnap is None: 
            raise ValueError()
        if self.subhalo_prop is None: 
            raise ValueError()
     
        # Cenque Object with assigned starformation properties
        file_spec_str = self._file_spec(
                subhalo_prop=self.subhalo_prop, 
                sfr_prop=self.sfr_prop
                )
        cgpop_file = ''.join([
            'dat/CenQue/', 
            'cgpop'
            '.snapshot', str(self.nsnap), 
            file_spec_str, 
            '.hdf5'
            ]) 

        return cgpop_file

    def _file_spec(self, subhalo_prop=None, sfr_prop=None): 
        ''' file specifier depending on whether subhlao properties or sfr properties 
        have been specified 
        '''
        if subhalo_prop is None: 
            raise ValueError
    
        subhalo_str = ''.join([
            '.subhalo_', 
            'scatter', str(subhalo_prop['scatter']), '_', 
            'ancestor', str(subhalo_prop['nsnap_ancestor']), '_',
            subhalo_prop['source']
            ])

        sfr_str = ''
        if sfr_prop is not None: 
            sfr_str = ''.join([
                '.sfr_', 
                'fq_', sfr_prop['fq']['name'], '_', 
                'sfr_', sfr_prop['sfr']['name']
                ])

        return ''.join([subhalo_str, sfr_str])

    def get_datacol(self, keys=None): 
        ''' Get the list of data columns available in the class object
        '''
        # all possible data columns
        all_data = ['mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index', 
                'pos', 'gal_type', 'sfr', 'ssfr']

        datacol = [] 
        for datum in all_data: 
            if getattr(self, datum) is not None: 
                datacol.append(datum)
        return datacol

    def sample_trim(self, npwhere, quiet = False): 
        ''' Given numpy.where condition, apply numpy where condition
        to object data columns. 
        '''
        for column in self.data_columns:         
            obj_attr = getattr(self, column) 

            new_attr = obj_attr[npwhere]     

            setattr(self, column, new_attr)

        return None 
    # -------------------------
    # Galaxy properties  
    # -------------------------
    def Ssfr(self, mass=None, ssfr=None, **ssfr_kwargs):
        ''' Calculate SSFR distribution for the CGPop object

        Returns 
        -------
        ssfr_bin_mid : 
            midpoints of SSFR bins
        ssfr_hist : 
            SSFR histogram 
        '''
        if mass is None and self.mass is None: 
            raise ValueError
        if self.mass is not None: 
            mass = self.mass
            ssfr = self.ssfr
            if len(mass) != len(ssfr): 
                raise ValueError("mass, ssfr mismatch")

        gal_prop = Ssfr()
        bin_mid, ssfr_dist = gal_prop.Calculate(mass, ssfr)

        return bin_mid, ssfr_dist 

    def plotSsfr(self, ssfr_list=None, **pltkwargs):
        ''' Plot SSFR using PlotSSFR from plotting.plots 

        Parameters
        ----------
        pltkwargs : 
            - quenching : (True/False) If true plot stacked SSFR distribution 
                for CenQue data that differentiates the quiescent, quenching, 
                and star-fomring populations.
        '''
        plt.close() # in case there's another plot
        ssfr_plot = plots.PlotSSFR()
        ssfr_plot.plot(SSFRdist=self.Ssfr(), **pltkwargs)

        if 'groupcat' in pltkwargs.keys(): 
            if pltkwargs['groupcat']: 
                ssfr_plot.GroupCat(Mrcut=18)    # Mr cut is hardcoded
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
        return smf.Obj(self, **smf_kwargs)

    def plotSMF(self, **pltkwargs): 
        '''
        Plot SMF of CenQue object using PlotSMF object 
        '''
        plt.close() # in case there's another plot

        if 'line_color' in pltkwargs.keys(): 
            lc = pltkwargs['line_color']
        else: 
            lc = self.nsnap

        smf_plot = plots.PlotSMF()
        smf_plot.plot(smf=self.SMF(), line_color=self.nsnap, **pltkwargs)
        if self.zsnap < 0.1: 
            redshift = 0.1
        else: 
            redshift = self.zsnap
        smf_plot.model(redshift, source=self.subhalo_prop['source'])
        smf_plot.CentralSubhalo(self.nsnap, subhalo_prop=self.subhalo_prop)
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
        if self.sfr_prop is None: 
            raise ValueError
        
        qfrac = Fq()
        # Star-forming or Quiescent    
        return qfrac.Calculate(mass=self.mass, sfr=self.sfr, z=self.zsnap, sfms_prop=self.sfr_prop['sfms'])
    
    def plotFq(self, **pltkwargs): 
        '''
        Plot f_q of CenQue object using PlotFq plot object
        '''
        plt.close() # in case there's another plot

        FQdist = self.Fq()
        fq_plot = plots.PlotFq() 
        fq_plot.plot(FQdist=FQdist, z=self.zsnap, 
                label=r'$\mathtt{z =} '+str(self.zsnap)+'$',
                line_color=fq_plot.pretty_colors[self.nsnap],
                **pltkwargs)

        if 'model' in pltkwargs.keys():
            # hardcoded for convenience 
            if isinstance(pltkwargs['model'], bool): 
                if pltkwargs['model']: 
                    fq_plot.model(line_color='k', label=None)
            else: 
                fq_plot.model(
                        fq_prop={
                            'name': pltkwargs['model']
                            }, 
                        line_color='k', label=None) 

        fq_plot.set_axes()
        if 'show' in pltkwargs.keys(): 
            if pltkwargs['show']: 
                fq_plot.show()

        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                fq_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')
            return None
        else: 
            return fq_plot

    def plotTau(self, tau_prop, **pltkwargs): 
        ''' Plot quenching timescale as a funciton of stellar mass 
        '''
        plt.close() 
        tau_dict = tau_prop.copy()
        tau_plot = plots.PlotTau()
        tau_plot.plot(tau_dict, **pltkwargs)


        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                tau_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return tau_plot 

    def plotSMHM(self, bovyplot=True, **pltkwargs): 
        '''
        Plot M* versus M_halo relationship

        Notes
        -----
        Uses bovy_plot.scatter_plot so things are a big clunky 
        '''

        plt.close()

        smhm_plot = plots.PlotSMHM()
        smhm_plot.plot(stellarmass=self.mass, halomass=self.halo_mass, 
                bovyplot=bovyplot, **pltkwargs)
        
        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                smhm_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return smhm_plot 

    def plotSFMS(self, **pltkwargs): 
        ''' Plot the Star Forming Main Sequence of the CGPop object 
        
        Notes
        -----
        Uses bovy_plot.scatter_plot so things are a big clunky 
        '''
        plt.close() 

        sfms_plot = plots.PlotSFMS()
        if 'scatter' in pltkwargs: 
            if not pltkwargs['scatter']:
                pass 
            else: 
                sfms_plot.plot(mass=self.mass, sfr=self.sfr, sfr_class=self.sfr_class, 
                        color=self.nsnap, **pltkwargs)
        else: 
            sfms_plot.plot(mass=self.mass, sfr=self.sfr, sfr_class=self.sfr_class, 
                    color=self.nsnap, **pltkwargs)
        if 'sfqcut' in pltkwargs.keys(): 
            if pltkwargs['sfqcut']: 
                qfrac = Fq()
                m_arr = np.arange(9.0, 12.5, 0.5)
                sfms_plot.sub.plot(m_arr, qfrac.SFRcut(m_arr, self.zsnap, sfms_prop=self.sfr_prop['sfms']), c='k', ls='--', lw=4)

        if 'model' in pltkwargs.keys(): 
            if pltkwargs['model']: 
                sfms_plot.model(z=self.zsnap) 
        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                sfms_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return sfms_plot 
    

def AssignSFR(mass, redshift, sfr_prop=None, ancestor=None, descendant=None):
    ''' Assign star-forming properties based on the input mass and redshift. 
    Return sfr_class (star-forming or quiescent), SFR, and sSFR given a set 
    of masses and redshifts.

    Parameters
    ----------
    mass : array
        Array of galaxy stellar mass
    redshift : array
        Array of corresponding galaxy redshift
    sfr_prop : dictionary
        Dictionary that specifies the star-forming properties.

    Notes
    -----
    * Takes ~1.5 seconds for assign_sfr_ancestor in lineage object
    '''
    start_time = time.time()     # time the code 

    if sfr_prop is None: 
        raise ValueError('Specify SFR Properties dictionary')
    fq_dict = sfr_prop['fq']       # dictionary that specifies fq properties
    sfms_dict = sfr_prop['sfms']     # dictionary that specifies sfms properties
    if 'subhalogrowth' in sfr_prop.keys(): 
        if descendant is None: 
            raise ValueError
        if ancestor is None: 
            raise ValueError
        for key in descendant.__dict__.keys(): 
            if 'ancestor' in key: 
                ancestor_index = getattr(descendant, key) 

    ngal = len(mass) # Ngal 
    sfr_class= np.repeat('', ngal).astype('|S16') 
    sfr      = np.repeat(-999., ngal) 
    ssfr     = np.repeat(-999., ngal)
    delta_sfr   = np.repeat(-999., ngal)
    avg_sfr     = np.repeat(-999., ngal)

    # massive galaxies
    massive = np.where(mass > 0.0)[0]
    ngal_massive = len(massive)

    # f_Q(M*, z) for each massive galaxy
    qfrac = Fq()
    qf = qfrac.model(mass[massive], redshift[massive], lit=fq_dict['name']) 

    # determine quiescent/star-forming based on random's 
    # comparison to quiescent fraction
    np.random.seed()
    rand = np.random.uniform(0., 1., ngal_massive)
    quiescent = np.where(rand <= qf)
    starforming = np.where(rand > qf)
    ngal_q = len(quiescent[0])
    ngal_sf = len(starforming[0])

    # Assign SFR to quiescent galaxies
    sfr_class[massive[quiescent]] = 'quiescent'
    mu_q_ssfr = AverageLogSSFR_q_peak(mass[massive[quiescent]])
    sig_q_ssfr = ScatterLogSSFR_q_peak(mass[massive[quiescent]])
    ssfr[massive[quiescent]] = sig_q_ssfr * np.random.randn(ngal_q) + mu_q_ssfr 
    sfr[massive[quiescent]]  = ssfr[massive[quiescent]] + mass[massive[quiescent]]
    
    # Assign SFR to star-forming galaxies 
    sfr_class[massive[starforming]] = 'star-forming'
    mu_sf_sfr = AverageLogSFR_sfms(
            mass[massive[starforming]], 
            redshift[massive[starforming]], 
            sfms_prop=sfms_dict)
    sigma_sf_sfr = ScatterLogSFR_sfms(
            mass[massive[starforming]], 
            redshift[massive[starforming]],
            sfms_prop=sfms_dict)
    avg_sfr[massive[starforming]] = mu_sf_sfr
    delta_sfr[massive[starforming]] = sigma_sf_sfr * np.random.randn(ngal_sf)

    if 'subhalogrowth' in sfr_prop.keys(): 
        # if SFR assignment is dependent on subhalo growth, 
        # loop through mass bins and assign SFR by abundance matching the 
        # SFR to the Delta M*_sham. 

        mass_range = np.arange(mass[massive].min(), mass[massive].max(), 0.1)
        m_low = mass_range[:-1]
        m_high = mass_range[1:]

        for im in range(len(m_low)):
            sf_mass_bin = np.where(
                    (mass[massive] > m_low[im]) & 
                    (mass[massive] <= m_high[im]) & 
                    (sfr_class[massive] == 'star-forming')
                    )[0]
            succession, will = util.intersection_index(
                    ancestor_index, 
                    ancestor.snap_index[massive[sf_mass_bin]])
            deltaM = descendant.mass[succession] - mass[massive[sf_mass_bin[will]]]
            dM_sort_index = np.argsort(deltaM)

            deltaSFR = sigma_sf_sfr * np.random.randn(len(deltaM))
            dSFR_sort_index = np.argsort(deltaSFR)
            # abundance matched.             
            delta_sfr[massive[sf_mass_bin[will[dM_sort_index]]]] = deltaSFR[dSFR_sort_index]

    sfr[massive[starforming]] = mu_sf_sfr + delta_sfr[massive[starforming]]
    ssfr[massive[starforming]] = sfr[massive[starforming]] - mass[massive[starforming]]
    
    print 'Assign SFR took ', time.time() - start_time  

    return [sfr_class, sfr, ssfr, delta_sfr, avg_sfr] 







"""
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
