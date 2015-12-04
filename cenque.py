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

#---- Local ----
from ssfr import Ssfr
from assign_sfr import assign_sfr 
from util.mass_bins import simple_mass_bin
from central_subhalo import CentralSubhalos
from quiescent_fraction import sfq_classify
# -- plotting --
from plotting.plot_fq import PlotFq
from plotting.plot_tau import PlotTau
from plotting.plot_ssfr import PlotSSFR
from plotting.plot_sfms import PlotSFMS
from plotting.plot_mstar_mhalo import PlotMstarMhalo

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

        if 'scatter' in self.kwargs.keys(): 
            self.mass_scatter = self.kwargs['scatter']
        else: 
            self.mass_scatter = None       # n_snapshot 

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

    def import_treepm(self, nsnap, scatter=0.0): 
        """ Import the following snapshot data from TreePM --> SHAM snapshots. Also imports 
        snapshot, redshift and t_cosmic metadata. 
        * SHAM stellar mass 
        * parent index 
        * child index 
        * halo mass 
        """
        #start_time = time.time()
        
        centsub = CentralSubhalos(scatter=scatter)
        snapshot_file = centsub.file(nsnap)

        f = h5py.File(snapshot_file, 'r')       # read in h5py file 
        grp = f['cenque_data'] 

        self.data_columns = ['mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index', 'pos']
        
        self.mass = grp['mass'][:] 
        self.parent = grp['parent'][:]
        self.child = grp['child'][:]
        self.ilk = grp['ilk'][:]
        self.snap_index = grp['index'][:]
        self.pos = grp['pos'][:]
        #self.halo_mass = grp['mass_halo'][:]

        # maximum halo mass from halo merger tree.
        # this halo mass is used for SHAM rather
        self.halo_mass = grp['halo.m.max'][:]       
    
        # Meta data of snapshot 
        self.metadata = [ 'nsnap', 'zsnap', 't_cosmic', 't_step', 'cenque_type', 'mass_scatter']
        # get snapshot redshift/cosmic time data using Andrew's table
        n_snaps, z_snap, t_snap, t_wid = np.loadtxt('snapshot_table.dat', 
                unpack=True, usecols=[0, 2, 3, 4])
    
        self.nsnap = nsnap 
        self.zsnap = z_snap[(n_snaps.tolist()).index(nsnap)]        # redshift of snapshot
        self.t_cosmic = t_snap[(n_snaps.tolist()).index(nsnap)] # t_cosmic of snapshot 
        self.t_step = t_wid[(n_snaps.tolist()).index(nsnap)]    # Gyrs until next snapshot
        self.cenque_type = 'treepm_import'
        self.mass_scatter = scatter

        #print "import_treepm takes ", time.time() - start_time

        return None

    def set_mass_bins(self): 
        """ Mass bin of 
        """

        massbins = simple_mass_bin() 
        # nothing else specified here but this should later be implemented 
        # to change based on kwargs

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
        """ File name of CenQue object. Only used when writeout method 
        is called
        """
        
        if self.nsnap == None: 
            raise ValueError()

        if self.cenque_type == None: 
            raise ValueError()

        if self.mass_scatter == None: 
            raise ValueError()
    
        # Cenque Object with assigned starformation properties
        if self.cenque_type == 'treepm_import':
            file_type_str = ''

        elif self.cenque_type == 'sf_assigned':

            # default properties
            if 'sf_prop' not in self.__dict__.keys():
                self.sf_prop = {'name': 'average'}
            if 'fq_prop' not in self.__dict__.keys():
                self.fq_prop = {'name': 'wetzelsmooth'}
            if 'tau_prop' not in self.__dict__.keys():
                self.tau_prop = {'name': 'instant'}
            if 'mass_evol' not in self.__dict__.keys():
                self.mass_evol = 'sham'

            sfr_str = '_'
            if self.sf_prop['name'] == 'average': 
                sfr_str += 'average'
            elif self.sf_prop['name'] == 'average_noscatter': 
                sfr_str += 'average_noscatter'
            else: 
                raise NotImplementedError()
            sfr_str += '_sfassign'
                
            fq_str = '_'
            if self.fq_prop['name'] == 'wetzelsmooth': 
                fq_str += 'wetzelsmooth' 
            fq_str +='_fq'

            #if self.tau_prop['name'] in ('instant', 'constant', 'satellite', 'long'): 
            #    tau_str = ''.join(['_', self.tau_prop['name'], 'tau'])
            #elif self.tau_prop['name'] in ('line'):
            #    tau_str = ''.join([
            #        '_', self.tau_prop['name'], 'tau', 
            #        '_Mfid', str(self.tau_prop['fid_mass']), 
            #        '_slope', str(round(self.tau_prop['slope'], 4)), 
            #        '_yint', str(round(self.tau_prop['yint'],4))
            #        ])

            # combine specifiers
            file_type_str = ''.join([fq_str, sfr_str])
            #file_type_str = ''.join([tau_str, fq_str, sfr_str])
                
        elif 'evol_from' in self.cenque_type:
            # Cenque Object evolved from some n_snap 

            original_nsnap = int(
                    (self.cenque_type.split('from'))[-1]
                    ) 
            
            # default properties
            if 'sf_prop' not in self.__dict__.keys():
                self.sf_prop = {'name': 'average'}
            if 'fq_prop' not in self.__dict__.keys():
                self.fq_prop = {'name': 'wetzelsmooth'}
            if 'tau_prop' not in self.__dict__.keys():
                self.tau_prop = {'name': 'instant'}
            if 'mass_evol' not in self.__dict__.keys():
                self.mass_evol = 'sham'

        
            if self.mass_evol == 'integrated': 
                mass_str = '_integ'
            elif self.mass_evol == 'sham': 
                mass_str = '_sham'
            else: 
                raise NotImplementedError('asdfalkjlkjasdf') 
            
            sfr_str = '_'
            if self.sf_prop['name'] == 'average': 
                sfr_str += 'average'
            elif self.sf_prop['name'] == 'average_noscatter': 
                sfr_str += 'average_noscatter'
            else: 
                raise NotImplementedError()
            sfr_str += '_sfassign'
                
            fq_str = '_'
            if self.fq_prop['name'] == 'wetzelsmooth': 
                fq_str += 'wetzelsmooth' 
            fq_str +='_fq'

            # Tau specifier
            if self.tau_prop['name'] in ('instant', 'constant', 'satellite', 'long'): 
                tau_str = ''.join(['_', self.tau_prop['name'], 'tau'])
            elif self.tau_prop['name'] in ('line'):
                tau_str = ''.join([
                    '_', self.tau_prop['name'], 'tau', 
                    '_Mfid', str(self.tau_prop['fid_mass']), 
                    '_slope', str(round(self.tau_prop['slope'], 4)), 
                    '_yint', str(round(self.tau_prop['yint'],4))
                    ])

            # combine specifiers
            file_type_str = ''.join([tau_str, fq_str, sfr_str, '_', self.cenque_type]) 

        else: 
            raise NameError() 

        cenque_filename = ''.join([
            'dat/central_quenching/', 
            'cenque_centrals_snapshot', 
            str(self.nsnap), 
            file_type_str, 
            '_mass_scatter',
            str(round(self.mass_scatter,1)),
            '.hdf5'
            ]) 

        return cenque_filename

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

        ssfr_plot = PlotSSFR()
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

    # --- Quiescent Fraction ---
    def Fq(self, **fq_kwargs):
        '''
        Calculate quiescent fraction of CenQue class object
        '''
        if cenque.zsnap is None: 
            raise ValueError
        if cenque.mass is None: 
            raise ValueError
        if cenque.sfr is None: 
            raise ValueError

        # Star-forming or Quiescent    
        sf_q = sfq_classify(self.mass, self.sfr, self.zsnap, Mrcut=18)

        masses, f_q = [], [] 
        for i_bin in xrange(self.mass_bins.nbins): 

            in_massbin = np.where( 
                    (cenque.mass > self.mass_bins.mass_low[i_bin]) &
                    (cenque.mass <= self.mass_bins.mass_high[i_bin])
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

        fq_plot = PlotFq() 
        fq_plot.cenque(self)

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

    def plotTau(self, **pltkwargs): 
        '''
        Plot quenching timescale 
        '''

        plt.close() 

        tau_plot = PlotTau(self.tau_prop)

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

        mm_plot = PlotMstarMhalo(self, **pltkwargs)
        
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

        sfms_plot = PlotSFMS
        sfms_plot.cenque(self, **pltkwargs)
        
        if 'savefig' in pltkwargs.keys():
            if isinstance(pltkwargs['savefig'], str): 
                sfms_plot.save_fig(pltkwargs['savefig'])
            else: 
                ValueError('savefig = figure_file_name')

            return None

        else: 
            return mm_plot 
        
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

"""
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
