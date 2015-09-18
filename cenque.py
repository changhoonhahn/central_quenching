'''

Central_Quenchign project


Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import h5py
import time

#---- Local ----
from assign_sfr import assign_sfr 
from util.mass_bins import simple_mass_bin

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
        self.zsnap = None           # z_snapshot
        self.t_cosmic = None    # t_cosmic for snapshot
        self.t_step = None      # t_cosmic step 

        self.mass_bins = self.set_mass_bins()
        self.data_columns = None 
        if 'cenque_type' in self.kwargs.keys(): 
            self.cenque_type = self.kwargs['cenque_type']
        else: 
            self.cenque_type = None

    def import_treepm(self, nsnap): 
        """ Import the following snapshot data from TreePM --> SHAM snapshots. Also imports 
        snapshot, redshift and t_cosmic metadata. 
        * SHAM stellar mass 
        * parent index 
        * child index 
        * halo mass 
        """
        start_time = time.time()

        snapshot_file = ''.join([
            'dat/wetzel_tree/', 
            'subhalo_sham_centrals_snapshot', str(nsnap), '.hdf5'
            ]) 

        f = h5py.File(snapshot_file, 'r')       # read in h5py file 
        grp = f['cenque_data'] 

        self.data_columns = ['mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index', 'pos']
        self.cenque_type = 'treepm_import'

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

        # get snapshot redshift/cosmic time data using Andrew's table
        n_snaps, z_snap, t_snap, t_wid = np.loadtxt('snapshot_table.dat', 
                unpack=True, usecols=[0, 2, 3, 4])

        self.nsnap = nsnap 
        self.zsnap = z_snap[(n_snaps.tolist()).index(nsnap)]        # redshift of snapshot
        self.t_cosmic = t_snap[(n_snaps.tolist()).index(nsnap)] # t_cosmic of snapshot 
        self.t_step = t_wid[(n_snaps.tolist()).index(nsnap)]    # Gyrs until next snapshot

        print "import_treepm takes ", time.time() - start_time

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
        f = h5py.File(cq_file, 'r') 
        print 'Reading ', cq_file 

        grp = f['cenque_data']

        # save meta data first
        for i_meta, metadatum in enumerate(grp.attrs.keys()): 
            setattr(self, metadatum, (grp.attrs.values())[i_meta]) 

        for i_col, column in enumerate(grp.keys()): 
            setattr(self, column, grp[column][:])
        
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
        metadata = [ 'nsnap', 'zsnap', 't_cosmic', 't_step' ]
        for metadatum in metadata: 
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
    
        # Cenque Object with assigned starformation properties
        if self.cenque_type == 'treepm_import':
            file_type_str = ''

        elif self.cenque_type == 'sf_assigned':

            # default properties
            if 'sf_prop' not in self.__dict__.keys():
                self.sf_prop = {'name': 'average'}
            if 'fq_prop' not in self.__dict__.keys():
                self.fq_prop = {'name': 'wetzelsmooth'}
            if 'mass_evol' not in self.__dict__.keys():
                self.mass_evol = 'sham'

            sfr_str = '_'
            if self.sf_prop['name'] == 'average': 
                sfr_str += 'average'
            else: 
                raise NotImplementedError()
            sfr_str += '_sfassign'
                
            fq_str = '_'
            if self.fq_prop['name'] == 'wetzelsmooth': 
                fq_str += 'wetzelsmooth' 
            fq_str +='_fq'

            # combine specifiers
            file_type_str = ''.join([fq_str, sfr_str])
                
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
            else: 
                raise NotImplementedError()
            sfr_str += '_sfassign'
                
            fq_str = '_'
            if self.fq_prop['name'] == 'wetzelsmooth': 
                fq_str += 'wetzelsmooth' 
            fq_str +='_fq'

            # Tau specifier
            if self.tau_prop['name'] in ('instant', 'constant'): 
                tau_str = ''.join(['_', self.tau_prop['name'], 'tau'])

            # combine specifiers
            file_type_str = ''.join([tau_str, fq_str, sfr_str, '_', self.cenque_type]) 

        else: 
            raise NameError() 

        cenque_filename = ''.join([
            'dat/central_quenching/', 
            'cenque_centrals_snapshot', 
            str(self.nsnap), 
            file_type_str, 
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
    
def build_cenque_importsnap(): 
    ''' Build CenQue snapshots with TreePM data imported
    
    Notes
    -----
    * pretty much hardcoded
    '''
    for i_snap in np.arange(1, 13, 1): 
        snap = CenQue() 
        snap.ImportSnap(nsnap=i_snap)
        snap.writeout(nsnap=i_snap)

def build_cenque_original(i_snap=13, **kwargs):
    snap = CenQue() 
    snap.AssignSFR(i_snap, **kwargs) 
    snap.writeout(nsnap=i_snap, file_type='sf assign', **kwargs)

if __name__=='__main__': 
    blah = CenQue()
    blah.import_treepm(13)
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
