'''

Central_Quenchign project


Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import h5py
import time

#---- Local ----
import cenque_utility as util
import sf_mainseq as sfms

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
        self.nsnap = None       # n_snapshot 
        self.zsnap = None           # z_snapshot
        self.t_cosmic = None    # t_cosmic for snapshot
        self.t_step = None      # t_cosmic step 

        self.mass_bins = self.set_mass_bins()
        self.data_columns = None 
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

        self.data_columns = ['mass', 'parent', 'child', 'ilk', 'snap_index', 'pos']
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

        massbins = util.simple_mass_bin() 
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
        # kwargs specifies the input file 
        input_file = util.cenque_file( **kwargs ) 
        f = h5py.File(input_file, 'r') 
        print input_file 

        grp = f['cenque_data']

        # save meta data first
        for i_meta, metadatum in enumerate(grp.attrs.keys()): 
            setattr(self, metadatum, (grp.attrs.values())[i_meta]) 
    
        for i_col, column in enumerate(grp.keys()): 
            setattr(self, column, grp[column][:])
        
        if 'sham_mass' not in grp.keys(): 
            setattr(self, 'sham_mass', grp['mass'][:])

        f.close() 

        return None 

    def writeout(self): 
        """ Write out self.data_columns to h5py file
        """
        
        output_file = self.file()  # output file 

        print 'Writing ', output_file 
        f = h5py.File(output_file, 'w')         # hdf5 file format (open file) 
        grp = f.create_group('cenque_data')     # create group 
    
        # set output columns 
        if self.sfr is None: 
            # basic 
            columns = ['mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index']
        elif 'columns' in kwargs.keys():
            # if columns are specified 
            columns = kwargs['columns']
        else:       
            # if SFR/SSFR have been assigned
            if kwargs['sfr'] == 'sfr_func': 
                columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'halo_mass', 
                        'sfr_phase', 'sfr_freq', 'sfr_amp', 
                        'parent', 'child', 'ilk', 'snap_index']
            elif kwargs['sfr'] == 'sfr_avg': 
                columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'halo_mass', 
                        'sfr_resid', 'parent', 'child', 'ilk', 'snap_index']
            else: 
                raise NotImplementedError('kasdflk')

        n_cols = len(columns)       # number of columns 
        col_fmt = []                # column format 
    
        # save each column to dataset within 'cenque_data' group 
        for column in columns:      
            column_attr = getattr(self, column)

            # write out file
            grp.create_dataset( column, data=column_attr ) 
        
        # save metadata 
        metadata = [ 'nsnap', 'zsnap', 't_cosmic', 't_step' ]
        for metadatum in metadata: 
            if self.nsnap is not None: 
                grp.attrs[metadatum] = getattr(self, metadatum) 

        f.close()  

        return None 
    
    def file(self): 
        """ File name of CenQue object. Only used to writeout 
        CenQue class object
        """

        if self.nsnap == None: 
            raise ValueError()

        if self.cenque_type == None: 
            raise ValueError()
    
        # Cenque Object with assigned starformation properties
        if self.cenque_type == 'sfassigned':

            #if kwargs['sfr'] == 'sfr_avg': 
            #    file_type_str = '_sfravg_sfpropassign'
            #elif kwargs['sfr'] == 'sfr_func': 
            #    file_type_str = '_sfrfunc_sfpropassign'
            #else: 
            #    raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 
                
            # Quenching Fraction specifier 
            if 'fqing_slope' in kwargs.keys(): 
                fqing_slope_str = str("%.2f" % kwargs['fqing_slope'])
            else: 
                fqing_slope_str = str("%.2f" % 0.63)

            if 'fqing_yint' in kwargs.keys(): 
                fqing_yint_str = str("%.2f" % kwargs['fqing_yint'])
            else: 
                fqing_yint_str = str("%.2f" % -6.04) 

            fqing_str = ''.join([fqing_slope_str, '_', fqing_yint_str, 'fqing']) 

            # combine specifiers
            file_type_str = ''.join(['_', fqing_str, file_type_str])
                
        elif 'evo_from' in self.cenque_type:
            # evolved from nsnap
            original_nsnap = int(((kwargs['file_type']).split('from'))[-1]) 
        
            if kwargs['stellmass'].lower() == 'integrated': 
                mass_str = '_integ'
            elif kwargs['stellmass'].lower() == 'sham': 
                mass_str = '_sham'
            else: 
                raise NotImplementedError('asdfalkjlkjasdf') 

            # star-formation assign (random about average or functional) 
            if kwargs['sfr'] == 'sfr_avg': 
                file_type_str = mass_str+'_sfravg_evol_from'+str(original_nsnap) 
            elif kwargs['sfr'] == 'sfr_func': 
                file_type_str = mass_str+'_sfrfunc_evol_from'+str(original_nsnap) 
            else: 
                raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 
           
            # Tau specifier
            if kwargs['tau'] == 'discrete': 
                tau_str = '_'+'_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
            elif kwargs['tau'] == 'linefit':
                tau_str = '_line'+'_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
            else: 
                tau_str = '_'+kwargs['tau']+'tau'

            # Quenching Fraction specifier 
            if 'fqing_slope' in kwargs.keys(): 
                fqing_slope_str = str(kwargs['fqing_slope'])
            else: 
                fqing_slope_str = str(0.63)

            if 'fqing_yint' in kwargs.keys(): 
                fqing_yint_str = str(kwargs['fqing_yint'])
            else: 
                fqing_yint_str = str(-6.04) 

            fqing_str = '_'+fqing_slope_str+'_'+fqing_yint_str+'fqing'

            # combine specifiers
            file_type_str = ''.join([tau_str, fqing_str, file_type_str]) 

        else: 
            raise NameError() 

        if 'sfms_slope' not in kwargs.keys(): 
            sfms_param_str = ''
        else: 
            slope_str = "%.2f" % kwargs['sfms_slope']
            yint_str = "%.2f" % kwargs['sfms_yint'] 
            sfms_param_str = '_sfms_slope'+slope_str+'_yint'+yint_str

        cenque_filename = ''.join(['dat/central_quenching/', 
            'cenque_centrals_snapshot', str(kwargs['nsnap']), file_type_str, sfms_param_str, 
            '.hdf5']) 

    return cenque_filename

    def sample_trim(self, npwhere, columns = self.data_columns, quiet = False): 
        ''' Given numpy.where condition, apply numpy where condition
        to object data columns. 
        '''

        if not quiet: 
            #n_remove = len(bool) - len(bool[bool == True])
            #print 'Removing ', n_remove, ' elements from ', len(bool)  
    
        if columns is None:         # if data columns aren't specified
            data_columns = ['mass', 'halo_mass', 'sfr', 'ssfr', 'gal_type', 
                    'parent', 'child', 'ilk', 'snap_index', 'pos']
        else: 
            data_columns = columns 

        n_list = len(bool) 
        
        for column in data_columns:         # loop through columns and only keep true values
            attr_list = getattr(self, column) 

            try: 
                n_attr_list = len(attr_list) 
            except TypeError: 
                n_attr_list = 0 
            
            if n_list != n_attr_list: 
                raise TypeError(column+": boolean does not match length!") 
            else:  
                # impose boolean so only "true" values are kept
                new_attr_list = attr_list[bool]     

                setattr(self, column, new_attr_list)    # save to class 
    

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
    #EvolveCenQue(13, 1, fqing_yint=-5.84, tau='instant')  
    #tau='linefit', tau_param=[-0.5, 0.4]) 
    #EvolveCenQue(13, 1, fqing_yint=-5.84, tau='linefit', tau_param=[-0.4, 0.2])
    build_cenque_importsnap() 
    build_cenque_original(sfr='sfr_func') 
    EvolveCenQue(13, 1, tau='linefit', tau_param=[-0.7, 0.4], 
            sfr='sfr_func', stellmass='sham') 
    EvolveCenQue(13, 1, tau='linefit', tau_param=[-0.7, 0.4], 
            sfr='sfr_func', stellmass='integrated') 
