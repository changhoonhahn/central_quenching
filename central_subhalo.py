"""



"""
import h5py
import numpy as np

# --- Local --- 
import pyfits
from treepm import sham
from treepm import subhalo_io 
from utilities import utility as wetzel_util

class CentralSubhalos: 

    def __init__(self, scatter = 0.0): 
        '''
        Class that describes the Central Subhalo Catalogs generated from Andrew Wetzel's TreePM merger tree code 

        '''
        self.scatter = scatter

    def file(self, snapshot): 
        ''' File name of output 
        '''
        # write to hdf5 file 
        subsham_cen_file = ''.join([ 
            'dat/wetzel_tree/', 
            'subhalo_sham_centrals_snapshot', 
            str(snapshot), 
            '_scatter', 
            str(self.scatter), 
            '.hdf5'
            ]) 

        return subsham_cen_file 

    def read(self, snapshot): 
        '''

        '''
        snapshot_file = self.file(snapshot)

        f = h5py.File(snapshot_file, 'r')
        grp = f['cenque_data']

        for i_meta, metadatum in enumerate(grp.attrs.keys()): 
            setattr(self, metadatum, (grp.attrs.values())[i_meta])

        for datum in grp.keys():
            setattr(self, datum, grp[datum][:])

        return None
    
    def build_catalogs(self, snapshots = range(1,21)): 
        '''
        Build catalogs
        '''
        self.snapshots = snapshots
    
        # read in TreePM Subhalo snapshots from z~0.0502 to z~1.0833
        # note: 250 here is the simulation box length (Mpc/h comoving)
        sub = subhalo_io.Treepm.read(
                'subhalo', 
                250, 
                zis = self.snapshots
                ) 

        # assign M* to subhalos using SHAM
        sham.assign(
                sub, 
                'm.star', 
                scat = self.scatter, 
                dis_mf = 0.0, 
                zis = self.snapshots
                ) 

        for i_snap in self.snapshots: 
            central_indices = wetzel_util.utility_catalog.indices_ilk(sub[i_snap], ilk = 'cen') 
            
            print 'Number of Central Subhalos in Snapshot ', i_snap, '=', len(central_indices) 
    
            # subhalo properties
            mhalo     = (sub[i_snap]['halo.m'])[central_indices]     # M_halo
            mhalo_max = (sub[i_snap]['m.max'])[central_indices]     # M_halo
            mstar     = (sub[i_snap]['m.star'])[central_indices]     # M* of subhalo
            pos       = (sub[i_snap]['pos'])[central_indices]        # position of subhalo
            ilk       = (sub[i_snap]['ilk'])[central_indices]        # classification flag in case we want to consider centrals and virtual centrals separately 

            # keep track of child indices for subhalo
            if i_snap == np.min(self.snapshots):     
                child = np.zeros(len(central_indices))
                child = child - 999.
                child.astype(int)
            else: 
                child = wetzel_util.utility_catalog.indices_tree(sub, i_snap, i_snap-1, central_indices) 

            if i_snap == np.max(self.snapshots):     # parent index of subhalo
                parent = np.zeros(len(central_indices))
                parent = parent - 999.
                parent.astype(int)
            else: 
                parent = wetzel_util.utility_catalog.indices_tree(sub, i_snap, i_snap+1, central_indices) 
            
            subsham_cen_file = self.file(i_snap)    # output name
            f = h5py.File(subsham_cen_file, 'w') 

            grp = f.create_group('cenque_data') 
            
            grp.attrs['snapshot'] = i_snap
            grp.create_dataset('index', data=central_indices)
            grp.create_dataset('pos', data=pos) 
            grp.create_dataset('ilk', data=ilk)
            grp.create_dataset('mass', data=mstar)
            grp.create_dataset('child', data=child)
            grp.create_dataset('parent', data=parent)
            grp.create_dataset('halo.m', data=mhalo)
            grp.create_dataset('halo.m.max', data=mhalo_max)
            
            print 'M*: min ', np.min(mstar), ' max ', np.max(mstar), ' mean ', np.mean(mstar)
            print 'Snapshot ', i_snap,' has ', len(child[child >= 0]), ' parents that have children' 
            print 'Snapshot ', i_snap,' has ', len(parent[parent >= 0]), ' children that have parents'

            f.close() 

        return None

if __name__=='__main__': 
    censub = CentralSubhalos(scatter = 0.2)
    censub.build_catalogs()
    censub.read(10)
    print censub.__dict__.keys()
    print censub.snapshot
    print getattr(censub, 'halo.m.max')
