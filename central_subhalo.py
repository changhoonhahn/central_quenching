""" 


"""
import h5py
import numpy as np
import os.path

# --- Local --- 
import sham_hack as sham
from treepm import subhalo_io 
from utilities import utility as wetzel_util

class CentralSubhalos(object): 
    def __init__(self, **kwargs): 
        ''' Class that describes the Central Subhalo Catalogs generated 
        from Andrew Wetzel's TreePM merger tree code. Essentially a wrapper 
        for the central galaxies from TreePM. 
        '''
        self.kwargs = kwargs

    def File(self, **spec_kwargs): 
        ''' File name of output file given snapshot 
        '''
        # write to hdf5 file 
        subsham_cen_file = ''.join([ 
            'dat/wetzel_tree/', 
            'subhalo_sham', 
            '.central',
            '.snapshot', 
            str(spec_kwargs['snapshot']), 
            self._file_spec(**spec_kwargs), 
            '.hdf5']) 
        return subsham_cen_file 

    def _file_spec(self, **spec_kwargs): 
        ''' file specifier string that describes the key choices of parameters in the file
        '''
        spec_str = ''.join([
            '.ancestor', str(spec_kwargs['nsnap_ancestor']), 
            '.scatter', str(spec_kwargs['scatter']),
            '.', spec_kwargs['source']])
        return spec_str 

    def Read(self, snapshot, scatter=0.0, source='li-drory-march', nsnap_ancestor=20): 
        ''' Read in the hdf5 file that contains the Central Subhalo data and import 
        it appropriately to the class object structure
        '''
        spec_kwargs = {
                'scatter': scatter, 
                'snapshot': snapshot, 
                'source' : source,
                'nsnap_ancestor': nsnap_ancestor
                }
        snapshot_file = self.File(**spec_kwargs)

        self.file_name = snapshot_file
        if not os.path.isfile(snapshot_file): 
            print snapshot_file, ' does NOT exist'
            print 'Now building'
            self.build_catalogs(
                    snapshots=range(1,21), 
                    scatter=scatter, 
                    source=source, 
                    nsnap_ancestor=nsnap_ancestor
                    )

        f = h5py.File(snapshot_file, 'r')
        grp = f['cenque_data']

        for i_meta, metadatum in enumerate(grp.attrs.keys()): 
            setattr(self, metadatum, (grp.attrs.values())[i_meta])
        self.metadata = [str(key) for key in grp.attrs.keys()]

        for datum in grp.keys():
            setattr(self, datum, grp[datum][:])
        self.data_columns = [str(key).encode('ascii') for key in grp.keys()]

        return None
    
    def build_catalogs(self, snapshots=range(1,21), scatter=0.0, nsnap_ancestor=20, source='li-drory-march'): 
        ''' Build central subhalo snapshot catalogs using TreePM code. The kwargs specify the 
        properties of the central subhalos
        '''
        self.source = source
        self.scatter = scatter
        self.snapshots = snapshots
        self.nsnap_ancestor = nsnap_ancestor
        if nsnap_ancestor > np.max(snapshots): 
            raise ValueError
    
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
                source = self.source,
                zis = self.snapshots
                ) 
        #ancestor_centrals = wetzel_util.utility_catalog.indices_ilk(sub[20], ilk = 'cen') 
        for i_snap in self.snapshots: 
            central_indices = wetzel_util.utility_catalog.indices_ilk(sub[i_snap], ilk = 'cen') 
            print 'Number of Central Subhalos in Snapshot ', i_snap, '=', len(central_indices)
            print 'out of ', len(sub[i_snap]['halo.m'])
            
            print sub[i_snap]['ilk']        
            print (sub[i_snap]['ilk'])[central_indices]        
            print np.min(sub[i_snap]['ilk']), np.max(sub[i_snap]['ilk'])
            print np.min((sub[i_snap]['ilk'])[central_indices]), np.max((sub[i_snap]['ilk'])[central_indices])

            # subhalo properties
            mhalo     = (sub[i_snap]['halo.m'])[central_indices]     # M_halo
            mhalo_max = (sub[i_snap]['m.max'])[central_indices]     # M_halo
            mstar     = (sub[i_snap]['m.star'])[central_indices]     # M* of subhalo
            pos       = (sub[i_snap]['pos'])[central_indices]        # position of subhalo
            ilk       = (sub[i_snap]['ilk'])[central_indices]        # classification flag in case we want to consider centrals and virtual centrals separately 

            # keep track of child indices for subhalo
            if i_snap == np.min(self.snapshots):     
                child = np.repeat(-999, len(central_indices))
            else: 
                child = wetzel_util.utility_catalog.indices_tree(
                        sub, i_snap, i_snap-1, central_indices) 

            if i_snap == np.max(self.snapshots):     # parent index of subhalo
                parent = np.repeat(-999, len(central_indices))
            else: 
                parent = wetzel_util.utility_catalog.indices_tree(
                        sub, i_snap, i_snap+1, central_indices) 
            
            spec_kwargs = {
                    'scatter': self.scatter, 
                    'snapshot': i_snap, 
                    'source': self.source,
                    'nsnap_ancestor': nsnap_ancestor
                    }
            subsham_cen_file = self.File(**spec_kwargs)    # output name

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

            # ancestor index of subhalo for ancestor in snapshot nsnap_ancestor
            if i_snap == nsnap_ancestor:     
                grp.create_dataset(
                        'ancestor'+str(np.max(self.snapshots)), 
                        data=np.repeat(-999, len(central_indices))
                        )
            else: 
                ancestor = wetzel_util.utility_catalog.indices_tree(
                        sub, i_snap, nsnap_ancestor, central_indices)
                grp.create_dataset('ancestor'+str(nsnap_ancestor), data=ancestor)
                #for ii_snap in range(i_snap+1, np.max(self.snapshots)+1): 
                #    print i_snap, ii_snap
                #    ancestor = wetzel_util.utility_catalog.indices_tree(
                #            sub, i_snap, ii_snap, central_indices)
                #    grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
            f.close() 

        return None


class Subhalos(CentralSubhalos): 
    def __init__(self, **kwargs): 
        ''' Class that describes the Subhalo Catalogs generated from 
        Andrew Wetzel's TreePM merger tree code. This class describes
        all subhalos (both centrals and satellites) 

        '''
        super(Subhalos, self).__init__(**kwargs) 

    def File(self, **spec_kwargs): 
        ''' File name of output file given snapshot 
        '''
        # write to hdf5 file 
        subsham_file = ''.join([ 
            'dat/wetzel_tree/', 
            'subhalo_sham', 
            '.all',
            '.snapshot', 
            str(spec_kwargs['snapshot']), 
            self._file_spec(**spec_kwargs), 
            '.hdf5']) 
        return subsham_file 

    def build_catalogs(self, snapshots=range(1,21), scatter = 0.0, source='li-drory-march', nsnap_ancestor=20): 
        '''
        Build subhalo snapshot catalogs
        '''
        self.source = source
        self.scatter = scatter
        self.snapshots = snapshots
        self.nsnap_ancestor = nsnap_ancestor

        if nsnap_ancestor > np.max(snapshots): 
            raise ValueError
    
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
                source = self.source, 
                zis = self.snapshots 
                ) 

        for i_snap in self.snapshots: 
            
            print 'Number of Subhalos in Snapshot ', i_snap, ' =', len(sub[i_snap]['halo.m'])

            # subhalo properties
            mhalo     = (sub[i_snap]['halo.m'])     # M_halo
            mhalo_max = (sub[i_snap]['m.max'])     # M_halo
            mstar     = (sub[i_snap]['m.star'])     # M* of subhalo
            pos       = (sub[i_snap]['pos'])        # position of subhalo
            ilk       = (sub[i_snap]['ilk'])        # classification flag in case we want to consider centrals and virtual centrals separately 

            # keep track of child indices for subhalo
            if i_snap == np.min(self.snapshots):     
                child = np.repeat(-999, len(mhalo))
            else: 
                child = wetzel_util.utility_catalog.indices_tree(sub, i_snap, i_snap-1, range(len(mhalo)))
                print len(mhalo), len(child)

            if i_snap == np.max(self.snapshots):     # parent index of subhalo
                parent = np.repeat(-999, len(mhalo))
            else: 
                parent = wetzel_util.utility_catalog.indices_tree(sub, i_snap, i_snap+1, range(len(mhalo)))
                print len(parent), len(parent)
            
            spec_kwargs = {
                    'scatter': self.scatter, 
                    'snapshot': i_snap, 
                    'source': self.source,
                    'nsnap_ancestor': nsnap_ancestor
                    }
            subsham_file = self.File(**spec_kwargs)    # output name
            f = h5py.File(subsham_file, 'w') 
            
            grp = f.create_group('cenque_data') 
            
            grp.attrs['snapshot'] = i_snap
            grp.create_dataset('index', data=range(len(mhalo)))
            grp.create_dataset('pos', data=pos) 
            grp.create_dataset('ilk', data=ilk)
            grp.create_dataset('mass', data=mstar)
            grp.create_dataset('child', data=child)
            grp.create_dataset('parent', data=parent)
            grp.create_dataset('halo.m', data=mhalo)
            grp.create_dataset('halo.m.max', data=mhalo_max)

            if i_snap == nsnap_ancestor:     # ancestor index of subhalo
                ancestor = np.repeat(-999, len(mhalo))
            else: 
                ii_snap = nsnap_ancestor
                ancestor = wetzel_util.utility_catalog.indices_tree(
                        sub, i_snap, ii_snap, range(len(mhalo)))
                grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
                #for ii_snap in range(i_snap+1, np.max(self.snapshots)+1): 
                #    ancestor = wetzel_util.utility_catalog.indices_tree(
                #            sub, i_snap, ii_snap, range(len(mhalo)))
                #    grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
            f.close() 

        return None




if __name__=='__main__': 
    for nsnap_ancestor in [10, 15, 20]: 
        for scat in [ 0.0, 0.2 ]: 
            subh = Subhalos() 
            subh.build_catalogs(
                    snapshots=range(1, nsnap_ancestor+1), 
                    scatter=scat, 
                    source='li-march', 
                    nsnap_ancestor=nsnap_ancestor)
            del subh
