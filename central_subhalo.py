"""



"""
import h5py
import numpy as np
import os.path

# --- Local --- 
from smf import SMF
#from treepm import sham
import sham_hack as sham
from treepm import subhalo_io 
from util.cenque_utility import get_z_nsnap
from util.cenque_utility import intersection_index
from utilities import utility as wetzel_util

from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

class CentralSubhalos(object): 
    def __init__(self, **kwargs): 
        '''
        Class that describes the Central Subhalo Catalogs generated 
        from Andrew Wetzel's TreePM merger tree code.

        '''
        self.kwargs = kwargs

    def file(self, **spec_kwargs): 
        ''' File name of output file given snapshot 
        '''
        self.spec_str = self.file_spec(**spec_kwargs)
        # write to hdf5 file 
        subsham_cen_file = ''.join([ 
            'dat/wetzel_tree/', 'subhalo_sham.centrals.snapshot', 
            str(spec_kwargs['snapshot']), self.spec_str, '.hdf5'
            ]) 

        return subsham_cen_file 

    def file_spec(self, **spec_kwargs): 
        '''
        file specifier string that describes the key choices of parameters in the file
        '''
        spec_str = ''.join([
            '.scatter', str(spec_kwargs['scatter']),
            '.source_', spec_kwargs['source']
            ])
        return spec_str 

    def read(self, snapshot, scatter = 0.0, source='li-drory-march' ): 
        '''
        '''
        spec_kwargs = {
                'scatter': scatter, 
                'snapshot': snapshot, 
                'source' : source
                }
        snapshot_file = self.file(**spec_kwargs)
        self.file_name = snapshot_file
        if not os.path.isfile(snapshot_file): 
            print snapshot_file, ' does NOT exist'
            print 'Now building'
            self.build_catalogs(snapshots=range(1,21), scatter = scatter, source=source)

        f = h5py.File(snapshot_file, 'r')
        grp = f['cenque_data']

        for i_meta, metadatum in enumerate(grp.attrs.keys()): 
            setattr(self, metadatum, (grp.attrs.values())[i_meta])
        self.metadata = [str(key) for key in grp.attrs.keys()]

        for datum in grp.keys():
            setattr(self, datum, grp[datum][:])
        self.data_columns = [str(key) for key in grp.keys()]

        return None
    
    def build_catalogs(self, snapshots = range(1,21), scatter = 0.0, source='li-drory-march'): 
        '''
        Build subhalo snapshot catalogs
        '''
        self.source = source
        self.scatter = scatter
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
                source = self.source,
                zis = self.snapshots
                ) 
        ancestor_centrals = wetzel_util.utility_catalog.indices_ilk(sub[20], ilk = 'cen') 
        for i_snap in self.snapshots: 
            central_indices = wetzel_util.utility_catalog.indices_ilk(sub[i_snap], ilk = 'cen') 
            print 'Number of Central Subhalos in Snapshot ', i_snap, '=', len(central_indices)
            print 'out of ', len(sub[i_snap]['halo.m'])
    
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
            
            if i_snap == np.max(self.snapshots):     # ancestor index of subhalo
                ancestor = np.repeat(-999, len(central_indices))
            else: 
                ancestor = wetzel_util.utility_catalog.indices_tree(
                        sub, i_snap, np.max(self.snapshots), central_indices) 
            
            spec_kwargs = {
                    'scatter': self.scatter, 
                    'snapshot': i_snap, 
                    'source' : self.source
                    }
            subsham_cen_file = self.file(**spec_kwargs)    # output name

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
            grp.create_dataset('ancestor'+str(np.max(self.snapshots)), data=ancestor)
            
            f.close() 

        return None


class Subhalos(CentralSubhalos): 
    def __init__(self, scatter = 0.0): 
        '''
        Class that describes the Subhalo Catalogs generated from 
        Andrew Wetzel's TreePM merger tree code. This class describes
        both central and satellite catalogs

        '''
        super(Subhalos, self).__init__(scatter = scatter) 

    def file(self, **spec_kwargs): 
        ''' File name of output 
        '''
        self.spec_str = self.file_spec(**spec_kwargs)
        # write to hdf5 file 
        subsham_cen_file = ''.join([ 
            'dat/wetzel_tree/', 'subhalo_sham.all.snapshot', 
            str(spec_kwargs['snapshot']), self.spec_str, '.hdf5'
            ]) 

        return subsham_cen_file 

    def build_catalogs(self, snapshots = range(1,21), scatter = 0.0, source='li-drory-march'): 
        '''
        Build subhalo snapshot catalogs
        '''
        self.source = source
        self.scatter = scatter
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
                    'source' : self.source
                    }
            subsham_file = self.file(**spec_kwargs)    # output name
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

            if i_snap == np.max(self.snapshots):     # ancestor index of subhalo
                ancestor = np.repeat(-999, len(mhalo))
            else: 
                for ii_snap in range(i_snap+1, np.max(self.snapshots)+1): 
                    ancestor = wetzel_util.utility_catalog.indices_tree(
                            sub, i_snap, ii_snap, range(len(mhalo)))
                    grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
            f.close() 

        return None


def test_subhalos_smf(scatter = 0.0, source='li-drory-march', type='all'): 
    '''
    Test subhalos by comparing their SMFs to analytic SMFs 

    test_subhalos_smf(scatter = 0.0, source='li-march', type='all')
    test_subhalos_smf(scatter = 0.0, source='li-march', type='central')
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    for i_snap in range(1, 21):
        smf = SMF()
        if type == 'all': 
            subh = Subhalos()
        elif type == 'central': 
            subh = CentralSubhalos()
        else: 
            raise ValueError

        subh.read(i_snap, scatter=scatter, source=source)
        
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        sub.plot(subh_mass, subh_phi, lw=4, c='gray') 

        analytic_mass, analytic_phi = smf.analytic(get_z_nsnap(i_snap), source=source) 
        sub.plot(analytic_mass, analytic_phi, lw=4, ls='--', c='k', 
                label=r"$ z = "+str(round(get_z_nsnap(i_snap),2))+"$") 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    #sub.legend(loc='upper right')
    if type == 'all': 
        subhalo_str = 'subhalo' 
    elif type == 'central': 
        subhalo_str = 'central_subhalo' 
    else: 
        raise ValueError
    fig.savefig(
            ''.join([
                'figure/'
                'qaplot_', subhalo_str, '_', source, '.png'
                ]), 
            bbox_inches='tight'
            )
    plt.close()


def test_orphan_subhalo_smf(nsnap_ancestor = 20, scatter = 0.0, source='li-drory-march'): 
    '''
    SMF for orphan subhalos
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    subh = CentralSubhalos()
    subh.read(nsnap_ancestor, scatter=scatter, source=source)
    ancestor_index = subh.index
    ancestor_mass = subh.mass

    for i_snap in [1,5,10, 15]:
        smf = SMF()
        subh = CentralSubhalos()
        subh.read(i_snap, scatter=scatter, source=source)
    
        # total central subhalo SMF
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        sub.plot(subh_mass, subh_phi, lw=4, c=pretty_colors[i_snap], alpha=0.5) 

        has_ancestor = np.where(subh.ancestor20 > 0)[0]

        # SMF of central subhalos with ancestors
        nonorph_mass, nonorph_phi = smf.smf(subh.mass[has_ancestor])
        sub.plot(nonorph_mass, nonorph_phi, lw=2, ls='--', c=pretty_colors[i_snap]) 

    anc_mass, anc_phi = smf.smf(ancestor_mass)
    sub.plot(anc_mass, anc_phi,  lw=4, c='gray') 
    #anc_mass, anc_phi = smf.smf(ancestor_mass[int_ancestors])
    #sub.plot(anc_mass, anc_phi, lw=2, ls='--', c='k') 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    #sub.legend(loc='upper right')
    subhalo_str = ''.join([
            'orphan_central_subhalo.scatter', 
            str(round(scatter, 2)), 
            '.source-', source
            ])

    fig.savefig(
            ''.join([
                'figure/'
                'qaplot_', subhalo_str, '.png'
                ]), 
            bbox_inches='tight'
            )
    #plt.show()
    plt.close()


if __name__=='__main__': 
    for scat in [ 0.0, 0.2 ]: 
        subh = CentralSubhalos() 
        subh.build_catalogs(snapshots = range(1, 21), scatter=scat, source='li-march')
        del subh
    #test_orphan_subhalo_smf(scatter=0.0, source='li-march')
