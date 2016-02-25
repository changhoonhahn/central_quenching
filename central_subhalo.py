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

    def read(self, snapshot, scatter=0.0, source='li-drory-march' ): 
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

            if i_snap == np.max(self.snapshots):     # ancestor index of subhalo
                ancestor = np.repeat(-999, len(central_indices))
                grp.create_dataset('ancestor'+str(np.max(self.snapshots)), data=ancestor)
            else: 
                ii_snap = np.max(self.snapshots) 
                ancestor = wetzel_util.utility_catalog.indices_tree(
                        sub, i_snap, ii_snap, central_indices)
                grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
                #for ii_snap in range(i_snap+1, np.max(self.snapshots)+1): 
                #    print i_snap, ii_snap
                #    ancestor = wetzel_util.utility_catalog.indices_tree(
                #            sub, i_snap, ii_snap, central_indices)
                #    grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
            
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
                ii_snap = np.max(self.snapshots)
                ancestor = wetzel_util.utility_catalog.indices_tree(
                        sub, i_snap, ii_snap, range(len(mhalo)))
                grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
                #for ii_snap in range(i_snap+1, np.max(self.snapshots)+1): 
                #    ancestor = wetzel_util.utility_catalog.indices_tree(
                #            sub, i_snap, ii_snap, range(len(mhalo)))
                #    grp.create_dataset('ancestor'+str(ii_snap), data=ancestor)
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

def test_descendant_subhalo_smf(nsnap_ancestor = 20, scatter = 0.0, source='li-drory-march'): 
    '''
    SMF for central subhalos that have ancestors
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # Central subhalo at nsnap_ancestor
    subh = CentralSubhalos()
    subh.read(nsnap_ancestor, scatter=scatter, source=source)
    ancestor_index = subh.index
    ancestor_mass = subh.mass

    nomass = np.where(subh.mass == 0.0)
    nomass_anc_index = subh.index[nomass]
    nomass_anc_mass = subh.mass[nomass]

    massive = np.where(subh.mass > 0.0)
    massive_anc_index = subh.index[massive]
    massive_anc_mass = subh.mass[massive]

    for i_snap in [1, 5, 10, 15]:

        if i_snap == 1: 
            massive_label = 'Ancestor M* > 0'
            nomass_label = 'Ancesotr M* = 0'
            label = 'Total'
        else: 
            massive_label = None 
            nomass_label = None
            label = None

        smf = SMF()

        # total central subhalo SMF
        subh = CentralSubhalos()
        subh.read(i_snap, scatter=scatter, source=source)
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        sub.plot(subh_mass, subh_phi, lw=4, c=pretty_colors[i_snap], alpha=0.5, label=label) 
        # SMF of central subhalos who's ancestors are massive centrals at snapshot 20
        if i_snap == 1: 
            has_ancestor, has_descendant = intersection_index(subh.ancestor20, massive_anc_index)
        else: 
            has_ancestor, has_d = intersection_index(subh.ancestor20, massive_anc_index)
        print 'Snapshot = ', i_snap, len(has_ancestor)
        mass, phi = smf.smf(subh.mass[has_ancestor])
        sub.plot(mass, phi, lw=2, ls='--', c=pretty_colors[i_snap], label=massive_label) 
        
        # SMF of central subhalos who's ancestors are non massive centrals at snapshot 20
        if i_snap == 1: 
            has_ancestor, has_descendant = intersection_index(subh.ancestor20, nomass_anc_index)
        else: 
            has_ancestor, has_d = intersection_index(subh.ancestor20, nomass_anc_index)
        mass, phi = smf.smf(subh.mass[has_ancestor])
        sub.plot(mass, phi, lw=2, ls='-.', c=pretty_colors[i_snap], label=nomass_label) 

    anc_mass, anc_phi = smf.smf(ancestor_mass)
    sub.plot(anc_mass, anc_phi, lw=4, c='gray') 
    print 'Total number of ancestors ', len(ancestor_mass)
    #print 'With descendants at snapshot 1', len(ancestor_mass[has_descendant])
    #anc_massive = np.where(ancestor_mass[has_descendant] > 0.)
    #print 'With mass greater than 0', len(ancestor_mass[has_descendant[anc_massive]])
    #anc_mass, anc_phi = smf.smf(ancestor_mass[has_descendant])
    #sub.plot(anc_mass, anc_phi, ls='--', lw=4, c='k') 
    #anc_mass, anc_phi = smf.smf(ancestor_mass[has_descendant[anc_massive]])
    #sub.plot(anc_mass, anc_phi, ls='--', lw=2, c='green') 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    sub.legend(loc='upper right')
    fig.savefig(
            ''.join([
                'figure/'
                'qaplot_',
                'descendant_central_subhalo.scatter', 
                str(round(scatter, 2)), 
                '.source-', source,
                '.png'
                ]), 
            bbox_inches='tight'
            )
    #plt.show()
    plt.close()

def test_nomass_descendant_subhalo_smf(nsnap_ancestor = 20, scatter = 0.0, source='li-drory-march'): 
    '''
    SMF for central subhalos that have ancestors with no stellar mass at nsnap=20
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # Central subhalo at nsnap_ancestor
    subh = CentralSubhalos()
    subh.read(nsnap_ancestor, scatter=scatter, source=source)
    nomass = np.where(subh.mass == 0.)
    ancestor_index = subh.index[nomass]
    ancestor_mass = subh.mass[nomass]
    print 'Total number of 0 mass ancestors ', len(ancestor_mass)

    for i_snap in [1, 5, 10, 15]:
        smf = SMF()
        subh = CentralSubhalos()
        subh.read(i_snap, scatter=scatter, source=source)
    
        # total central subhalo SMF
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        sub.plot(subh_mass, subh_phi, lw=4, c=pretty_colors[i_snap], alpha=0.5) 

        # SMF of central subhalos with ancestors
        # has_ancestor = np.where(subh.ancestor20 > 0)[0]
        # SMF of central subhalos who's ancestors are centrals at snapshot 20
        if i_snap == 1: 
            has_ancestor, has_descendant = intersection_index(subh.ancestor20,ancestor_index)
        else: 
            has_ancestor, has_d = intersection_index(subh.ancestor20,ancestor_index)
        #print subh.ancestor20[has_ancestor]
        #print ancestor_index[has_descendant]
        #massive = np.where(subh.mass[has_ancestor] > 0.0)
        #print 'Snapshot = ', i_snap, len(has_ancestor), len(subh.mass[has_ancestor[massive]])

        nonorph_mass, nonorph_phi = smf.smf(subh.mass[has_ancestor])
        sub.plot(nonorph_mass, nonorph_phi, lw=2, ls='-.', c=pretty_colors[i_snap]) 
        
        mass, phi = smf.smf(subh.mass)
        sub.plot(mass, phi, lw=2, c=pretty_colors[i_snap]) 

    #anc_mass, anc_phi = smf.smf(ancestor_mass)
    #sub.plot(anc_mass, anc_phi, lw=4, c='gray') 
    print 'With descendants at snapshot 1', len(ancestor_mass[has_descendant])
    anc_massive = np.where(ancestor_mass[has_descendant] > 0.)
    print 'With mass greater than 0', len(ancestor_mass[has_descendant[anc_massive]])
    #anc_mass, anc_phi = smf.smf(ancestor_mass[has_descendant])
    #sub.plot(anc_mass, anc_phi, ls='--', lw=4, c='k') 
    #anc_mass, anc_phi = smf.smf(ancestor_mass[has_descendant[anc_massive]])
    #sub.plot(anc_mass, anc_phi, ls='--', lw=2, c='green') 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    #sub.legend(loc='upper right')
    fig.savefig(
            ''.join([
                'figure/'
                'qaplot_',
                'nomass_descendant_central_subhalo.scatter', 
                str(round(scatter, 2)), 
                '.source-', source,
                '.png'
                ]), 
            bbox_inches='tight'
            )
    #plt.show()
    plt.close()

def test_orphan_subhalo_smf(nsnap_ancestor = 20, scatter = 0.0, source='li-drory-march'): 
    '''
    SMF for orphan central subhalos that do not have ancestors at nsnap_ancestor snapshot
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # Central subhalo at nsnap_ancestor
    subh = CentralSubhalos()
    subh.read(nsnap_ancestor, scatter=scatter, source=source)
    ancestor_index = subh.index
    ancestor_mass = subh.mass

    for i_snap in [1, 5, 10, 15]:
        smf = SMF()
        subh = CentralSubhalos()
        subh.read(i_snap, scatter=scatter, source=source)
    
        # total central subhalo SMF
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        sub.plot(subh_mass, subh_phi, lw=4, c=pretty_colors[i_snap], alpha=0.5) 

        # SMF of central subhalos who's ancestors are centrals at snapshot 20
        if i_snap == 1: 
            label = 'Orphans'
        else: 
            label = None
        orphan = np.invert(np.in1d(subh.ancestor20, ancestor_index))
        orph_mass, orph_phi = smf.smf(subh.mass[orphan])
        sub.plot(orph_mass, orph_phi, lw=2, ls='--', c=pretty_colors[i_snap], label=label) 

    anc_mass, anc_phi = smf.smf(ancestor_mass)
    sub.plot(anc_mass, anc_phi,  lw=4, c='gray', label='All Centrals') 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    sub.legend(loc='upper right')
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
    plt.close()

def test_ancestors_without_descendant(nsnap_ancestor = 20, scatter = 0.0, source='li-drory-march'): 
    '''
    What happens to ancestors that do not have descendants at nsnap = 1


    Notes
    -----
    * More than 50% of the subhalos at nsnap=20 do not have descendants at nsnap = 1. 
        What happens to them? 
    * Do they become satellites? --> Some become satellites, others 
        with no mass subhalos don't stay in the catalog at all
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # Central subhalo at nsnap_ancestor (ancesotr)
    anc = CentralSubhalos()
    anc.read(nsnap_ancestor, scatter=scatter, source=source)
    
    # Centrals subhalo at nsnap = 1 (descendant)
    des = CentralSubhalos()
    des.read(1, scatter=scatter, source=source)
    
    no_descendants = np.invert(np.in1d(anc.index, des.ancestor20))
    massive_nodescendants = np.where(anc.mass[no_descendants] > 0.0)
    print 'N_SH no descendants ', len(anc.mass[no_descendants])
    print 'N_SH total ', len(anc.mass)
    print np.float(len(anc.mass[no_descendants]))/np.float(len(anc.mass))

    no_des_index = anc.index[no_descendants][massive_nodescendants][:5]
    print anc.mass[no_descendants][massive_nodescendants][:5]
    for isnap in range(1, nsnap_ancestor)[::-1]: 
        i_des = CentralSubhalos()
        i_des.read(isnap, scatter=scatter, source=source)
    
        #print isnap, np.in1d(no_des_index, i_des.ancestor20)

        #if not np.in1d(no_des_index, i_des.ancestor20)[0]: 
        des_sh = Subhalos()
        des_sh.read(isnap, scatter=scatter, source=source)
        #if np.sum(np.in1d(des_sh.ancestor20, no_des_index)) != len(no_des_index): 
        #    raise ValueError
        print des_sh.ilk[np.in1d(des_sh.ancestor20, no_des_index)][np.argsort(des_sh.ancestor20[np.in1d(des_sh.ancestor20, no_des_index)])]
        print 'M* ', des_sh.mass[np.in1d(des_sh.ancestor20, no_des_index)][np.argsort(des_sh.ancestor20[np.in1d(des_sh.ancestor20, no_des_index)])]
        print 'M_halo ', getattr(des_sh, 'halo.m')[np.in1d(des_sh.ancestor20, no_des_index)][np.argsort(des_sh.ancestor20[np.in1d(des_sh.ancestor20, no_des_index)])]
        #print des_sh.mass[np.in1d(des_sh.index, no_des_index)]
        #np.in1d(i_des.ancestor20, no_des_index)

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    sub.legend(loc='upper right')
    #subhalo_str = ''.join([
    #        'orphan_central_subhalo.scatter', 
    #        str(round(scatter, 2)), 
    #        '.source-', source
    #        ])

    #fig.savefig(
    #        ''.join([
    #            'figure/'
    #            'qaplot_', subhalo_str, '.png'
    #            ]), 
    #        bbox_inches='tight'
    #        )
    #plt.show()
    #plt.close()

if __name__=='__main__': 
    #test_subhalos_smf(scatter=0.0, source='li-march', type='all')
    #test_subhalos_smf(scatter=0.0, source='li-march', type='central')
    #for scat in [ 0.0, 0.2 ]: 
    #    subh = Subhalos() 
    #    subh.build_catalogs(snapshots = range(1, 21), scatter=scat, source='li-march')
    #    del subh
    #test_orphan_subhalo_smf(scatter=0.0, source='li-march')
    test_descendant_subhalo_smf(scatter=0.0, source='li-march')
    #test_nomass_descendant_subhalo_smf(scatter=0.0, source='li-march')
    #test_ancestors_without_descendant(scatter=0.0, source='li-march')

