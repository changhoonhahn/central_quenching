'''


Code for analyzing the satellit quenching timescales


'''
import h5py
import pickle 
import numpy as np
import os.path
from scipy.interpolate import interp1d

# --- Local --- 
import sham_hack as sham
import util.util as Util
from util.util import code_dir
import sfr_evol 
from gal_prop import Fq
from gal_prop import Ssfr

from abcee import PlotABC
from abcee import ReadABCrun
from inherit import ABCInherit
from observations import GroupCat

#from treepm import subhalo_io 
import subhalo_io_hack as subhalo_io
from utilities import utility as wetzel_util
from plotting.plots import PlotFq
from plotting.plots import PlotSSFR

import matplotlib.pyplot as plt 


class SatelliteSubhalos(object): 
    def __init__(self, **kwargs): 
        ''' Class that describes the Satellite Subhalo Catalogs generated from 
        Andrew Wetzel's TreePM merger tree code. This class describes
        all subhalos (both centrals and satellites) 

        '''
        self.kwargs = kwargs

    def File(self, **spec_kwargs): 
        ''' File name of output file given snapshot 
        '''
        # write to hdf5 file 
        subsham_file = ''.join([ 
            code_dir(), 'dat/wetzel_tree/', 
            'subhalo_sham', 
            '.satellite',
            '.snapshot', 
            str(spec_kwargs['snapshot']), 
            self._file_spec(**spec_kwargs), 
            '.hdf5']) 
        return subsham_file 
    
    def _file_spec(self, **spec_kwargs): 
        ''' file specifier string that describes the key choices of parameters in the file
        '''
        spec_str = ''.join([
            '.scatter', str(spec_kwargs['scatter']),
            '.', spec_kwargs['source']])
        return spec_str 

    def build_catalogs(self, scatter = 0.0, source='li-drory-march'): 
        ''' Build central subhalo snapshot catalogs using TreePM code. The kwargs specify the 
        properties of the central subhalos
        '''
        self.source = source
        self.scatter = scatter
    
        # read in TreePM Subhalo snapshots from z~0.0502 to z~1.0833
        # note: 250 here is the simulation box length (Mpc/h comoving)
        sub = subhalo_io.Treepm.read(
                'subhalo', 
                250, 
                zis = range(1,33)
                ) 

        # assign M* to subhalos using SHAM
        sham.assign(
                sub, 
                'm.star', 
                scat = self.scatter, 
                dis_mf = 0.0, 
                source = self.source,
                zis = range(1,33) 
                ) 
        # satellite indices 
        sat_indices = wetzel_util.utility_catalog.indices_ilk(sub[1], ilk = 'sat') 
        print 'Number of Satellite Subhalos in Snapshot 1 =', len(sat_indices)
        print 'out of ', len(sub[1]['halo.m'])
        # subhalo properties
        mhalo     = (sub[1]['halo.m'])[sat_indices]     # M_halo
        mhalo_max = (sub[1]['m.max'])[sat_indices]     # M_halo
        mstar     = (sub[1]['m.star'])[sat_indices]     # M* of subhalo
        pos       = (sub[1]['pos'])[sat_indices]        # position of subhalo
        ilk       = (sub[1]['ilk'])[sat_indices]        # classification flag in case we want to consider centrals and virtual centrals separately 
        zi_infall  = (sub[1]['inf.first.zi'])[sat_indices] # infall snapshot 

        spec_kwargs = {
                'scatter': self.scatter, 
                'snapshot': 1, 
                'source': self.source
                }
        subsham_sat_file = self.File(**spec_kwargs)    # output name

        f = h5py.File(subsham_sat_file, 'w') 
        grp = f.create_group('cenque_data') 
        
        grp.attrs['snapshot'] = 1 
        grp.create_dataset('index', data=sat_indices)
        grp.create_dataset('pos', data=pos) 
        grp.create_dataset('ilk', data=ilk)
        grp.create_dataset('mass', data=mstar)
        grp.create_dataset('halo.m', data=mhalo)
        grp.create_dataset('halo.m.max', data=mhalo_max)
        grp.create_dataset('inf.first.zi', data=zi_infall)

        t_infall = np.repeat(-999., len(zi_infall))
        z_infall = np.repeat(-999., len(zi_infall))
        mass_infall = np.repeat(-999., len(zi_infall))

        for i_snap in range(2, 33): 
            # ancestor index of subhalo for ancestor in snapshot nsnap_ancestor
            ancestor = wetzel_util.utility_catalog.indices_tree(
                    sub, 1, i_snap, sat_indices)
            grp.create_dataset('ancestor'+str(i_snap), data=ancestor)
            
            # ancestor SHAM stellar mass and halo masses  
            has_ancestor, has_descendant = Util.intersection_index(
                    ancestor, np.arange(len(sub[i_snap]['m.star'])))
            ancestor_mass = np.repeat(-999., len(sat_indices))
            ancestor_mass[has_ancestor] = sub[i_snap]['m.star'][has_descendant]
            grp.create_dataset('ancestor'+str(i_snap)+'_m.star', data=ancestor_mass)
            
            # keep track of 
            infalls_now = np.where(zi_infall == i_snap)[0]
            t_infall[infalls_now] = Util.get_t_nsnap(i_snap)
            z_infall[infalls_now] = Util.get_z_nsnap(i_snap) 

            inf_ancestor, inf_descendant = Util.intersection_index(
                    ancestor[infalls_now], np.arange(len(sub[i_snap]['m.star'])))
            mass_infall[infalls_now[inf_ancestor]] = sub[i_snap]['m.star'][inf_descendant]
            
        grp.create_dataset('inf.first.t', data=t_infall)
        grp.create_dataset('inf.first.z', data=z_infall)
        grp.create_dataset('inf.first.mass', data=mass_infall)
        f.close() 
        return None
    
    def Read(self, scatter=0.0, source='li-drory-march'): 
        ''' Read in the hdf5 file that contains the Central Subhalo data and import 
        it appropriately to the class object structure
        '''
        spec_kwargs = {
                'scatter': scatter, 
                'snapshot': 1, 
                'source' : source
                }
        snapshot_file = self.File(**spec_kwargs)

        self.file_name = snapshot_file
        if not os.path.isfile(snapshot_file): 
            print snapshot_file, ' does NOT exist'
            print 'Now building'
            self.build_catalogs(
                    scatter=scatter, 
                    source=source
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


class SGPop(object): 
    def __init__(self, **kwargs): 
        ''' Class object that descirbes the satellite galaxy population generated 
        from the SatelliteGalaxy class, which is a hacked together wrapper for 
        Andrew Wetzel's TreePM Subhalo Catalog. 
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
        self.first_infall_nsnap = None     # first infall snapshot 
        self.first_infall_z = None     # first infall snapshot 
        self.first_infall_t = None     # first infall snapshot 

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

        self.zsnap = None           # z_snapshot
        self.t_cosmic = None    # t_cosmic for snapshot
        self.t_step = None      # t_cosmic step 

        self.data_columns = None    # default has no data

        self.ssfr_dist = None

    def ImportSubhalo(self, subhalo_prop=None): 
        ''' Import  Satellite Subhalo Snapshots with SHAM M* through the CentralSubhalos class 
        object.
        
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
        
        satsub = SatelliteSubhalos()
        satsub.Read(
                scatter=self.subhalo_prop['scatter'],   # scatter
                source=self.subhalo_prop['source'],      # SMF source
                )
        print satsub.file_name
    
        # pass SatalliteSubhalos class attributes to 
        self.data_columns = self.get_datacol() 
        for col in satsub.data_columns: 
            if col == 'halo.m.max': 
                newcol = 'halo_mass'
            elif col == 'index': 
                newcol = 'snap_index'
            elif col == 'halo.m': 
                continue 
            elif col == 'inf.first.zi': 
                newcol = 'first_infall_nsnap'
            elif col == 'inf.first.t': 
                newcol = 'first_infall_t'
            elif col == 'inf.first.z': 
                newcol = 'first_infall_z'
            elif col == 'inf.first.mass': 
                newcol = 'first_infall_mass'
            else: 
                newcol = col
            setattr(self, newcol, getattr(satsub, col))
            if newcol not in self.data_columns: 
                self.data_columns.append(newcol)

        # Meta data 
        self.metadata = ['nsnap', 'zsnap', 't_cosmic', 't_step', 'subhalo_prop']
        # get snapshot redshift/cosmic time data using Andrew's table
        n_snaps, z_snap, t_snap, t_wid = np.loadtxt(
                Util.snapshottable(), 
                unpack=True, 
                usecols=[0, 2, 3, 4]
                )
        self.nsnap = 1 
        self.zsnap = z_snap[(n_snaps.tolist()).index(self.nsnap)]        # redshift of snapshot
        self.t_cosmic = t_snap[(n_snaps.tolist()).index(self.nsnap)] # t_cosmic of snapshot 
        self.t_step = t_wid[(n_snaps.tolist()).index(self.nsnap)]    # Gyrs until next snapshot
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
                subhalo_prop=self.subhalo_prop
                )
        cgpop_file = ''.join([Util.code_dir(), 
            'dat/galpop/', 
            'sgpop'
            '.snapshot', str(self.nsnap), 
            file_spec_str, 
            '.hdf5'
            ]) 

        return cgpop_file

    def _file_spec(self, subhalo_prop=None): 
        ''' file specifier depending on whether subhlao properties or sfr properties 
        have been specified 
        '''
        if subhalo_prop is None: 
            raise ValueError
    
        subhalo_str = ''.join([
            '.subhalo_', 
            'scatter', str(subhalo_prop['scatter']), '_', 
            'ancestor32_',
            subhalo_prop['source']
            ])
        return subhalo_str

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
   
    def _clean_initialize(self): 
        ''' Convenience method for sf_inherit.py where all the attributes are initialized
        '''
        if self.snap_index is None: 
            raise ValueError("nothing to clean here") 
    
        n_gal = len(self.snap_index)
        self.sfr      = np.repeat(-999., n_gal)
        self.ssfr     = np.repeat(-999., n_gal)
        self.min_ssfr = np.repeat(-999., n_gal)
        self.sfr_class = np.chararray(n_gal, itemsize=16)
        self.sfr_class[:] = ''

        return None
    

def AssignCenSFR(tf, abcrun=None, prior_name=None):
    ''' Given SGPop object, assign SFRs at infalling time based on the 
    median ABC Run of the Central Galaxy evolution 
    '''
    if abcrun is None: 
        raise ValueError
    if prior_name is None: 
        raise ValueError
    abcinh = ABCInherit(tf, abcrun=abcrun, prior_name=prior_name) 
    
    inh = abcinh.PickleLoad('all') 
    
    sg_obj = SGPop()
    sg_obj.ImportSubhalo(subhalo_prop=abcinh.sim_kwargs['subhalo_prop']) 
    sg_obj.sfms_prop = abcinh.sim_kwargs['sfr_prop']['sfms']
    sg_obj.abcrun = abcrun
    sg_obj.prior_name = prior_name 
    
    sg_obj.first_infall_sfr = np.repeat(-999., len(sg_obj.first_infall_nsnap))
    for i_snap in range(2, abcinh.sim_kwargs['nsnap_ancestor']): 
        infall_now = np.where(
                (sg_obj.first_infall_nsnap == i_snap) & 
                (sg_obj.first_infall_mass > 0.))[0]

        # satellite galaxies that end up massive but have 0. infalling mass 
        #massive_infall_zero =np.where(
        #        (sg_obj.first_infall_mass[infall_now] == 0.) & 
        #        (sg_obj.mass[infall_now] > 0.)
        #        ) 
        
        des_snap = inh[str(i_snap)]
        m_bins = np.arange(des_snap.mass.min()-0.1, des_snap.mass.max()+0.2, 0.2) 

        for i_m in range(len(m_bins)-1): 
            in_massbin = np.where(
                    (des_snap.mass >= m_bins[i_m]) & 
                    (des_snap.mass < m_bins[i_m+1]))

            infall_massbin = np.where(
                    (sg_obj.first_infall_nsnap == i_snap) & 
                    (sg_obj.first_infall_mass >= m_bins[i_m]) & 
                    (sg_obj.first_infall_mass < m_bins[i_m+1]))

            if len(infall_massbin[0]) == 0: 
                continue 

            des_ssfrs = des_snap.sfr[in_massbin] - des_snap.mass[in_massbin] 
            pdf_ssfr_i, ssfr_binedges = np.histogram(des_ssfrs, 
                    bins=np.arange(-13.0, -8.0, 0.2), normed=True) 

            cdf_ssfr_i = np.zeros(len(ssfr_binedges))
            cdf_ssfr_i[1:] = np.cumsum(pdf_ssfr_i * np.diff(ssfr_binedges))
            cdf_interp = interp1d(cdf_ssfr_i, ssfr_binedges, kind='linear') 

            rand_i = np.random.uniform(size=len(infall_massbin[0]))
            sg_obj.first_infall_sfr[infall_massbin] = cdf_interp(rand_i) + \
                    sg_obj.first_infall_mass[infall_massbin]
    return sg_obj


def EvolveSatSFR(sg_obj): 
    ''' Evolve the infalling SFR assigned based on central galaxy SFRs
    '''
    # initial SFR, SSFR
    sg_obj.sfr = np.repeat(-999., len(sg_obj.mass))
    sg_obj.ssfr = np.repeat(-999., len(sg_obj.mass))

    # First classify the galaxies into SF, Q, and Qing
    # this is so that t_Q,start = t_inf for Qing galaxies, 
    # Q galaxies are left alone
    # SF galaxies are affected by the delay time 
    infall = np.where(
            (sg_obj.first_infall_mass > 0.) &
            (sg_obj.first_infall_sfr > -999.) & 
            (sg_obj.first_infall_t >  5.7))
    print len(infall[0]), ' infall galaxies'

    fq_obj = Fq()  
    sfq_class = fq_obj.Classify(
            sg_obj.first_infall_mass[infall], 
            sg_obj.first_infall_sfr[infall], 
            sg_obj.first_infall_z[infall], 
            sg_obj.sfms_prop)
    
    sf_infall = (infall[0])[np.where(sfq_class == 'star-forming')] 
    q_infall = (infall[0])[np.where(sfq_class == 'quiescent')]

    q_mass = sg_obj.first_infall_mass[q_infall]
    q_sfr = sg_obj.first_infall_sfr[q_infall]
    q_ssfr = q_sfr - q_mass 

    q_cut_SSFR = sfr_evol.AverageLogSSFR_q_peak(q_mass) + \
            1.5 * sfr_evol.ScatterLogSSFR_q_peak(q_mass)

    qing_infall = q_infall[np.where(q_ssfr > q_cut_SSFR)]   # quenching 
    qq_infall = q_infall[np.where(q_ssfr <= q_cut_SSFR)]    # quenched

    # Quiescent @ infall 
    ''' The SSFR is preserved from infall, so effectively their 
    SFR decreases
    '''
    print len(qq_infall), ' quiescent galaxies' 
    sg_obj.ssfr[qq_infall] = sg_obj.first_infall_sfr[qq_infall] - \
            sg_obj.first_infall_mass[qq_infall]
    sg_obj.sfr[qq_infall] = sg_obj.mass[qq_infall] + sg_obj.ssfr[qq_infall]
    
    #print (sg_obj.mass[qq_infall])[np.where(sg_obj.ssfr[qq_infall] > -11.)]
    #plt.hist(sg_obj.ssfr[qq_infall], bins=20)
    #plt.show() 

    sg_obj.t_qstart = np.repeat(-999., len(sg_obj.mass))
    sg_obj.z_qstart = np.repeat(-999., len(sg_obj.mass))
    sg_obj.q_ssfr = np.repeat(-999., len(sg_obj.mass))

    # Quenching @ infall  
    ''' Quenching satellite galaxies skip the delay phase and immediately 
    start quenching when the simulation begins. 
    '''
    sg_obj.t_qstart[qing_infall] = sg_obj.first_infall_t[qing_infall]
    sg_obj.z_qstart[qing_infall] = sg_obj.first_infall_z[qing_infall]
    sg_obj.q_ssfr[qing_infall] = sfr_evol.AverageLogSSFR_q_peak(sg_obj.mass[qing_infall]) + \
            sfr_evol.ScatterLogSSFR_q_peak(sg_obj.mass[qing_infall]) * np.random.randn(len(qing_infall))
    
    dlogSFR_MS_M = 0.53 * (sg_obj.mass[qing_infall] - sg_obj.first_infall_mass[qing_infall]) 
    
    dlogSFR_MS_z = sg_obj.sfms_prop['zslope'] * (sg_obj.zsnap - sg_obj.first_infall_z[qing_infall])

    dlogSFR_Q = sfr_evol.DeltaLogSFR_quenching(
            sg_obj.t_qstart[qing_infall], 
            sg_obj.t_cosmic, 
            M_q=sg_obj.mass[qing_infall],
            tau_prop={'name': 'satellite'})

    sg_obj.sfr[qing_infall] = sg_obj.first_infall_sfr[qing_infall] + \
            dlogSFR_MS_M + dlogSFR_MS_z + dlogSFR_Q
    sg_obj.ssfr[qing_infall] = sg_obj.sfr[qing_infall] - sg_obj.mass[qing_infall]
    
    # deal with over quenching 
    overquenched = np.where(sg_obj.ssfr[qing_infall] < sg_obj.q_ssfr[qing_infall])
    sg_obj.ssfr[qing_infall[overquenched]] = sg_obj.q_ssfr[qing_infall[overquenched]]
    
    # Star-forming @ infall 
    sg_obj.t_qstart[sf_infall] = sg_obj.first_infall_t[sf_infall] + tDelay(sg_obj.mass[sf_infall]) 
    # if t_qstart > t_cosmic, then it does not quenching during the simualtion 
    sf_q = np.where(sg_obj.t_qstart[sf_infall] < sg_obj.t_cosmic) 
    sf_noq = np.where(sg_obj.t_qstart[sf_infall] >= sg_obj.t_cosmic) 
    sf_infall_q = sf_infall[sf_q]
    sf_infall_noq = sf_infall[sf_noq]

    sg_obj.z_qstart[sf_infall_q] = Util.get_zsnap(sg_obj.t_qstart[sf_infall_q]) 
    sg_obj.z_qstart[sf_infall_noq] = 999.
    
    # SF galaxies that quench
    sg_obj.q_ssfr[sf_infall_q] = \
            sfr_evol.AverageLogSSFR_q_peak(sg_obj.mass[sf_infall_q]) + \
            sfr_evol.ScatterLogSSFR_q_peak(sg_obj.mass[sf_infall_q]) * \
            np.random.randn(len(sf_infall_q))

    dlogSFR_MS_M = 0.53 * (sg_obj.mass[sf_infall_q] - sg_obj.first_infall_mass[sf_infall_q]) 
    
    dlogSFR_MS_z = sg_obj.sfms_prop['zslope'] * (sg_obj.zsnap - sg_obj.first_infall_z[sf_infall_q])
    dlogSFR_Q = sfr_evol.DeltaLogSFR_quenching(
            sg_obj.t_qstart[sf_infall_q], 
            sg_obj.t_cosmic, 
            M_q=sg_obj.mass[sf_infall_q],
            tau_prop={'name': 'satellite'})

    sg_obj.sfr[sf_infall_q] = sg_obj.first_infall_sfr[sf_infall_q] + \
            dlogSFR_MS_M + dlogSFR_MS_z + dlogSFR_Q
    sg_obj.ssfr[sf_infall_q] = sg_obj.sfr[sf_infall_q] - sg_obj.mass[sf_infall_q]
    # deal with over quenching 
    overquenched = np.where(sg_obj.ssfr[sf_infall_q] < sg_obj.q_ssfr[sf_infall_q])
    sg_obj.ssfr[sf_infall_q[overquenched]] = sg_obj.q_ssfr[sf_infall_q[overquenched]]


    # SF galaxies that do NOT quench
    dlogSFR_MS_M = 0.53 * (sg_obj.mass[sf_infall_noq] - sg_obj.first_infall_mass[sf_infall_noq]) 
    
    dlogSFR_MS_z = sg_obj.sfms_prop['zslope'] * (sg_obj.zsnap - sg_obj.first_infall_z[sf_infall_noq])
    
    sg_obj.sfr[sf_infall_noq] = sg_obj.first_infall_sfr[sf_infall_noq] + \
            dlogSFR_MS_M + dlogSFR_MS_z
    sg_obj.ssfr[sf_infall_noq] = sg_obj.sfr[sf_infall_noq] - sg_obj.mass[sf_infall_noq]
    return sg_obj 


def PlotEvolvedSat(sg_obj): 
    ''' Plot stuff for the evolved satellite population 
    '''
    infall = np.where(
            (sg_obj.first_infall_mass > 0.) &
            (sg_obj.first_infall_sfr > -999.) & 
            (sg_obj.first_infall_t >  5.7))
    
    # SSFR plot
    ssfr_plot = PlotSSFR() 
    ssfr_plot.plot(mass=sg_obj.mass[infall], ssfr=sg_obj.ssfr[infall], line_color=3)

    ssfr_plot.GroupCat(position='satellite') 
    ssfr_plot.set_axes() 

    ssfr_fig_name = ''.join(['figure/', 
        'SSFR.Satellite',
        '.ABC_posterior',
        '.', sg_obj.abcrun, 
        '.', sg_obj.prior_name, '_prior', 
        '.png'])
    ssfr_plot.save_fig(ssfr_fig_name)
    plt.close()

    # Quiescent fraction plot  
    fq_plot = PlotFq()
    fq_plot.plot(mass=sg_obj.mass[infall], sfr=sg_obj.sfr[infall], z=sg_obj.zsnap, line_color='r',
            sfms_prop=sg_obj.sfms_prop, label='SHAM Sat. Simulation')
    grpcat = GroupCat(Mrcut=18, position='satellite')
    grpcat.Read()
    fq_plot.plot(mass=grpcat.mass, sfr=grpcat.sfr, z=np.median(grpcat.z), line_color='k', line_style='--',
            sfms_prop=sg_obj.sfms_prop, label='Satellite Group Catalog')
    fq_plot.set_axes()
    fq_fig_name = ''.join(['figure/', 
        'Fq.Satellite',
        '.ABC_posterior',
        '.', sg_obj.abcrun, 
        '.', sg_obj.prior_name, '_prior', 
        '.png'])
    fq_plot.save_fig(fq_fig_name)
    plt.close() 

    return None 


def tDelay(Mstars): 
    ''' Hacked parameterizatiion of t_delay in Wetzel plot
    '''
    m_arr = np.array([6.51782e9, 1.13089e10, 2.18647e10, 3.90733e10, 6.74621e10, 1.08724e11, 1.56475e11])
    t_del = np.array([3.30321e0 , 3.30102e0 , 3.19772e0 , 2.91165e0 , 2.14061e0 , 1.27830e0 , 2.15064e-1])

    func = interp1d(np.log10(m_arr), t_del, kind='linear')
    
    within = np.where(
            (Mstars >= np.log10(m_arr.min())) & 
            (Mstars <= np.log10(m_arr.max())))
    
    output = np.zeros(len(Mstars))
    output[within] = func(Mstars[within])
    
    greater = np.where(Mstars > np.log10(m_arr.max())) 
    output[greater] = (t_del[-1] - t_del[-2])/(m_arr[-1] - m_arr[-2]) * (Mstars[greater] - m_arr[-1]) + t_del[-1]

    less = np.where(Mstars < np.log10(m_arr.min())) 
    output[less] = (t_del[0] - t_del[1])/(m_arr[0] - m_arr[1]) * (Mstars[less] - m_arr[0]) + t_del[0]

    neg = np.where(output < 0.) 
    output[neg] = 0.

    return output




if __name__=='__main__': 
    #for scat in [0.0, 0.2]: 
    #    subh = SatelliteSubhalos() 
    #    subh.build_catalogs(scatter=scat, source='li-march')
    #    del subh

    sat_test_tmp = AssignCenSFR(7, abcrun='multifq_wideprior', prior_name='updated')
    #pickle.dump(sat_test_tmp, open('/data1/hahn/pmc_abc/pickle/sat_test_tmp.p', 'wb'))
    #sat_test_tmp = pickle.load(open('/data1/hahn/pmc_abc/pickle/sat_test_tmp.p', 'rb'))
    evol_sat = EvolveSatSFR(sat_test_tmp)
    PlotEvolvedSat(evol_sat)
