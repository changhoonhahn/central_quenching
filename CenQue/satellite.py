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

import abcpmc
from abcpmc import mpi_util

#from treepm import subhalo_io 
import subhalo_io_hack as subhalo_io
from utilities import utility as wetzel_util
from plotting.plots import PlotFq
from plotting.plots import PlotSSFR

import corner
import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors 


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

    # assign SFR to galaxies that fell in before nsnap = 15
    old_infall = np.where(
            (sg_obj.first_infall_nsnap >= abcinh.sim_kwargs['nsnap_ancestor']) &
            (sg_obj.first_infall_mass > 0.))[0]
    
    Pq = np.random.uniform(size=len(old_infall))
    qfrac = Fq() 
    fq_oldinfall = qfrac.model(sg_obj.first_infall_mass[old_infall], sg_obj.first_infall_z[old_infall], lit='wetzel')

    # star-forming
    sfing = np.where(Pq > fq_oldinfall) 
    sig_sfms = sfr_evol.ScatterLogSFR_sfms(
            sg_obj.first_infall_mass[old_infall[sfing]], 
            sg_obj.first_infall_z[old_infall[sfing]], 
            sfms_prop=sg_obj.sfms_prop)
    mu_sfms = sfr_evol.AverageLogSFR_sfms(
            sg_obj.first_infall_mass[old_infall[sfing]], 
            sg_obj.first_infall_z[old_infall[sfing]], 
            sfms_prop=sg_obj.sfms_prop)
    sg_obj.first_infall_sfr[old_infall[sfing]] = mu_sfms + np.random.randn(len(sfing[0])) * sig_sfms

    # quiescent
    qed = np.where(Pq <= fq_oldinfall) 
    mu_ssfr_q = sfr_evol.AverageLogSSFR_q_peak(sg_obj.first_infall_mass[old_infall[qed]])
    sig_ssfr_q = sfr_evol.ScatterLogSSFR_q_peak(sg_obj.first_infall_mass[old_infall[qed]])

    ssfr_q = mu_ssfr_q + np.random.randn(len(qed[0])) * sig_ssfr_q 

    sg_obj.first_infall_sfr[old_infall[qed]] = ssfr_q + sg_obj.first_infall_mass[old_infall[qed]]

    return sg_obj


def EvolveSatSFR(sg_in, tqdelay_dict=None): 
    ''' Evolve the infalling SFR assigned based on central galaxy SFRs
    '''
    if tqdelay_dict == None: 
        raise ValueError
    
    # initalize
    sg_obj = SGPop()
    sg_obj.ilk = sg_in.ilk.copy()
    sg_obj.pos = sg_in.pos.copy()
    sg_obj.snap_index = sg_in.snap_index.copy()
    sg_obj.zsnap = sg_in.zsnap.copy() 
    sg_obj.t_cosmic = sg_in.t_cosmic.copy() 
    sg_obj.sfms_prop = sg_in.sfms_prop.copy()
    sg_obj.abcrun = sg_in.abcrun
    sg_obj.prior_name = sg_in.prior_name
    sg_obj.mass = sg_in.mass.copy()  
    sg_obj.sfr =  np.repeat(-999., len(sg_obj.mass))
    sg_obj.ssfr = np.repeat(-999., len(sg_obj.mass))
    sg_obj.halo_mass = sg_in.halo_mass.copy()
    sg_obj.first_infall_nsnap = sg_in.first_infall_nsnap.copy()
    sg_obj.first_infall_z = sg_in.first_infall_z.copy()  
    sg_obj.first_infall_t = sg_in.first_infall_t.copy()  
    sg_obj.first_infall_sfr = sg_in.first_infall_sfr.copy()  
    sg_obj.first_infall_mass = sg_in.first_infall_mass.copy()  
    sg_obj.t_qstart = np.repeat(-999., len(sg_obj.mass))
    sg_obj.z_qstart = np.repeat(-999., len(sg_obj.mass))
    sg_obj.q_ssfr = np.repeat(-999., len(sg_obj.mass))

    # First classify the galaxies into SF, Q, and Qing
    # this is so that t_Q,start = t_inf for Qing galaxies, 
    # Q galaxies are left alone
    # SF galaxies are affected by the delay time 
    infall = np.where(
            (sg_obj.first_infall_mass > 0.) &
            (sg_obj.first_infall_sfr > -999.))

    fq_obj = Fq()  
    sfq_class = fq_obj.Classify(
            sg_obj.first_infall_mass[infall], 
            sg_obj.first_infall_sfr[infall], 
            sg_obj.first_infall_z[infall], 
            sg_obj.sfms_prop)
     
    sf_infall = (infall[0])[np.where(sfq_class == 'star-forming')]  # starforming @ infall
    q_infall = (infall[0])[np.where(sfq_class == 'quiescent')]      # quiescent @ infall

    ssfr_final = sfr_evol.AverageLogSSFR_q_peak(sg_obj.mass[infall]) + \
            sfr_evol.ScatterLogSSFR_q_peak(sg_obj.mass[infall]) * np.random.randn(len(infall[0]))
    
    # sub-divide the quiescent @ infall population to those that are 
    # quenching @ infall + quiescent @ infall based on simple SSFR cut
    q_mass_infall = sg_obj.first_infall_mass[q_infall]
    q_sfr_infall = sg_obj.first_infall_sfr[q_infall]
    q_ssfr_infall = q_sfr_infall - q_mass_infall

    q_cut_SSFR = sfr_evol.AverageLogSSFR_q_peak(q_mass_infall) + \
            1.5 * sfr_evol.ScatterLogSSFR_q_peak(q_mass_infall)
    
    qing_infall = q_infall[np.where(q_ssfr_infall > q_cut_SSFR)]   # quenching 
    qq_infall = q_infall[np.where(q_ssfr_infall <= q_cut_SSFR)]    # quenched

    # Quiescent @ infall-----
    # The SSFR is preserved from infall, so effectively their SFR decreases
    sg_obj.ssfr[qq_infall] = sg_obj.first_infall_sfr[qq_infall] - \
            sg_obj.first_infall_mass[qq_infall]
    sg_obj.sfr[qq_infall] = sg_obj.mass[qq_infall] + sg_obj.ssfr[qq_infall]

    # Quenching @ infall-----
    # Quenching satellite galaxies skip the delay phase and immediately 
    # start quenching when the simulation begins. 
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
    #overquenched = np.where(sg_obj.ssfr[qing_infall] < sg_obj.q_ssfr[qing_infall])
    #sg_obj.ssfr[qing_infall[overquenched]] = sg_obj.q_ssfr[qing_infall[overquenched]]
    
    # Star-forming @ infall 
    if tqdelay_dict['name'] == 'hacked': 
        sg_obj.t_qstart[sf_infall] = sg_obj.first_infall_t[sf_infall] + \
                tDelay(sg_obj.mass[sf_infall], **tqdelay_dict) 
    else: 
        sg_obj.t_qstart[sf_infall] = sg_obj.first_infall_t[sf_infall] + \
                tDelay(sg_obj.first_infall_mass[sf_infall], **tqdelay_dict) 
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
    #overquenched = np.where(sg_obj.ssfr[sf_infall_q] < sg_obj.q_ssfr[sf_infall_q])
    #sg_obj.ssfr[sf_infall_q[overquenched]] = sg_obj.q_ssfr[sf_infall_q[overquenched]]

    # SF galaxies that do NOT quench
    dlogSFR_MS_M = 0.53 * (sg_obj.mass[sf_infall_noq] - sg_obj.first_infall_mass[sf_infall_noq]) 
    
    dlogSFR_MS_z = sg_obj.sfms_prop['zslope'] * (sg_obj.zsnap - sg_obj.first_infall_z[sf_infall_noq])
    
    sg_obj.sfr[sf_infall_noq] = sg_obj.first_infall_sfr[sf_infall_noq] + \
            dlogSFR_MS_M + dlogSFR_MS_z
    sg_obj.ssfr[sf_infall_noq] = sg_obj.sfr[sf_infall_noq] - sg_obj.mass[sf_infall_noq]


    # deal with overquenching all at once
    overquench = np.where(sg_obj.ssfr[infall] < ssfr_final)
    sg_obj.ssfr[infall[0][overquench]] = ssfr_final[overquench] 
    sg_obj.sfr[infall[0][overquench]] = ssfr_final[overquench] + sg_obj.mass[infall[0][overquench]]

    return sg_obj 


def PlotEvolvedSat(sg_obj): 
    ''' Plot stuff for the evolved satellite population 
    '''
    #infall = np.where(
    #        (sg_obj.first_infall_mass > 0.) &
    #        (sg_obj.first_infall_sfr > -999.) & 
    #        (sg_obj.first_infall_t >  5.7))
    infall = np.arange(len(sg_obj.first_infall_mass))
    
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


def tDelay(Mstars, **tqdelay_dict): 
    ''' Hacked parameterizatiion of t_delay in Wetzel plot
    '''
    if tqdelay_dict['name'] == 'hacked': 
        m_arr = np.array([6.59215e+9, 1.14512e+10, 2.21628e+10, 3.96133e+10, 6.82792e+10, 1.09802e+11, 1.57529e+11])
        t_del = np.array([3.48711e+0, 3.48546e+0, 3.40789e+0, 3.19310e+0, 2.61415e+0, 1.96669e+0, 1.16836e+0])

        func = interp1d(np.log10(m_arr), t_del, kind='linear')
        
        within = np.where(
                (Mstars >= np.log10(m_arr.min())) & 
                (Mstars <= np.log10(m_arr.max())))
        
        output = np.zeros(len(Mstars))
        output[within] = func(Mstars[within])

        m_arr = np.log10(m_arr)
        
        greater = np.where(Mstars > m_arr.max()) 
        output[greater] = (t_del[-1] - t_del[-2])/(m_arr[-1] - m_arr[-2]) * (Mstars[greater] - m_arr[-1]) + t_del[-1]

        less = np.where(Mstars < m_arr.min()) 
        output[less] = (t_del[0] - t_del[1])/(m_arr[0] - m_arr[1]) * (Mstars[less] - m_arr[0]) + t_del[0]

        neg = np.where(output < 0.) 
        output[neg] = 0.

    elif tqdelay_dict['name'] == 'explin': 
         
        output = -1.*10.**(tqdelay_dict['m'] + Mstars) + tqdelay_dict['b'] 
        
        neg = np.where(output < 0.) 
        output[neg] = 0.

    return output


def tRapid(Mstars): 
    '''
    '''
    m_arr = np.array([6.37119e+9, 1.05173e+10, 1.54219e+10, 2.55939e+10, 3.98025e+10, 6.06088e+10, 9.95227e+10, 1.57374e+11])
    t_rap = np.array([8.04511e-1, 7.06767e-1, 6.69173e-1, 5.48872e-1, 3.23308e-1, 2.85714e-1, 2.25564e-1, 1.80451e-1]) 
    func = interp1d(np.log10(m_arr), t_rap, kind='linear')
    within = np.where(
            (Mstars >= np.log10(m_arr.min())) & 
            (Mstars <= np.log10(m_arr.max())))
    output = np.zeros(len(Mstars))
    output[within] = func(Mstars[within])
        
    m_arr = np.log10(m_arr)
        
    greater = np.where(Mstars > m_arr.max()) 
    output[greater] = (t_rap[-1] - t_rap[-2])/(m_arr[-1] - m_arr[-2]) * (Mstars[greater] - m_arr[-1]) + t_rap[-1]

    less = np.where(Mstars < m_arr.min()) 
    output[less] = (t_rap[0] - t_rap[1])/(m_arr[0] - m_arr[1]) * (Mstars[less] - m_arr[0]) + t_rap[0]

    neg = np.where(output < 0.) 
    output[neg] = 0.

    return output


def Plot_tRapid(): 
    '''
    '''
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(1)
    sub = fig.add_subplot(111)

    m_val = np.array([6.37119e+9, 1.05173e+10, 1.54219e+10, 2.55939e+10, 3.98025e+10, 6.06088e+10, 9.95227e+10, 1.57374e+11])
    t_rap = np.array([8.04511e-1, 7.06767e-1, 6.69173e-1, 5.48872e-1, 3.23308e-1, 2.85714e-1, 2.25564e-1, 1.80451e-1]) 
    sub.scatter(m_val, t_rap, c='k', s=16)

    m_arr = np.arange(9.5, 12.0, 0.1) 
    sub.plot(10**m_arr, tRapid(m_arr), c=pretty_colors[1])
    sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop={'name':'satellite'}), c=pretty_colors[3])
    
    # x-axis
    sub.set_xscale('log') 
    sub.set_xlim([5*10**9, 2*10**11.]) 
    # y-axis
    sub.set_ylim([0.0, 4.0])
    sub.set_ylabel(r'$\mathtt{t_{Q, rapid}}$', fontsize=25) 

    plt.show()
    plt.close()
    return None



def ABC(T, eps_input, Npart=1000, cen_tf=None, cen_prior_name=None, cen_abcrun=None):
    ''' ABC-PMC implementation. 

    Parameters
    ----------
    T : (int) 
        Number of iterations

    eps_input : (float)
        Starting epsilon threshold value 

    N_part : (int)
        Number of particles

    prior_name : (string)
        String that specifies what prior to use.

    abcrun : (string)
        String that specifies abc run information 
    '''
    abcinh = ABCInherit(cen_tf, abcrun=cen_abcrun, prior_name=cen_prior_name) 

    # Data (Group Catalog Satellite fQ)
    grpcat = GroupCat(Mrcut=18, position='satellite')
    grpcat.Read()
    qfrac = Fq() 
    m_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
    M_mid = 0.5 * (m_bin[:-1] + m_bin[1:]) 

    sfq = qfrac.Classify(grpcat.mass, grpcat.sfr, np.median(grpcat.z), 
            sfms_prop=abcinh.sim_kwargs['sfr_prop']['sfms'])
    ngal, dum = np.histogram(grpcat.mass, bins=m_bin)
    ngal_q, dum = np.histogram(grpcat.mass[sfq == 'quiescent'], bins=m_bin)
    data_sum = [M_mid, ngal_q.astype('float')/ngal.astype('float')]
    
    # Simulator 
    cen_assigned_sat_file = ''.join(['/data1/hahn/pmc_abc/pickle/', 
        'satellite', '.cenassign', '.', cen_abcrun, '_ABC', '.', cen_prior_name, '_prior', '.p'])
    
    if not os.path.isfile(cen_assigned_sat_file): 
        sat_cen = AssignCenSFR(cen_tf, abcrun=cen_abcrun, prior_name=cen_prior_name)
        pickle.dump(sat_cen, open(cen_assigned_sat_file, 'wb'))
    else: 
        sat_cen = pickle.load(open(cen_assigned_sat_file, 'rb'))

    def Simz(tt):       # Simulator (forward model) 
        tqdel_dict = {'name': 'explin', 'm': tt[0] , 'b': tt[1]}

        sat_evol = EvolveSatSFR(sat_cen, tqdelay_dict=tqdel_dict)

        sfq_sim = qfrac.Classify(sat_evol.mass, sat_evol.sfr, sat_evol.zsnap, 
                sfms_prop=sat_evol.sfms_prop)
        ngal_sim, dum = np.histogram(sat_evol.mass, bins=m_bin)
        ngal_q_sim, dum = np.histogram(sat_evol.mass[sfq_sim == 'quiescent'], bins=m_bin)
        sim_sum = [M_mid, ngal_q_sim.astype('float')/ngal_sim.astype('float')]
        return sim_sum

    # Priors
    prior_min = [-11.75, 2.]
    prior_max = [-10.25, 4.]
    prior = abcpmc.TophatPrior(prior_min, prior_max)    # ABCPMC prior object

    def rho(simum, datum):  
        datum_dist = datum[1]
        simum_dist = simum[1]
        drho = np.sum((datum_dist - simum_dist)**2)
        return drho
    
    abcrun_flag = cen_abcrun + '_central'

    theta_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'Satellite.tQdelay.theta_t', str(pewl), '_', abcrun_flag, '.dat']) 
    w_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'Satellite.tQdelay.w_t', str(pewl), '_', abcrun_flag, '.dat']) 
    dist_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'Satellite.tQdelay.dist_t', str(pewl), '_', abcrun_flag, '.dat']) 
    eps_file = ''.join([code_dir(), 
        'dat/pmc_abc/Satellite.tQdelay.epsilon_', abcrun_flag, '.dat'])
   
    eps = abcpmc.ConstEps(T, eps_input)
    try:
        mpi_pool = mpi_util.MpiPool()
        abcpmc_sampler = abcpmc.Sampler(
                N=Npart,                # N_particles
                Y=data_sum,             # data
                postfn=Simz,            # simulator 
                dist=rho,           # distance function  
                pool=mpi_pool)  
    except AttributeError: 
        abcpmc_sampler = abcpmc.Sampler(
                N=Npart,                # N_particles
                Y=data_sum,             # data
                postfn=Simz,            # simulator 
                dist=rho)           # distance function  
    abcpmc_sampler.particle_proposal_cls = abcpmc.ParticleProposal

    pools = []
    f = open(eps_file, "w")
    f.close()
    eps_str = ''
    for pool in abcpmc_sampler.sample(prior, eps, pool=None):
        print '----------------------------------------'
        print 'eps ', pool.eps
        new_eps_str = '\t'+str(pool.eps)+'\n'
        if eps_str != new_eps_str:  # if eps is different, open fiel and append 
            f = open(eps_file, "a")
            eps_str = new_eps_str
            f.write(eps_str)
            f.close()

        print("T:{0},ratio: {1:>.4f}".format(pool.t, pool.ratio))
        print eps(pool.t)

        # write theta, weights, and distances to file 
        np.savetxt(theta_file(pool.t), pool.thetas, 
            header='tQdelay_slope, tQdelay_offset')
        np.savetxt(w_file(pool.t), pool.ws)
        np.savetxt(dist_file(pool.t), pool.dists)
    
        # update epsilon based on median thresholding 
        eps.eps = np.median(pool.dists)
        pools.append(pool)

    return pools 


def PlotABC_Corner(tf, cen_abcrun=None):
    '''
    '''
    abcrun_flag = cen_abcrun + '_central'
    theta_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'Satellite.tQdelay.theta_t', str(pewl), '_', abcrun_flag, '.dat']) 
   
    theta = np.loadtxt(theta_file(tf)) 
    med_theta = [np.median(theta[:,i]) for i in range(len(theta[0]))]

    params = [ 'tqdelay_slope', 'tqdelay_offset']


    prior_min = [-11.75, 2.]
    prior_max = [-10.25, 4.]
    #prior_min = [-12.0, 3.]
    #prior_max = [-11.0, 5.]
    prior_range = [(prior_min[i], prior_max[i]) for i in range(len(prior_min))]

    fig = corner.corner(
            theta,
            truths=med_theta,
            truth_color='#ee6a50',
            labels=['$t_{Q,\;delay}$ Slope', '$t_{Q,\;delay}$ offset'],
            label_kwargs={'fontsize': 15},
            range=prior_range,
            quantiles=[0.16,0.5,0.84],
            show_titles=True,
            title_args={"fontsize": 12},
            plot_datapoints=True,
            fill_contours=True,
            levels=[0.68, 0.95],
            color='b',
            bins=20,
            smooth=1.0)
    
    fig_file = ''.join(['figure/',
        'Satellite',
        '.Corner',
        '.tQdelayABC',
        '.', cen_abcrun, '_central.png'])
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def PlotABC_tQdelay(tf, cen_abcrun=None): 
    '''
    '''
    abcrun_flag = cen_abcrun + '_central'
    theta_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'Satellite.tQdelay.theta_t', str(pewl), '_', abcrun_flag, '.dat']) 
   
    theta = np.loadtxt(theta_file(tf)) 

    m_arr = np.arange(7.5, 12.0, 0.1) 
    tdels = [] 
    for i in xrange(len(theta)): 
        i_tqdelay_dict = {'name': 'explin', 'm': theta[i][0], 'b': theta[i][1]}
        tdels.append(tDelay(m_arr, **i_tqdelay_dict))
    
    tdels = np.array(tdels)
    a, b, c, d, e = np.percentile(tdels, [2.5, 16, 50, 84, 97.5], axis=0)

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(7,7))
    sub = fig.add_subplot(111)

    mstar = np.array([1.56781e+11, 1.01089e+11, 6.27548e+10, 3.98107e+10, 2.49831e+10, 1.57633e+10, 9.89223e+9, 6.37831e+9])
    tqdel_high = np.array([1.36456e+0 , 2.53148e+0 , 3.06686e+0 , 3.52693e+0 , 3.62625e+0 , 3.99574e+0, 3.99915e+0 ,3.99915e+0])
    tqdel_low = np.array([8.15728e-1, 1.98257e+0, 2.29240e+0, 2.82772e+0, 2.99466e+0, 2.72548e+0, 2.99016e+0, 2.68333e+0])
    sub.fill_between(mstar, tqdel_low, tqdel_high, color=pretty_colors[0], edgecolor="none", label='Roughly Wetzel+(2013)')

    #sub.plot(10**m_arr, tDelay(m_arr, **abc_tqdelay_dict), c=pretty_colors[3], label='ABC Best-fit')
    
    sub.fill_between(10**m_arr, a, e, color=pretty_colors[3], alpha=0.5, edgecolor="none", 
            label='ABC Best-fit')
    sub.fill_between(10**m_arr, b, d, color=pretty_colors[3], alpha=1., edgecolor="none")
    sub.plot(10**m_arr, tDelay(m_arr, **{'name': 'hacked'}))

    sub.set_ylim([0.0, 4.0])
    sub.set_ylabel(r'$\mathtt{t_{Q, delay}}$', fontsize=25) 
    sub.set_xlim([5*10**9, 2*10**11.]) 
    sub.set_xscale("log") 
    sub.set_xlabel(r'$\mathtt{M}_*$', fontsize=25) 
    sub.legend(loc='lower left') 
    
    fig.savefig('figure/sallite_tqdelay_'+cen_abcrun+'.png', bbox_inches='tight') 
    return None


def PlotABC_EvolvedSat(tf, cen_tf=7, cen_abcrun=None, cen_prior_name=None): 
    ''' Plot stuff for the evolved satellite population 
    '''
    # Simulator 
    cen_assigned_sat_file = ''.join(['/data1/hahn/pmc_abc/pickle/', 
        'satellite', '.cenassign', 
        '.', cen_abcrun, '_ABC', 
        '.', cen_prior_name, '_prior', '.p'])
    if not os.path.isfile(cen_assigned_sat_file): 
        sat_cen = AssignCenSFR(cen_tf, abcrun=cen_abcrun, prior_name=cen_prior_name)
        pickle.dump(sat_cen, open(cen_assigned_sat_file, 'wb'))
    else: 
        sat_cen = pickle.load(open(cen_assigned_sat_file, 'rb'))

    abcrun_flag = cen_abcrun + '_central'
    theta_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'Satellite.tQdelay.theta_t', str(pewl), '_', abcrun_flag, '.dat']) 
   
    theta = np.loadtxt(theta_file(tf)) 
    med_theta = [np.median(theta[:,i]) for i in range(len(theta[0]))]

    abc_tqdelay_dict = {'name': 'explin', 'm': med_theta[0], 'b': med_theta[1]}
    sg_obj = EvolveSatSFR(sat_cen, tqdelay_dict=abc_tqdelay_dict)

    #infall = np.arange(len(sg_obj.first_infall_mass)) 

    infall = np.where(
            (sg_obj.first_infall_mass > 0.) &
            (sg_obj.first_infall_sfr > -999.))# & (sg_obj.first_infall_t >  5.7))
    print len(sg_obj.first_infall_mass)
    print len(infall[0]) 
    
    # SSFR plot
    ssfr_plot = PlotSSFR() 
    ssfr_plot.plot(mass=sg_obj.mass[infall], ssfr=sg_obj.ssfr[infall], line_color=3)

    ssfr_plot.GroupCat(position='satellite') 
    ssfr_plot.set_axes() 

    ssfr_fig_name = ''.join(['figure/', 
        'SSFR.Satellite', 
        '.tQdelayABC',
        '.', sg_obj.abcrun, '.', sg_obj.prior_name, '_prior', 
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
        '.tQdelayABC',
        '.', sg_obj.abcrun, '.', sg_obj.prior_name, '_prior', 
        '.png'])
    fq_plot.save_fig(fq_fig_name)
    plt.close() 

    return None 



def PlotABC_InfallCondition(tf, cen_tf=7, cen_abcrun=None, cen_prior_name=None): 
    ''' 
    '''
    # Simulator 
    cen_assigned_sat_file = ''.join(['/data1/hahn/pmc_abc/pickle/', 
        'satellite', '.cenassign', 
        '.', cen_abcrun, '_ABC', 
        '.', cen_prior_name, '_prior', '.p'])
    if not os.path.isfile(cen_assigned_sat_file): 
        sat_cen = AssignCenSFR(cen_tf, abcrun=cen_abcrun, prior_name=cen_prior_name)
        pickle.dump(sat_cen, open(cen_assigned_sat_file, 'wb'))
    else: 
        sat_cen = pickle.load(open(cen_assigned_sat_file, 'rb'))

    abcrun_flag = cen_abcrun + '_central'
    theta_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'Satellite.tQdelay.theta_t', str(pewl), '_', abcrun_flag, '.dat']) 
   
    theta = np.loadtxt(theta_file(tf)) 
    med_theta = [np.median(theta[:,i]) for i in range(len(theta[0]))]

    abc_tqdelay_dict = {'name': 'explin', 'm': med_theta[0], 'b': med_theta[1]}
    sg_obj = EvolveSatSFR(sat_cen, tqdelay_dict=abc_tqdelay_dict)

    #infall = np.where(
    #        (sg_obj.first_infall_mass > 0.) &
    #        (sg_obj.first_infall_sfr > -999.) & (sg_obj.first_infall_t >  5.7))
    infallmass0 = np.where(sg_obj.first_infall_mass == 0.) 
    infallmassG0 = np.where(sg_obj.first_infall_mass > 0.) 

    infall_nosfr = np.where(sg_obj.first_infall_sfr == -999.) 
    infall_longago = np.where(sg_obj.first_infall_t < 5.7) 

    print sg_obj.mass[infall_longago].min(), sg_obj.mass[infall_longago].max()
    print sg_obj.sfr[infall_longago].min(), sg_obj.sfr[infall_longago].max()
    print sg_obj.first_infall_sfr[infall_longago].min(), sg_obj.first_infall_sfr[infall_longago].max()

    
    m_arr = np.arange(9.5, 11.6, 0.1) 
    m_mid = 0.5 * (m_arr[:-1] + m_arr[1:])
    
    n_all, dum = np.histogram(sg_obj.mass, bins=m_arr)
    n_mass0, dum = np.histogram(sg_obj.mass[infallmass0], bins=m_arr)
    n_nosfr, dum = np.histogram(sg_obj.mass[infall_nosfr], bins=m_arr)
    n_longago, dum = np.histogram(sg_obj.mass[infall_longago], bins=m_arr)

    plt.plot(m_mid, n_mass0.astype('float')/n_all.astype('float')) 
    plt.plot(m_mid, n_nosfr.astype('float')/n_all.astype('float')) 
    plt.plot(m_mid, n_longago.astype('float')/n_all.astype('float')) 
    plt.show() 
    plt.close()

    ssfr_plot = PlotSSFR() 
    ssfr_plot.plot(mass=sg_obj.mass[infallmassG0], ssfr=sg_obj.ssfr[infallmassG0], line_color=5)
    ssfr_plot.plot(mass=sg_obj.mass, ssfr=sg_obj.ssfr, line_color='k', line_style='--')
    ssfr_plot.set_axes() 

    plt.show()
    return None 



if __name__=='__main__': 
    #for scat in [0.0, 0.2]: 
    #    subh = SatelliteSubhalos() 
    #    subh.build_catalogs(scatter=scat, source='li-march')
    #    del subh
    #ABC(20, [100.], Npart=200, cen_tf=6, cen_prior_name='updated', cen_abcrun='multifq_wideprior_nosmfevo')
    #PlotABC_tQdelay(12, cen_abcrun='multifq_wideprior_nosmfevo')
    #PlotABC_InfallCondition(10, cen_abcrun='multifq_wideprior', cen_prior_name='updated')

    #PlotABC_EvolvedSat(12, cen_abcrun='multifq_wideprior_nosmfevo', cen_prior_name='updated')
    #Plot_tRapid()
    #PlotABC_Corner(12, cen_abcrun='multifq_wideprior_nosmfevo')

    cen_assigned_sat_file = ''.join(['/data1/hahn/pmc_abc/pickle/', 
        'satellite.cenassign.multifq_wideprior_ABC.updated_prior.p'])
    sat_test_tmp = AssignCenSFR(7, abcrun='multifq_wideprior', prior_name='updated')
    pickle.dump(sat_test_tmp, open(cen_assigned_sat_file, 'wb'))
    #sat_cen = pickle.load(open(cen_assigned_sat_file, 'rb'))
    #evol_sat = EvolveSatSFR(sat_cen, tqdelay_dict={'name': 'hacked'})
    #PlotEvolvedSat(evol_sat)

    #sat_test_tmp = AssignCenSFR(7, abcrun='multifq_wideprior', prior_name='updated')
    #pickle.dump(sat_test_tmp, open('/data1/hahn/pmc_abc/pickle/sat_test_tmp.p', 'wb'))
    #sat_test_tmp = pickle.load(open('/data1/hahn/pmc_abc/pickle/sat_test_tmp.p', 'rb'))
    #evol_sat = EvolveSatSFR(sat_test_tmp)
    #PlotEvolvedSat(evol_sat)
