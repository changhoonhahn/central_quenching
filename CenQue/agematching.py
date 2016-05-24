'''

Test age matching 


'''
import h5py
import pickle 
import numpy as np
import os.path
from scipy.interpolate import interp1d

# --- Local --- 
import util.util as Util
from observations import GroupCat
from satellite import SGPop 
from satellite import EvolveSatSFR

import sham_hack as sham
import subhalo_io_hack as subhalo_io
from utilities import utility as wetzel_util



class AllSubhalos(object): 
    def __init__(self, scatter=0.2, source='li-march'): 
        ''' Class that describes the All Subhalo Catalogs generated from 
        Andrew Wetzel's TreePM merger tree code. 
        '''
        self.scatter = scatter
        self.source = source 
        self.ReadSHAM()

    def Build(self, age_indicator=None):
        '''
        '''
        self.CentralSFR()
        self.SatelliteSFR()

        if age_indicator is not None: 
            self.CentralAgeMatch(age_indicator=age_indicator)

        self.WriteASCII(age_indicator=age_indicator)

    def CentralSFR(self):
        ''' Assign SFRs to central galaxies based on the SDSS DR7 observed central 
        SFRs.
        '''
        # central galaxy indices
        cen_index = np.where( (self.ilk >= 1.0) & (self.ilk < 2.01) )
        self.cen_index = cen_index
        cen_mass = getattr(self, 'm.star')[cen_index]
    
        # import central SDSS DR7 galaxies 
        sdss = GroupCat(Mrcut=18, position='central')
        sdss.Read()

        dlogm = 0.1
        m_arr = np.arange(sdss.mass.min()-0.5*dlogm, 11.2+dlogm, dlogm)
        m_low = m_arr[:-1]
        m_high = m_arr[1:]
        
        self.sfr = np.repeat(-999., len(getattr(self, 'm.star')))
        self.ssfr = np.repeat(-999., len(getattr(self, 'm.star')))
        self.position = np.repeat(-1, len(getattr(self, 'm.star'))) 

        self.position[cen_index] = 0 

        for im in range(len(m_low)):
            if im < len(m_low) - 1: 
                sdss_bin = np.where((sdss.mass >= m_low[im]) & (sdss.mass < m_high[im]))
                subh_bin = np.where((cen_mass >= m_low[im]) & (cen_mass < m_high[im]))
            else: 
                sdss_bin = np.where(sdss.mass >= m_low[im])
                subh_bin = np.where(cen_mass >= m_low[im])

            sdss_sfr = sdss.sfr[sdss_bin]
            dist, bin_edges = np.histogram(
                    sdss_sfr, 
                    range=[-3., 3.], 
                    bins=50, 
                    normed=True)

            sfr_cdf = np.zeros(len(bin_edges))
            sfr_cdf[1:] = np.cumsum(dist * np.diff(bin_edges))
            cdf_interp = interp1d(sfr_cdf, bin_edges, kind='linear') 
            
            rand_i = np.random.uniform(size=len(subh_bin[0]))
            cen_sfr = cdf_interp(rand_i)
            #print cen_sfr.min(), cen_sfr.max()
            self.sfr[cen_index[0][subh_bin]] = cen_sfr 
            self.ssfr[cen_index[0][subh_bin]] = self.sfr[cen_index[0][subh_bin]] - getattr(self, 'm.star')[cen_index[0][subh_bin]]

        return None

    def CentralAgeMatch(self, age_indicator=None): 
        '''
        '''
        try: 
            self.cen_index
        except AttributeError: 
            self.CentralSFR() 
            
        dlogm = 0.1
        m_arr = np.arange(9.7-0.5*dlogm, 11.2+dlogm, dlogm)
        m_low = m_arr[:-1]
        m_high = m_arr[1:]
        
        # central Mhalo evolution 
        mhalo_evol = self.Mhalo_evol[self.cen_index]
        mhalo_evol[:,0] = getattr(self, 'm.halo')[self.cen_index].copy() 
        cen_mhalo = getattr(self, 'm.halo')[self.cen_index]

        mhalo_ratio = 10**(mhalo_evol - cen_mhalo[:,None])
        
        # whatever metric for age 
        if age_indicator == 'acc80': 
            age = np.abs(
                    mhalo_ratio - np.tile(0.8, (mhalo_ratio.shape[0], mhalo_ratio.shape[1]))
                    ).argmin(axis=1) + 1
        self.central_age = age 
        cen_mass = getattr(self, 'm.star')[self.cen_index]
        cen_sfr = self.sfr[self.cen_index] 
        for im in range(len(m_low)):
            if im < len(m_low) - 1: 
                subh_bin = np.where((cen_mass >= m_low[im]) & (cen_mass < m_high[im]))[0]
            else: 
                subh_bin = np.where(cen_mass >= m_low[im])[0]
            
            # sort indices of halo ages that have M* within this mass bin 
            # low snapshot -> high snapshot
            isort_age = np.argsort(age[subh_bin])

            # sort SFRs of galaxies within this mass bin 
            sfr_sorted = np.sort(cen_sfr[subh_bin])

            # younger halos get higher sfr!
            self.sfr[self.cen_index[0][subh_bin[isort_age]]] = sfr_sorted[::-1]

        return None

    def SatelliteSFR(self, cen_abcrun='multifq_wideprior', cen_prior_name='updated'):
        ''' 
        '''
        # Satellites with central SFR properties assigned at infall
        cen_assigned_sat_file = ''.join([
            '/data1/hahn/pmc_abc/pickle/', 
            'satellite', 
            '.cenassign', 
            '.', cen_abcrun, '_ABC', 
            '.', cen_prior_name, '_prior', 
            '.p'])
        sat_cen = pickle.load(open(cen_assigned_sat_file, 'rb'))
        # Evolve the satellite SFRs based on hacked tqdelay and tq_rapid 
        # from Wetzel et al. (2013)
        evol_sat = EvolveSatSFR(sat_cen, tqdelay_dict={'name': 'hacked'})
        self.sat_index = evol_sat.snap_index
        self.sfr[evol_sat.snap_index] = evol_sat.sfr
        self.ssfr[evol_sat.snap_index] = evol_sat.sfr - evol_sat.mass 
        self.position[evol_sat.snap_index] = 1 
        
        masses = getattr(self, 'm.star')
        masses[evol_sat.snap_index] = evol_sat.mass
        setattr(self, 'm.star', masses)

        sat_mass = getattr(self, 'm.star')[evol_sat.snap_index]
        return None 

    def WriteASCII(self, age_indicator=None):
        '''
        '''
        if age_indicator is None: 
            age_str = '.NOagematch'
        else: 
            age_str = '.agematch_'+age_indicator

        output_file = ''.join([
            '/data1/hahn/agematch/'
            'subhalo_sham.all', 
            '.scatter', str(self.scatter), 
            '.source_', self.source, 
            age_str, 
            '.dat'])

        x = self.pos[:,0]
        y = self.pos[:,1]
        z = self.pos[:,2]
        vx = self.vel[:,0]
        vy = self.vel[:,1]
        vz = self.vel[:,2]
        mstar = getattr(self, 'm.star')
        masslim = np.where(mstar > 9.7)

        data_list = [x[masslim], y[masslim], z[masslim], vx[masslim], vy[masslim], vz[masslim], 
                mstar[masslim], self.sfr[masslim], self.ssfr[masslim], self.position[masslim]]
        data_hdrs = ' x, y, z, vx, vy, vz, mstar, sfr, ssfr, position (0=cen, 1=sat)'

        np.savetxt(
                output_file, 
                np.array(data_list).T, 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%i'],
                delimiter='\t', 
                header=data_hdrs
                ) 
        return None
    
    def ReadSHAM(self): 
        ''' Read in the SHAMed halo catalog 
        '''
        spec_kwargs = {
                'scatter': self.scatter, 
                'snapshot': 1, 
                'source' : self.source
                }
        snapshot_file = self._FileSHAM(**spec_kwargs)

        if not os.path.isfile(snapshot_file): 
            print snapshot_file, ' does NOT exist'
            print 'Now building'
            self.BuildSHAM()

        f = h5py.File(snapshot_file, 'r')
        grp = f['data']

        for i_meta, metadatum in enumerate(grp.attrs.keys()): 
            setattr(self, metadatum, (grp.attrs.values())[i_meta])
        self.metadata = [str(key) for key in grp.attrs.keys()]

        for datum in grp.keys():
            setattr(self, datum, grp[datum][:])
        self.data_columns = [str(key).encode('ascii') for key in grp.keys()]

        return None

    def BuildSHAM(self): 
        ''' Build central subhalo snapshot catalogs using TreePM code. The kwargs specify the 
        properties of the central subhalos
        '''
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
        # subhalo properties
        mhalo     = (sub[1]['halo.m'])      # M_halo
        mhalo_max = (sub[1]['m.max'])       # M_halo_max (Mpeak)
        mstar     = (sub[1]['m.star'])      # M* of subhalo
        pos       = (sub[1]['pos'])         # position of subhalo
        vel       = (sub[1]['vel'])         # velocity of subhalo
        ilk       = (sub[1]['ilk'])         # classification flag in case we want to consider centrals and virtual centrals separately 
        zi_infall  = (sub[1]['inf.first.zi']) # infall snapshot 

        spec_kwargs = {
                'scatter': self.scatter, 
                'snapshot': 1, 
                'source': self.source
                }
        subsham_file = self._FileSHAM(**spec_kwargs)    # output name

        f = h5py.File(subsham_file, 'w') 
        grp = f.create_group('data') 
        
        grp.attrs['snapshot'] = 1 
        grp.create_dataset('index', data=range(len(mhalo)))
        grp.create_dataset('pos', data=pos) 
        grp.create_dataset('vel', data=vel) 
        grp.create_dataset('ilk', data=ilk)
        grp.create_dataset('m.star', data=mstar)
        grp.create_dataset('m.halo', data=mhalo)
        grp.create_dataset('m.halo.max', data=mhalo_max)
        grp.create_dataset('inf.first.zi', data=zi_infall)

        Mstar_evol = np.tile(-999., (len(mhalo),32))
        Mhalo_evol = np.tile(-999., (len(mhalo),32))
        Mhalomax_evol = np.tile(-999., (len(mhalo),32))

        for i_snap in range(2, 33): 
            # ancestor index of subhalo for ancestor in snapshot nsnap_ancestor
            ancestor = wetzel_util.utility_catalog.indices_tree(
                    sub, 1, i_snap, range(len(mhalo)))
            grp.create_dataset('ancestor'+str(i_snap), data=ancestor)
            
            # ancestor SHAM stellar mass and halo masses  
            has_ancestor, has_descendant = Util.intersection_index(
                    ancestor, np.arange(len(sub[i_snap]['m.star'])))
            
            # keep track of halo mass  
            (Mstar_evol[:,i_snap-1])[has_ancestor] = sub[i_snap]['m.star'][has_descendant]
            (Mhalo_evol[:,i_snap-1])[has_ancestor] = sub[i_snap]['halo.m'][has_descendant]
            (Mhalomax_evol[:,i_snap-1])[has_ancestor] = sub[i_snap]['m.max'][has_descendant]
        
        grp.create_dataset('Mstar_evol', data=Mstar_evol)
        grp.create_dataset('Mhalo_evol', data=Mhalo_evol)
        grp.create_dataset('Mhalomax_evol', data=Mhalomax_evol)
        f.close() 
        return None

    def _FileSHAM(self, **spec_kwargs): 
        ''' File name of output file given snapshot 
        '''
        # write to hdf5 file 
        subsham_file = ''.join([ 
            '/data1/hahn/agematch/'
            'subhalo_sham', 
            '.all',
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




if __name__=='__main__':
    owl = AllSubhalos()
    #owl.BuildSHAM()
    owl.Build(age_indicator='acc80')
