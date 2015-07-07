'''

Central_Quenchign project


Author(s): ChangHoon Hahn

'''

import numpy as np
import random
import h5py

#---- Local ----
import cenque_utility as util
import sf_mainseq as sfms

class CenQue: 
    ''' Central quenching (CenQue) galaxy catalogue
    '''
    def __init__(self): 
        self.mass = None
        self.sfr = None
        self.ssfr = None 
        self.parent = None
        self.child = None
        self.ilk = None
        self.snap_index = None
        self.pos = None     # position 
        self.gal_type = None    # quiescent/star-forming 
        self.halo_mass = None 
        self.sham_mass = None

        self.nsnap = None       # n_snapshot 
        self.zsnap = None           # z_snapshot
        self.t_cosmic = None    # t_cosmic for snapshot
        self.t_step = None      # t_cosmic step 

    def ImportSnap(self, nsnap): 
        ''' Import snapshot data from TreePM --> SHAM snapshots
        '''

        snapshot_file = ''.join(['dat/wetzel_tree/', 
            'subhalo_sham_centrals_snapshot', str(nsnap), '.hdf5']) 
        f = h5py.File(snapshot_file, 'r')       # read in h5py file 
        grp = f['cenque_data'] 

        self.mass = grp['mass'][:] 
        self.parent = grp['parent'][:]
        self.child = grp['child'][:]
        self.ilk = grp['ilk'][:]
        self.snap_index = grp['index'][:]
        self.pos = grp['pos'][:]
        self.halo_mass = grp['mass_halo'][:]
            
        # Old code that reads in fits rather than hdf5 files 
        '''
        snapshot_dir = 'dat/wetzel_tree/'
        snapshot_file = ''.join([ snapshot_dir, 
            'subhalo_sham_centrals_snapshot', str(nsnap), '.fits' 
            ]) 
        snapshot = util.mrdfits(snapshot_file)   # read snapshot file
        
        # import snapshot values
        self.mass = snapshot.mass 
        self.parent = snapshot.parent
        self.child = snapshot.child
        self.parent = snapshot.parent
        self.ilk = snapshot.ilk
        self.snap_index = snapshot.index
        self.pos = snapshot.pos
        self.halo_mass = snapshot.mass_halo 
        '''

        # get snapshot redshift/cosmic time data using Andrew's table
        n_snaps, z_snap, t_snap, t_wid = np.loadtxt('snapshot_table.dat', 
                unpack=True, usecols=[0, 2, 3, 4])

        self.nsnap = nsnap 
        self.zsnap = z_snap[(n_snaps.tolist()).index(nsnap)]        # redshift of snapshot
        self.t_cosmic = t_snap[(n_snaps.tolist()).index(nsnap)] # t_cosmic of snapshot 
        self.t_step = t_wid[(n_snaps.tolist()).index(nsnap)]    # Gyrs until next snapshot

    def AssignSFR(self, origin_nsnap, mass_bin=None, **kwargs): 
        ''' Assign SFR based on galaxy mass for snapshot 

        Assign galaxies as star-forming or quiescent
        based on quiescent fraction. 
        Then based on the assigned galaxy type 
        assign SSFR & SFR to the stellar masses
        
        Parameters
        ----------
        self : data class
        origin_nsnap : # of snapshot  
        mass_bin : mass bin for SFR assignment
        sfr : 'sfr_avg' or 'sfr_func'

        Notes
        -----
        * kwargs['sfr_avg'] : SFR assignment sampled about the SF-MS average SFR with 
        randomly sampled SFR residual 
            
        * kwargs['sfr_func'] : SFR assignment sampled about the SF-MS average SFR and an 
            arbitrary sampling of amplitude, phase, and frequency 

        '''
        # if AssignSFR is called 
        if self.mass is None:           # snapshot not yet imported
            # import snapshot 
            self.ImportSnap(origin_nsnap) 
        else: 
            pass 

        # get mass bins
        if mass_bin is None:     # if mass bin isn't specified 
            mass_bins = util.simple_mass_bin()  # use simplest mass bins
        else: 
            raise NameError("not yet coded") 

        # if slope/yint of the SF-MS is specified in kwargs
        if 'sfms_slope' in kwargs.keys(): 
            sfms_sfr_fit = sfms.get_sfmsfit_sfr(kwargs['sfms_slope'], kwargs['sfms_yint'])
        else:   
            pass
            ## otherwise, just use the fit from the group catalogs (hardcoded)
            #groupcat_slope, groupcat_yint = sfms.get_bestfit_groupcat_sfms(Mrcut=18)
            #sfms_sfr_fit = sfms.get_sfmsfit_sfr(groupcat_slope, groupcat_yint)
        
        if (min(self.mass) < min(mass_bins.mass_low)) or (max(self.mass) > max(mass_bins.mass_high)): 
            # remove galaxies below the minimum and maximum mass
            # and galaxies without children 
            mass_child_limit = (self.mass > min(mass_bins.mass_low)) & \
                    (self.mass <= max(mass_bins.mass_high)) & \
                    (self.child >= 0)  
    
            trim_column = [] 
            for col in self.__dict__.keys(): 
                col_data = getattr(self, col)

                if isinstance( col_data, type(self.mass) ):  
                    if len(col_data) == len(self.mass): 
                        trim_column.append(col) 
                    else: 
                        pass

            # only keep galaxies within the min and max mass
            self.sample_select(mass_child_limit, 
                    columns = trim_column)          

        if self.gal_type is None:         
            self.gal_type = np.array(['' for i in range(len(self.mass))], dtype='|S16') 
            self.sfr = np.array([-999. for i in range(len(self.mass))]) 
            self.ssfr = np.array([-999. for i in range(len(self.mass))]) 
    
            # SFR parameters
            if kwargs['sfr'] == 'sfr_func': 
                # periodic SFR 
                self.sfr_amp = np.array([-999. for i in range(len(self.mass))])
                self.sfr_freq = np.array([-999. for i in range(len(self.mass))])
                self.sfr_phase = np.array([-999. for i in range(len(self.mass))])

            elif kwargs['sfr'] == 'sfr_avg': 
                self.sfr_resid = np.array([-999. for i in range(len(self.mass))])
                
            else: 
                raise NotImplementedError('asdlkfj') 

        # loop through mass bins and assign SFRs ---------------------------------------------
        for i_m in range(mass_bins.nbins):              
        
            mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
            mass_bin_high = round(mass_bins.mass_high[i_m], 2) 
            # boolean list for mass range
            mass_bin_bool = (self.mass > mass_bin_low) & \
                    (self.mass <= mass_bin_high) &\
                    (self.gal_type == '')           # only assign for galaxies without galtype

            # indices of galaxies within mass range
            mass_bin_index = np.array(range(len(self.mass)))[mass_bin_bool] 
        
            mass_bin_ngal = np.sum(mass_bin_bool)   # Ngal in mass bin 
            if mass_bin_ngal == 0:                  # if there are no galaxies within mass bin
                continue 
    
            # get quiescent fraction for mass bin at z_snapshot 
            #mass_bin_qf = util.get_fquenching(mass_bins.mass_mid[i_m], self.zsnap, 
            #        yint=kwargs['fqing_yint'])
            mass_bin_qf = util.get_fq(mass_bins.mass_mid[i_m], self.zsnap, lit='wetzelsmooth') 

            # Ngal,quiescent in mass bin
            mass_bin_n_q = int( np.rint(mass_bin_qf * np.float(mass_bin_ngal)) ) 

            # Ngal,active in mass bin 
            mass_bin_n_sf = mass_bin_ngal - mass_bin_n_q

            #print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]
            #print 'fQ = ', mass_bin_qf, ' Ngal = ', mass_bin_ngal, \
            #        ' Ngal,Q = ', mass_bin_n_q, ' Ngal,SF = ', mass_bin_n_sf
            
            #### Quiescent ####
            if mass_bin_n_q > 0: 
                # randomly sample mass_bin_n_q galaxies from the mass bin 
                try: 
                    mass_bin_q_index = random.sample(mass_bin_index, mass_bin_n_q) 
                except ValueError: 
                    mass_bin_q_index = mass_bin_index 

                self.gal_type[mass_bin_q_index] = 'quiescent'   # label galaxy type 

                '''
                sample SSFR from log-normal distribution with 0.2dex scatter 
                centered about some determined mean
                '''
                q_ssfr_mean = util.get_q_ssfr_mean(self.mass[mass_bin_q_index]) 
                self.ssfr[mass_bin_q_index] = 0.18 * np.random.randn(mass_bin_n_q) + q_ssfr_mean
                # calculate SFR by SSFR+Mass
                self.sfr[mass_bin_q_index] = self.ssfr[mass_bin_q_index] + self.mass[mass_bin_q_index]

            #### Star-Forming ####
            if mass_bin_n_sf > 0: 
                try: 
                    mass_bin_sf_index = [x for x in mass_bin_index if x not in mass_bin_q_index]
                except NameError:       # if there aren't any quiescent galaxies
                    mass_bin_sf_index = mass_bin_index
                
                self.gal_type[mass_bin_sf_index] = 'star-forming'   # label galaxy type 
            
                if kwargs['sfr'] == 'sfr_avg': 
                    #get average and scatter of SF main sequence 
                    [sf_avg_sfr, sf_sig_sfr] = util.get_sfr_mstar_z_bestfit(
                            mass_bins.mass_mid[i_m], self.zsnap, Mrcut=18) 

                    #print 'SF Average(SFR) = ', sf_avg_sfr, ' sigma_SFR = ', sf_sig_sfr

                    self.sfr_resid[mass_bin_sf_index] = sf_sig_sfr * np.random.randn(mass_bin_n_sf)
                    self.sfr[mass_bin_sf_index] = util.sfr_avg_residual(
                            mass_bins.mass_mid[i_m], self.zsnap, 
                            resid = self.sfr_resid[mass_bin_sf_index]) 
                    #self.sfr[mass_bin_sf_index] = sf_sig_sfr * np.random.randn(mass_bin_n_sf) + sf_avg_sfr 

                elif kwargs['sfr'] == 'sfr_func': 
                    # Fluctuating SFR with random amplitude, frequency, and phase 

                    # average and 1sigma of SF-MS 
                    [sf_avg_sfr, sf_sig_sfr] = util.get_sfr_mstar_z_bestfit(
                            mass_bins.mass_mid[i_m], self.zsnap, Mrcut=18) 

                    # sfr amplitude, frequency and phase respectively 
                    self.sfr_amp[mass_bin_sf_index] = \
                            sf_sig_sfr * np.random.randn(mass_bin_n_sf) 
                    self.sfr_freq[mass_bin_sf_index] = \
                            (2.0 * np.pi)/np.random.uniform(0.01, 0.1, mass_bin_n_sf)  
                    self.sfr_phase[mass_bin_sf_index] = \
                            np.random.uniform(0.0, 1.0, mass_bin_n_sf)

                    self.sfr[mass_bin_sf_index] = util.sfr_squarewave(
                            self.mass[mass_bin_sf_index], self.t_cosmic, 
                            amp = self.sfr_amp[mass_bin_sf_index], 
                            freq = self.sfr_freq[mass_bin_sf_index], 
                            phase = self.sfr_phase[mass_bin_sf_index])
                else: 
                    raise NotImplementedError('asdfkjlkjasdf') 

                self.ssfr[mass_bin_sf_index] = \
                        self.sfr[mass_bin_sf_index] - self.mass[mass_bin_sf_index]

        # check for SFR/SSFR assign fails
        assign_fail = (self.sfr == -999.0) | (self.ssfr == -999.0)
        if np.sum(assign_fail) > 0: 
            raise NameError('asdfasdfasdfasdf')

    def readin(self, **kwargs): 
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

    def writeout(self, **kwargs): 
        ''' simply outputs specified columns to file
        '''
        output_file = util.cenque_file( **kwargs )  # output file 
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

    def sample_select(self, bool, columns=None, silent=False): 
        ''' Given boolean list remove False items from data
        '''

        if not silent: 
            n_remove = len(bool) - len(bool[bool == True])
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
    
    def get_index(self, bool): 
        ''' Given boolean from conditional statement, give indices of boolean
        '''
        return np.array(range(len(self.mass)))[bool]

def EvolveCenQue(origin_nsnap, final_nsnap, mass_bin=None, silent=True, **kwargs): 
    ''' Evolve SF properties from origin_nsnap to final_nsnap 
    
    Parameters
    ----------

    Notes
    -----
    * Quenching of SF galaxies now involve quenching-efold *and* overall SF MS SFR decrease

    '''

    if mass_bin is None:            # mass bins
        mass_bins = util.simple_mass_bin()  # use simplest mass bins
    else: 
        raise NotImplementedError("not yet coded") 

    # SF-MS fits 
    if 'sfms_slope' in kwargs.keys(): 
        sfms_sfr_fit = sfms.get_sfmsfit_sfr(kwargs['sfms_slope'], kwargs['sfms_yint'])
    else: 
        pass
        #groupcat_slope, groupcat_yint = sfms.get_bestfit_groupcat_sfms(Mrcut=18)
        #sfms_sfr_fit = sfms.get_sfmsfit_sfr(groupcat_slope, groupcat_yint)
  
    # import original snap SF prop 
    parent_cq = CenQue()
    parent_cq.readin(nsnap=origin_nsnap, file_type='sf assign', **kwargs)   
    
    if not silent: 
        print 'Quiescent Fraction = ', np.float(len(parent_cq.gal_type[parent_cq.gal_type == 'quiescent']))/np.float(len(parent_cq.gal_type)) 
    
    # evolve snapshot by snapshot ------------------------------------------------------
    for i_step in range(0, origin_nsnap - final_nsnap):    

        i_snap = origin_nsnap - i_step  # current Nsnap 
        child_snap = i_snap-1

        child_cq = CenQue() 
        child_cq.readin(nsnap=child_snap)         # read in snapshot SHAM masses
        
        # remove galaxies beyond mass bins
        child_mass_limit = (child_cq.mass > min(mass_bins.mass_low)) & \
                (child_cq.mass <= max(mass_bins.mass_high))
        child_cq.sample_select( child_mass_limit, 
                columns = ['mass', 'sham_mass', 'halo_mass', 'parent', 'child', 'ilk', 'snap_index']  
                )                   
        n_child = len(child_cq.mass)     # number of children left 

        # set up columns for assignment 
        child_cq.gal_type = np.array(['' for i in range(n_child)], dtype='|S16') 
        child_cq.sfr = np.array([-999. for i in range(n_child)]) 
        child_cq.ssfr = np.array([-999. for i in range(n_child)]) 
        child_cq.q_ssfr = np.array([-999. for i in range(n_child)]) 
        child_cq.tau = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_sfr = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_mass = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_halo_mass = np.array([-999. for i in range(n_child)]) 

        if kwargs['sfr'] == 'sfr_func': 
            child_cq.sfr_amp = np.array([-999. for i in range(n_child)]) 
            child_cq.sfr_freq = np.array([-999. for i in range(n_child)]) 
            child_cq.sfr_phase = np.array([-999. for i in range(n_child)]) 

        elif kwargs['sfr'] == 'sfr_avg': 
            child_cq.sfr_resid = np.array([-999. for i in range(n_child)]) 

        else: 
            raise NotImplementedError('asdlkfj') 

        # parent children match
        # assign indices into dictionaries and then get dictionary values
        # not very good -----------------------------------------------------------------
        parent_index_dict = dict( (k, i) for i, k in enumerate(parent_cq.snap_index) )
        child_index_dict = dict( (k, i) for i, k in enumerate(child_cq.parent) )

        parent_child_intersect = np.intersect1d( parent_cq.snap_index, child_cq.parent ) 
        parent_indx = [ parent_index_dict[x] for x in parent_child_intersect ] 
        child_indx = [ child_index_dict[x] for x in parent_child_intersect ]  
        # --------------------------------------------------------------------------------
        
        # children inherit galaxy type and store parent SFR and mass 
        (child_cq.gal_type)[child_indx] = [(parent_cq.gal_type)[i] for i in parent_indx]
        (child_cq.parent_sfr)[child_indx] = [(parent_cq.sfr)[i] for i in parent_indx]
        (child_cq.parent_mass)[child_indx] = [(parent_cq.mass)[i] for i in parent_indx]
        (child_cq.parent_halo_mass)[child_indx] = [(parent_cq.halo_mass)[i] for i in parent_indx]
       
        if kwargs['sfr'] == 'sfr_func': 
            (child_cq.sfr_amp)[child_indx] = [(parent_cq.sfr_amp)[i] for i in parent_indx]
            (child_cq.sfr_freq)[child_indx] = [(parent_cq.sfr_freq)[i] for i in parent_indx]
            (child_cq.sfr_phase)[child_indx] = [(parent_cq.sfr_phase)[i] for i in parent_indx]

        elif kwargs['sfr'] == 'sfr_avg': 
            (child_cq.sfr_resid)[child_indx] = [(parent_cq.sfr_resid)[i] for i in parent_indx]

        else: 
            raise NotImplementedError('lkajsdfkj')
        
        # Star-forming Children ----------------------------------------
        sf_child = (child_cq.gal_type[child_indx] == 'star-forming')
        sf_child_indx = np.array(child_indx)[sf_child]                  # index stuff
        sf_child_parent_indx = np.array(parent_indx)[sf_child]
        
        # SFR evolution amount
        #child_sfr, child_sig_sfr = util.get_sfr_mstar_z(child_cq.mass[sf_child_indx], 
        #        child_cq.zsnap, lit='primusfit')
        #parent_sfr, parent_sig_sfr = util.get_sfr_mstar_z( parent_cq.mass[sf_child_parent_indx], 
        #        parent_cq.zsnap, lit='primusfit') 
        #child_sfr, child_sig_sfr = util.get_sfr_mstar_z_flex(child_cq.mass[sf_child_indx], 
        #        child_cq.zsnap, sfms_sfr_fit)
        #parent_sfr, parent_sig_sfr = util.get_sfr_mstar_z_flex(
        #        parent_cq.mass[sf_child_parent_indx], parent_cq.zsnap, sfms_sfr_fit) 
    
        if kwargs['sfr'] == 'sfr_avg': 
            child_sfr, child_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    child_cq.mass[sf_child_indx], child_cq.zsnap, Mrcut=18)
            parent_sfr, parent_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    parent_cq.mass[sf_child_parent_indx], parent_cq.zsnap, Mrcut=18) 
        
            ''' SFR evolution is assumed to be equal to the overall change in SFR  
            '''
            dSFRt = child_sfr - parent_sfr          # overall change in SFR
            
            # evolve sf children SFR
            #child_cq.sfr[sf_child_indx] = child_cq.parent_sfr[sf_child_indx] + dSFRt
            child_cq.sfr[sf_child_indx] = util.sfr_avg_residual(
                    child_cq.mass[sf_child_indx], child_cq.zsnap, 
                    resid = child_cq.sfr_resid[sf_child_indx] ) 
        
        elif kwargs['sfr'] == 'sfr_func': 

            child_cq.sfr[sf_child_indx] = util.sfr_squarewave(
                    child_cq.mass[sf_child_indx], child_cq.t_cosmic, 
                    amp = child_cq.sfr_amp[sf_child_indx], 
                    freq = child_cq.sfr_freq[sf_child_indx], 
                    phase = child_cq.sfr_phase[sf_child_indx]) 

        else: 
            raise NotImplementedError('asdflkjadf') 

        child_cq.ssfr[sf_child_indx] = child_cq.sfr[sf_child_indx] - child_cq.mass[sf_child_indx]
        
        # --------------------------------------------------------------------------------
        # mass evolution of star-forming galaxies (only star-forming galaxies) 
        # should quiescent galaxies remain the same mass or adopt SHAM masses? 

        if kwargs['stellmass'].lower() == 'integrated':

            # integrated stellar mass
            if kwargs['sfr'] == 'sfr_func': 
                integrated_mass = util.integrated_mass_rk4(
                        util.sfr_squarewave, (child_cq.parent_mass)[sf_child_indx], 
                        parent_cq.t_cosmic, child_cq.t_cosmic, 
                        amp = (child_cq.sfr_amp)[sf_child_indx], 
                        freq = (child_cq.sfr_freq)[sf_child_indx], 
                        phase = (child_cq.sfr_phase)[sf_child_indx])
                (child_cq.mass)[sf_child_indx] = integrated_mass

            elif kwargs['sfr'] == 'sfr_avg': 

                (child_cq.mass)[sf_child_indx] = util.integrated_mass_rk4(
                        util.sfr_avg_residual, (child_cq.parent_mass)[sf_child_indx], 
                        parent_cq.t_cosmic, child_cq.t_cosmic, 
                        resid = (child_cq.sfr_resid)[sf_child_indx])
                if not silent: 
                    print 'Integrated Mass vs SHAM mass' 
                    print (child_cq.mass)[sf_child_indx] - (child_cq.sham_mass)[sf_child_indx]

        elif kwargs['stellmass'].lower() == 'sham': 
            # leave stellar mass 
            pass
        else: 
            raise NotImplementedError('asdfasdflkjasdflkj;') 

        # Quenching Children ------------------------------------------------        
        try: 
            parent_cq.tau
        except AttributeError:
            pass
        else: 
            # children inherit parent tau 
            (child_cq.tau)[child_indx] = [(parent_cq.tau)[i] for i in parent_indx]
            has_tau = (child_cq.tau)[child_indx] > 0.0  # has tau 
            still_quenching_indx = np.array(child_indx)[has_tau]
            still_quenching_parent_indx = np.array(parent_indx)[has_tau]
            
            # for galaxies with tau, keep them quenching!
            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[still_quenching_indx]
                ))                                              # SFR quenching amount  

            c_sfr, c_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    child_cq.mass[still_quenching_indx], child_cq.zsnap, Mrcut=18)
            p_sfr, p_sig_sfr = util.get_sfr_mstar_z_bestfit(
                    child_cq.parent_mass[still_quenching_indx], parent_cq.zsnap, Mrcut=18) 

            child_cq.sfr[still_quenching_indx] = \
                    child_cq.parent_sfr[still_quenching_indx] + tau_quench 
            child_cq.ssfr[still_quenching_indx] = \
                    child_cq.sfr[still_quenching_indx] - child_cq.mass[still_quenching_indx]
            
            # children inherit final quenched SSFR 
            (child_cq.q_ssfr)[child_indx] = [(parent_cq.q_ssfr)[i] for i in parent_indx]

        # Quiescent Children --------------------------------------------
        q_child = (child_cq.gal_type[child_indx] == 'quiescent') & (child_cq.tau[child_indx] == -999.0) # not quenching
        q_child_indx = np.array(child_indx)[q_child]
        q_child_parent_indx = np.array(parent_indx)[q_child]

        # keep SSFR same 
        child_cq.ssfr[q_child_indx] = parent_cq.ssfr[q_child_parent_indx]
        child_cq.sfr[q_child_indx] = child_cq.ssfr[q_child_indx] + \
                child_cq.mass[q_child_indx]
        
        quenching_fractions = [] 
        quenching_fractionss = [] 

        for i_m in range(mass_bins.nbins):              
            # boolean list for mass range
            mass_bin_bool = (child_cq.mass > mass_bins.mass_low[i_m]) & \
                    (child_cq.mass <= mass_bins.mass_high[i_m]) & \
                    (child_cq.gal_type != '') 
            if not silent: 
                print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]

            # indices of galaxies within mass range
            mass_bin_index = child_cq.get_index(mass_bin_bool) 
        
            mbin_ngal = np.sum(mass_bin_bool)   # Ngal in mass bin 

            if mbin_ngal == 0:              
                # if there are no galaxies within mass bin
                print 'No Galaxies in ', mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]
                continue 
            
            # observed quiescent fraction at M* and z
            mbin_qf = util.get_fq(mass_bins.mass_mid[i_m], child_cq.zsnap, lit='wetzelsmooth') 
            
            if not silent: 
                print 'nsnap = ', child_cq.nsnap, ' z = ', child_cq.zsnap, ' M* = ', mass_bins.mass_mid[i_m], ' fq = ', mbin_qf
                
            # Number of expected quiescent galaxies (Ngal,Q_exp)
            mbin_exp_n_q = int( np.rint(mbin_qf * np.float(mbin_ngal)) )    

            if mbin_ngal < mbin_exp_n_q: 
                # if for some reason expected number of quiescent galaxies
                # exceed the number of galaxies in the bin, all of them are 
                # quiescent
                mbin_exp_n_q = mbin_ngal 
            
            child_gal_type = child_cq.gal_type[mass_bin_bool] 
            mbin_sf_ngal = np.sum(child_gal_type == 'star-forming') 
            mbin_q_ngal = np.sum(child_gal_type == 'quiescent')  

            # number of SF galaxies that need to be quenched 
            ngal_2quench = mbin_exp_n_q - mbin_q_ngal 

            if not silent: 
                print 'sf', mbin_sf_ngal, 'q', mbin_q_ngal
                print mbin_exp_n_q, ' expected quiescent galaxies' 
                print 'current fq = ', np.float(mbin_q_ngal)/np.float(mbin_ngal)
                print ngal_2quench, ' SF galaxies need to be quenched'

            if ngal_2quench <= 0: 
                # quenching, overshot in previous snapshot  
                # move onto next mass bin 
                if not silent: 
                    print 'Quenching overshot in previous snapshot' 
                pass
                #continue 

            if mbin_sf_ngal == 0: 
                continue
            
            # Star-forming children 
            mbin_sf_index = child_cq.get_index( 
                    [ mass_bin_bool & (child_cq.gal_type == 'star-forming')] 
                    ) 
                
            # number of SF galaxies to quench will be determined by first determining how many galaxies will be quenched
            # if the entire SF population was quenched, taking the ratio of that number over the 'mbin_exp_n_q'. 
            # we use that factor as the "quenching fraction" (the fraction of SF galaxies that will be quenched for this mass
            # bin and redshift) 
            if 'tau' in kwargs.keys(): 
                if 'tau_param' in kwargs.keys(): 
                    taus = util.get_quenching_efold(
                            child_cq.parent_mass[mbin_sf_index], 
                            type=kwargs['tau'], param=kwargs['tau_param'])
                else: 
                    taus = util.get_quenching_efold(
                            child_cq.parent_mass[mbin_sf_index], 
                            type=kwargs['tau'])
            else: 
                raise TypeError('specify quenching e-fold: tau = instant, constant, linear') 
            
            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / taus 
                ))      
            #print child_cq.parent_sfr[mbin_sf_index] + tau_quench - child_cq.mass[mbin_sf_index]
            #print tau_quench
            #print min(child_cq.sfr[mbin_sf_index] + tau_quench - child_cq.mass[mbin_sf_index])
            #print max(child_cq.sfr[mbin_sf_index] + tau_quench - child_cq.mass[mbin_sf_index])

            sfqs = util.sfq_classify(child_cq.mass[mbin_sf_index], child_cq.sfr[mbin_sf_index] + tau_quench, child_cq.zsnap)
            ngal_totalq = np.float(np.sum(sfqs == 'quiescent'))
            #print min(child_cq.sfr[mbin_sf_index[sfqs == 'quiescent']] + tau_quench[sfqs == 'quiescent'] - child_cq.mass[mbin_sf_index[sfqs == 'quiescent']])
            #print max(child_cq.sfr[mbin_sf_index[sfqs == 'quiescent']] + tau_quench[sfqs == 'quiescent'] - child_cq.mass[mbin_sf_index[sfqs == 'quiescent']])
    
            if ngal_totalq == 0: 
                raise NameError("What the fuck") 
            #print 'ngal_totalq', ngal_totalq

            alpha = 1.5
            quenching_fraction = 0.025 * (mass_bins.mass_mid[i_m] - 9.5) * (1.8 - child_cq.zsnap)**alpha + 0.15
            quenching_fractions.append(quenching_fraction)

            quenching_fraction = np.float(ngal_2quench)/np.float(ngal_totalq)
            quenching_fractionss.append(quenching_fraction)
            
            alpha = 1.5
            fqing1 = 0.025 * (mass_bins.mass_mid[i_m] - 9.5) * (1.8 - child_cq.zsnap)**alpha + 0.15
                
            fqing2 = np.float(ngal_2quench)/np.float(ngal_totalq)

            if (mass_bins.mass_mid[i_m] < 11.0) and (fqing2 > fqing1): 
                quenching_fraction = fqing1
                if not silent: 
                    print '####################################### FQING 2 > FQING 1'
            else: 
                quenching_fraction = fqing2

            if quenching_fraction < 0.0: 
                quenching_fraction = 0.0
            elif quenching_fraction > 1.0: 
                quenching_fraction = 1.0
            #print 'quenching fraction ', quenching_fraction
            if quenching_fraction == 0.0: 
                continue

            #if quenching_fraction > 1.0: 
            #    #raise NameError('asdfasdfasdf')
            #    #quench_index = mbin_sf_index
            #    print 'QUENCHING FRACTION TOO LARGE ******************************************************************'
            #    quench_index = random.sample(mbin_sf_index, ngal_2quench) 
            #else: 
            quench_index = random.sample(mbin_sf_index, int( np.rint(quenching_fraction * np.float(mbin_sf_ngal) ) ) )
            
            #quench_index = random.sample(mbin_sf_index, ngal_2quench) 
            
            child_cq.gal_type[quench_index] = 'quiescent'  # boom quenched 
            
            # assign quenching e-folds for quenching galaxies
            if 'tau' in kwargs.keys(): 
                if 'tau_param' in kwargs.keys(): 
                    child_cq.tau[quench_index] = util.get_quenching_efold(
                            child_cq.parent_mass[quench_index], 
                            type=kwargs['tau'], param=kwargs['tau_param'])
                else: 
                    child_cq.tau[quench_index] = util.get_quenching_efold(
                            child_cq.parent_mass[quench_index], 
                            type=kwargs['tau'])
            else: 
                raise TypeError('specify quenching e-fold: tau = instant, constant, linear') 
    
            # SFR quenching amount
            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[quench_index]
                ))      

            child_cq.sfr[quench_index] = child_cq.parent_sfr[quench_index] + tau_quench 
            child_cq.ssfr[quench_index] = child_cq.sfr[quench_index] - child_cq.mass[quench_index]

            q_ssfr_mean = util.get_q_ssfr_mean(child_cq.mass[quench_index]) 
            child_cq.q_ssfr[quench_index] = 0.18 * np.random.randn(len(quench_index)) + q_ssfr_mean 
        
        if not silent: 
            print child_cq.zsnap
            print mass_bins.mass_mid
            print quenching_fractions
            print quenching_fractionss

        # deal with orphans ------------------------------------------------------------

        sfing = np.where(child_cq.gal_type == 'star-forming') 
        print max(child_cq.mass[sfing]) 
        print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
        child_cq.AssignSFR(child_cq.nsnap, **kwargs) 

        if len(child_cq.gal_type[child_cq.gal_type == '']) != 0:
            print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
            raise NameError('asdflkjasdlfj') 
        
        # deal with galaxies that are being quenched beyond the designated final quenching  
        over_quenched = child_cq.ssfr < child_cq.q_ssfr 
        child_cq.ssfr[over_quenched] = child_cq.q_ssfr[over_quenched]
        child_cq.tau[over_quenched] = -999.0            # done quenching 

        #if child_cq.nsnap == final_nsnap: 
    
        if kwargs['stellmass'] == 'sham': 
            if kwargs['sfr'] == 'sfr_avg': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_resid', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                        **kwargs)  
            elif kwargs['sfr'] == 'sfr_func': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_amp', 'sfr_freq', 'sfr_phase', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                    **kwargs)  

        elif kwargs['stellmass'] == 'integrated': 
            if kwargs['sfr'] == 'sfr_avg': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_resid', 'sham_mass', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                        **kwargs)  
            elif kwargs['sfr'] == 'sfr_func': 
                child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap),
                        columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                            'halo_mass', 'sfr_amp', 'sfr_freq', 'sfr_phase', 'sham_mass', 
                            'parent_sfr', 'parent_mass', 'parent_halo_mass', 'parent', 
                            'child', 'ilk', 'snap_index'], 
                    **kwargs)  
        else: 
            raise NotImplementedError('asdflkjasdkf') 

        parent_cq = child_cq

        if not silent: 
            print 'Quiescent Fraction = ', np.float(len(parent_cq.gal_type[parent_cq.gal_type == 'quiescent']))/np.float(len(parent_cq.gal_type)) 

def build_cenque_importsnap(**kwargs): 
    ''' 
    '''
    for i_snap in [13]: 
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
    #build_cenque_original(sfr='sfr_avg') 
    EvolveCenQue(13, 1, tau='linefit', tau_param=[-0.7, 0.4], 
            sfr='sfr_avg', stellmass='sham') 
