'''

Documentation here 

Author(s): ChangHoon Hahn 

'''

import numpy as np
import random

#---- Local ----
import cenque_utility as util

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

    def ImportSnap(self, nsnap): 
        ''' Import snapshot data from TreePM --> SHAM snapshots
        '''
        snapshot_dir = '/data1/hahn/wetzel_tree/'
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
        self.pos = snapshot.index 

        # get snapshot redshift/cosmic time data using Andrew's table
        n_snaps, z_snap, t_snap, t_wid = np.loadtxt('snapshot_table.dat', 
                unpack=True, usecols=[0, 2, 3, 4])

        self.nsnap = nsnap 
        self.z = z_snap[(n_snaps.tolist()).index(nsnap)]        # redshift of snapshot
        self.t_cosmic = t_snap[(n_snaps.tolist()).index(nsnap)] # t_cosmic of snapshot 
        self.t_step = t_wid[(n_snaps.tolist()).index(nsnap)]    # Gyrs until next snapshot

    def AssignSFR(self, origin_nsnap, mass_bin=None): 
        ''' Assign SFR based on galaxy mass for snapshot 

        Assign galaxies as star-forming or quiescent
        based on quiescent fraction. 
        Then based on the assigned galaxy type 
        assign SSFR & SFR to the stellar masses
        '''
        # if AssignSFR is called 
        if self.mass == None:                   # snapshot not yet imported
            # import snapshot 
            ImportSnapshot(self, origin_nsnap) 
        else: 
            pass 

        # get mass bins
        if mass_bin == None:     # if mass bin isn't specified 
            mass_bins = util.simple_mass_bin()  # use simplest mass bins
        else: 
            raise NameError("not yet coded") 

        # remove galaxies below the minimum and maximum mass
        # and galaxies without children 
        mass_child_limit = (self.mass > min(mass_bins.mass_low)) & \
                (self.mass <= max(mass_bins.mass_high)) & \
                (self.child >= 0)  

        # only keep galaxies within the min and max mass
        self.sample_select(mass_child_limit, 
                columns = ['mass', 'parent', 'child', 'ilk', 'snap_index', 'pos'])          
        
        self.gal_type = np.array(['' for i in range(len(self.mass))]) 
        self.sfr = np.array([-999. for i in range(len(self.mass))]) 
        self.ssfr = np.array([-999. for i in range(len(self.mass))]) 
        # loop through mass bins and assign SFRs
        for i_m in range(mass_bins.nbins): 
            mass_bin_bool = (self.mass > mass_bins.mass_low[i_m]) & \
                    (self.mass <= mass_bins.mass_high[i_m])     # boolean list for mass range

            # indices of galaxies within mass range
            mass_bin_index = np.array(range(len(self.mass)))[mass_bin_bool] 
        
            mass_bin_ngal = len(self.mass[mass_bin_bool])       # Ngal in mass bin 
            if mass_bin_ngal == 0:                  # if there are no galaxies within mass bin
                continue 
    
            # get quiescent fraction for mass bin at z_snapshot 
            mass_bin_qf = util.get_fq(mass_bins.mass_mid[i_m], self.z, 
                    lit='cosmosinterp') 

            # Ngal,quiescent in mass bin
            mass_bin_n_q = int( np.rint(mass_bin_qf * np.float(mass_bin_ngal)) ) 

            # Ngal,active in mass bin 
            mass_bin_n_sf = mass_bin_ngal - mass_bin_n_q

            print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]
            print 'fQ = ', mass_bin_qf, ' Ngal = ', mass_bin_ngal, \
                    ' Ngal,Q = ', mass_bin_n_q, ' Ngal,SF = ', mass_bin_n_sf
            
            #### Quiescent ####
            # randomly sample mass_bin_n_q galaxies from the mass bin 
            mass_bin_q_index = random.sample(mass_bin_index, mass_bin_n_q) 

            self.gal_type[mass_bin_q_index] = 'quiescent'   # label galaxy type 
            '''
            sample SSFR from log-normal distribution with 0.2dex scatter centered at -12.0 
            '''
            self.ssfr[mass_bin_q_index] = 0.2 * np.random.randn(mass_bin_n_q) - 12.0
            # calculate SFR by SSFR+Mass
            self.sfr[mass_bin_q_index] = self.ssfr[mass_bin_q_index] + self.mass[mass_bin_q_index]

            #### Star-Forming ####
            mass_bin_sf_index = [x for x in mass_bin_index if x not in mass_bin_q_index]
            
            self.gal_type[mass_bin_sf_index] = 'star-forming'   # label galaxy type 
            '''
            get average and scatter of SF main sequence 
            '''
            [sf_avg_sfr, sf_sig_sfr] = util.get_sfr_mstar_z(mass_bins.mass_mid[i_m], self.z) 
            self.sfr[mass_bin_sf_index] = sf_sig_sfr * np.random.randn(mass_bin_n_q) - sf_avg_sfr 
            self.ssfr[mass_bin_sf_index] = \
                    self.sfr[mass_bin_sf_index] - self.mass[mass_bin_sf_index]
        
    def writeout(self): 
        ''' simply outputs specified columns to file
        '''
        if self.sfr == None:
            # output columns 
            columns = ['mass', 'parent', 'child', 'ilk', 'snap_index']
        else:       # if SFR/SSFR have been assigned
            # output columns 
            columns = ['mass', 'sfr', 'ssfr', 'gal_type', 
                    'parent', 'child', 'ilk', 'snap_index']

        col_fmt = []                # column format 
        for column in columns:      # combine data into list for output  
            column_attr = getattr(self, column)
            try: 
                column_data
            except NameError: 
                column_data = [column_attr] 
            else: 
                column_data.append( column_attr ) 
                
            # column format
            if column in ('mass', 'sfr', 'ssfr'):
                col_fmt.append('%10.5f')    # float
            elif column in ('parent', 'child', 'ilk', 'snap_index'): 
                col_fmt.append('%d')        # integer
            elif column in ('gal_type'): 
                col_fmt.append('%s')        # string
            
            print column, len(column_attr), type(column_attr[0])

        snap_info = ''.join(['n_snap = ', str(self.nsnap), ' z = ', str(self.z), 
            ' t_cosmic = ', str(self.t_cosmic), ' \n'])
        output_header = snap_info+', '.join(columns)    # list the columns
        
        if self.sfr != None: 
            output_file = ''.join(['/data1/hahn/central_quenching/', 
                'cenque_centrals_snapshot', str(self.nsnap), '_sfpropassign.dat']) 
        else: 
            output_file = ''.join(['/data1/hahn/central_quenching/', 
                'cenque_centrals_snapshot', str(self.nsnap), '.dat']) 
        
        print col_fmt
        np.savetxt(output_file, np.column_stack(column_data), 
            fmt=col_fmt, delimiter='\t', header=output_header)

    def sample_select(self, bool, columns=None, silent=False): 
        ''' Given boolean list remove False items from data
        '''

        if not silent: 
            n_remove = len(bool) - len(bool[bool == True])
            print 'Removing ', n_remove, ' elements from ', len(bool)  
    
        if columns == None:         # if data columns aren't specified
            data_columns = ['mass', 'sfr', 'ssfr', 'gal_type', 
                    'parent', 'child', 'ilk', 'snap_index', 'pos']
        else: 
            data_columns = columns 

        n_list = len(bool) 
        
        for column in data_columns:         # loop through columns and only keep true values
            attr_list = getattr(self, column) ;

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

if __name__=='__main__': 
    snap13 = CenQue() 
    snap13.ImportSnap(13) 
    print snap13.__dict__.keys() 
    print snap13.z
    print snap13.t_cosmic
    print snap13.sfr
    snap13.AssignSFR(13)
    print snap13.__dict__.keys() 
    print snap13.z
    print snap13.t_cosmic
    print snap13.gal_type
    print snap13.ssfr
    print snap13.sfr
    print snap13.mass
    snap13.writeout()
