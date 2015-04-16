'''

Central_Quenchign project


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

        self.nsnap = None       # n_snapshot 
        self.zsnap = None           # z_snapshot
        self.t_cosmic = None    # t_cosmic for snapshot
        self.t_step = None      # t_cosmic step 

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
        self.zsnap = z_snap[(n_snaps.tolist()).index(nsnap)]        # redshift of snapshot
        self.t_cosmic = t_snap[(n_snaps.tolist()).index(nsnap)] # t_cosmic of snapshot 
        self.t_step = t_wid[(n_snaps.tolist()).index(nsnap)]    # Gyrs until next snapshot

    def AssignSFR(self, origin_nsnap, mass_bin=None, **kwargs): 
        ''' Assign SFR based on galaxy mass for snapshot 

        Assign galaxies as star-forming or quiescent
        based on quiescent fraction. 
        Then based on the assigned galaxy type 
        assign SSFR & SFR to the stellar masses
        '''
        # if AssignSFR is called 
        if self.mass == None:                   # snapshot not yet imported
            # import snapshot 
            self.ImportSnap(origin_nsnap) 
        else: 
            pass 

        # get mass bins
        if mass_bin == None:     # if mass bin isn't specified 
            mass_bins = util.simple_mass_bin()  # use simplest mass bins
        else: 
            raise NameError("not yet coded") 
        
        if (min(self.mass) < min(mass_bins.mass_low)) or (max(self.mass) > max(mass_bins.mass_high)): 
            # remove galaxies below the minimum and maximum mass
            # and galaxies without children 
            mass_child_limit = (self.mass > min(mass_bins.mass_low)) & \
                    (self.mass <= max(mass_bins.mass_high)) & \
                    (self.child >= 0)  

            # only keep galaxies within the min and max mass
            self.sample_select(mass_child_limit, 
                    columns = ['mass', 'parent', 'child', 'ilk', 'snap_index'])          

        if self.gal_type == None:         
            self.gal_type = np.array(['' for i in range(len(self.mass))], dtype='|S16') 
            self.sfr = np.array([-999. for i in range(len(self.mass))]) 
            self.ssfr = np.array([-999. for i in range(len(self.mass))]) 

        for i_m in range(mass_bins.nbins):              # loop through mass bins and assign SFRs
        
            mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
            mass_bin_high = round(mass_bins.mass_high[i_m], 2) 
            # boolean list for mass range
            mass_bin_bool = (self.mass > mass_bin_low) & \
                    (self.mass <= mass_bin_high) &\
                    (self.gal_type == '')           # only assign for galaxies without galtype

            # indices of galaxies within mass range
            mass_bin_index = np.array(range(len(self.mass)))[mass_bin_bool] 
        
            mass_bin_ngal = len(self.mass[mass_bin_bool])       # Ngal in mass bin 
            if mass_bin_ngal == 0:                  # if there are no galaxies within mass bin
                continue 
    
            # get quiescent fraction for mass bin at z_snapshot 
            mass_bin_qf = util.get_fq(mass_bins.mass_mid[i_m], self.zsnap, 
                    lit=kwargs['fq']) 

            # Ngal,quiescent in mass bin
            mass_bin_n_q = int( np.rint(mass_bin_qf * np.float(mass_bin_ngal)) ) 

            # Ngal,active in mass bin 
            mass_bin_n_sf = mass_bin_ngal - mass_bin_n_q

            print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]
            print 'fQ = ', mass_bin_qf, ' Ngal = ', mass_bin_ngal, \
                    ' Ngal,Q = ', mass_bin_n_q, ' Ngal,SF = ', mass_bin_n_sf
            
            #### Quiescent ####
            if mass_bin_n_q > 0: 
                # randomly sample mass_bin_n_q galaxies from the mass bin 
                try: 
                    mass_bin_q_index = random.sample(mass_bin_index, mass_bin_n_q) 
                except ValueError: 
                    mass_bin_q_index = mass_bin_index 

                self.gal_type[mass_bin_q_index] = 'quiescent'   # label galaxy type 
                '''
                sample SSFR from log-normal distribution with 0.2dex scatter centered at -12.0 
                '''
                self.ssfr[mass_bin_q_index] = 0.2 * np.random.randn(mass_bin_n_q) - 12.0
                # calculate SFR by SSFR+Mass
                self.sfr[mass_bin_q_index] = self.ssfr[mass_bin_q_index] + self.mass[mass_bin_q_index]

            #### Star-Forming ####
            if mass_bin_n_sf > 0: 
                try: 
                    mass_bin_sf_index = [x for x in mass_bin_index if x not in mass_bin_q_index]
                except NameError:       # if there aren't any quiescent galaxies
                    mass_bin_sf_index = mass_bin_index
                
                self.gal_type[mass_bin_sf_index] = 'star-forming'   # label galaxy type 
                '''
                get average and scatter of SF main sequence 
                '''
                [sf_avg_sfr, sf_sig_sfr] = util.get_sfr_mstar_z(mass_bins.mass_mid[i_m], self.zsnap) 

                print 'SF Average(SFR) = ', sf_avg_sfr, ' sigma_SFR = ', sf_sig_sfr

                self.sfr[mass_bin_sf_index] = sf_sig_sfr * np.random.randn(mass_bin_n_sf) + sf_avg_sfr 
                self.ssfr[mass_bin_sf_index] = \
                        self.sfr[mass_bin_sf_index] - self.mass[mass_bin_sf_index]

        # check for SFR/SSFR assign fails
        assign_fail = (self.sfr == -999.0) | (self.ssfr == -999.0)
        if len(self.sfr[assign_fail]) > 0: 
            raise NameError('asdfasdfasdfasdf')

    def readin(self, **kwargs): 
        ''' Read in cenque data written by CenQue

        Get snapshot and column info from header then import column 
        into CenQue data 

        overly convoluted way to read in data -.-
        '''
        # kwargs specifies the input file 
        input_file = util.cenque_file( **kwargs ) 

        with open(input_file, 'r') as f:
            for i, line in enumerate(f): 
                if i == 0: 
                    header_data = line.split(',') 
                    header_column = ['nsnap', 'zsnap', 't_cosmic', 't_step'] 

                    for i_head, head in enumerate(header_column): 
                        if head == 'nsnap': 
                            header_datum = int((header_data[i_head]).split('=')[1])
                        else: 
                            header_datum = float((header_data[i_head]).split('=')[1])

                        # save snapshot information from header 
                        setattr(self, head, header_datum)       
                elif i == 1: 
                    data_columns = map(str.strip, (line.strip()).split(','))
                else: 
                    break 
    
        flt_cols = [] 
        int_cols = []
        str_cols = [] 
        col_fmt = [] 
        for i_col, column in enumerate(data_columns): 
            if column in ('mass', 'sfr', 'ssfr', 'tau', 
                    'parent_sfr', 'parent_mass', 'q_ssfr'):                       # floats 
                flt_cols.append(i_col)
                col_fmt.append('float64')
            elif column in ('parent', 'child', 'ilk', 'snap_index'):    # ints
                int_cols.append(i_col)
                col_fmt.append('int64')
            elif column in ('gal_type'):                                # strings
                str_cols.append(i_col) 
                col_fmt.append('|S16')
            else: 
                raise NameError('Something went wrong here') 
    
        data_float = np.loadtxt(input_file, skiprows=2, delimiter=',', 
                unpack=True, usecols=flt_cols) 
        data_string = np.loadtxt(input_file, skiprows=2, delimiter=',', 
                unpack=True, usecols=str_cols, dtype='|S16') 
        data_int = np.loadtxt(input_file, skiprows=2, delimiter=',', 
                unpack=True, usecols=int_cols, dtype='int64') 

        # care taken to keep order
        i_flt = 0
        i_int = 0
        for i_col, column in enumerate(data_columns): 
            if column in ('mass', 'sfr', 'ssfr', 'tau',
                    'parent_sfr', 'parent_mass'):                       # floats 
                if len(flt_cols) == 1:  
                    setattr(self, column, data_float) 
                else:
                    setattr(self, column, data_float[i_flt]) 
                    i_flt += 1
            elif column in ('parent', 'child', 'ilk', 'snap_index'):    # ints
                if len(int_cols) == 1: 
                    setattr(self, column, data_int) 
                else: 
                    setattr(self, column, data_int[i_int]) 
                    i_int += 1
            elif column in ('gal_type'):                                # strings
                setattr(self, column, map(str.strip, data_string))       # strip string

    def writeout(self, **kwargs): 
        ''' simply outputs specified columns to file
        '''
        if self.sfr == None:
            # output columns 
            columns = ['mass', 'parent', 'child', 'ilk', 'snap_index']
        elif 'columns' in kwargs.keys(): 
            columns = kwargs['columns']
        else:       # if SFR/SSFR have been assigned
            # output columns 
            columns = ['mass', 'sfr', 'ssfr', 'gal_type', 
                    'parent', 'child', 'ilk', 'snap_index']
            
        n_cols = len(columns)       # number of columns 
        col_fmt = []                # column format 

        for column in columns:      # combine data into list for output  
            column_attr = getattr(self, column)
            try:                # save column data 
                column_data
            except NameError: 
                column_data = [column_attr] 
            else: 
                column_data.append( column_attr ) 
                
            if column in ('mass', 'sfr', 'ssfr', 'tau', 
                    'parent_sfr', 'parent_mass', 'q_ssfr'):        # column format
                col_fmt.append('%10.5f')    # float
            elif column in ('parent', 'child', 'ilk', 'snap_index'): 
                col_fmt.append('%d')        # integer
            elif column in ('gal_type'): 
                col_fmt.append('%s')        # string
            else: 
                raise NameError("column doesn't exist") 
        
        # header information 
        snap_info = ''.join(['n_snap = ', str(self.nsnap), ', z = ', str(self.zsnap), 
            ', t_cosmic = ', str(self.t_cosmic), ', t_step = ', str(self.t_step), ' \n'])
        output_header = snap_info+', '.join(columns)+'\n'    # list the columns
        
        output_file = util.cenque_file( **kwargs ) 
        
        # write out file
        f = open( output_file, 'w' )
        f.write( output_header )
        for i in range(len(column_data[0])): 
            f.write(', \t'.join( [ str(column_data[j][i]) for j in range(n_cols) ] )+'\n' ) 
        f.close() 

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

def EvolveCenQue(origin_nsnap, final_nsnap, mass_bin=None, **kwargs): 
    ''' Evolve SF properties from origin_nsnap to final_nsnap 
    '''

    if mass_bin == None: 
        mass_bins = util.simple_mass_bin()  # use simplest mass bins
    else: 
        raise NotImplementedError("not yet coded") 
  
    # import original snap SF prop 
    parent_cq = CenQue()
    parent_cq.readin(nsnap=origin_nsnap, file_type='sf assign', **kwargs)   

    print 'Quiescent Fraction = ', np.float(len(parent_cq.gal_type[parent_cq.gal_type == 'quiescent']))/np.float(len(parent_cq.gal_type)) 

    for i_step in range(0, origin_nsnap - final_nsnap):    # evolve snapshot by snapshot 

        i_snap = origin_nsnap - i_step          # current Nsnap 
        child_snap = i_snap-1

        child_cq = CenQue() 
        child_cq.readin(nsnap=child_snap)         # read in snapshot SHAM masses
        
        # remove galaxies beyond mass bins
        child_mass_limit = (child_cq.mass > min(mass_bins.mass_low)) & \
                (child_cq.mass <= max(mass_bins.mass_high))
        child_cq.sample_select( child_mass_limit, 
                columns = ['mass', 'parent', 'child', 'ilk', 'snap_index']  
                )                   
        n_child = len(child_cq.mass)     # number of children left 
        print 'Snapshot ', child_snap, ' has ', n_child, ' Galaxies'

        # set up columns for assignment 
        child_cq.gal_type = np.array(['' for i in range(n_child)], dtype='|S16') 
        child_cq.sfr = np.array([-999. for i in range(n_child)]) 
        child_cq.ssfr = np.array([-999. for i in range(n_child)]) 
        child_cq.q_ssfr = np.array([-999. for i in range(n_child)]) 
        child_cq.tau = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_sfr = np.array([-999. for i in range(n_child)]) 
        child_cq.parent_mass = np.array([-999. for i in range(n_child)]) 

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
        print len(child_cq.gal_type[child_cq.gal_type == '']), ' out of ',\
                len(child_cq.gal_type), ' child galaxies are orphans' 
        (child_cq.parent_sfr)[child_indx] = [(parent_cq.sfr)[i] for i in parent_indx]
        (child_cq.parent_mass)[child_indx] = [(parent_cq.mass)[i] for i in parent_indx]
        
        try: 
            parent_cq.tau
        except AttributeError:
            pass
        else: 
            # children inherit tau 
            (child_cq.tau)[child_indx] = [(parent_cq.tau)[i] for i in parent_indx]
            has_tau = (child_cq.tau)[child_indx] != -999.0               # no tau 
            still_quenching_indx = np.array(child_indx)[has_tau]
            
            # for galaxies with tau, keep them quenching!
            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[still_quenching_indx]
                ))                                              # SFR quenching amount  
            child_cq.sfr[still_quenching_indx] = \
                    child_cq.parent_sfr[still_quenching_indx] + tau_quench 
            child_cq.ssfr[still_quenching_indx] = \
                    child_cq.sfr[still_quenching_indx] - child_cq.mass[still_quenching_indx]
            print child_cq.gal_type[still_quenching_indx]
            
            # children inherit final quenched SSFR 
            (child_cq.q_ssfr)[child_indx] = [(parent_cq.q_ssfr)[i] for i in parent_indx]

        # Star-forming Children ----------------------------------------
        sf_child = (child_cq.gal_type[child_indx] == 'star-forming')
        sf_child_indx = np.array(child_indx)[sf_child]                  # index stuff
        sf_child_parent_indx = np.array(parent_indx)[sf_child]
        
        # SFR evolution amount
        child_sfr, child_sig_sfr = util.get_sfr_mstar_z(child_cq.mass[sf_child_indx], 
                child_cq.zsnap)
        parent_sfr, parent_sig_sfr = util.get_sfr_mstar_z(parent_cq.mass[sf_child_parent_indx], 
                parent_cq.zsnap) 
        dSFR = child_sfr - parent_sfr
        #util.get_sfr_mstar_z([(parent_cq.mass)[i] for i in sf_child_parent_indx], parent_cq.z) 
        
        # evolve sf children SFR
        child_cq.sfr[sf_child_indx] = child_cq.parent_sfr[sf_child_indx] + dSFR   
        child_cq.ssfr[sf_child_indx] = child_cq.sfr[sf_child_indx] - child_cq.mass[sf_child_indx]

        # Quiescent Children --------------------------------------------
        q_child = (child_cq.gal_type[child_indx] == 'quiescent') 
        q_child_indx = np.array(child_indx)[q_child]
        q_child_parent_indx = np.array(parent_indx)[q_child]

        # keep SSFR same 
        child_cq.ssfr[q_child_indx] = parent_cq.ssfr[q_child_parent_indx]
        child_cq.sfr[q_child_indx] = child_cq.ssfr[q_child_indx] + \
                child_cq.mass[q_child_indx]
       
        # loop through mass bin and evolved active children with parents to match 
        # mass bin fq
        for i_m in range(mass_bins.nbins):              
            
            # boolean list for mass range
            mass_bin_bool = (child_cq.mass > mass_bins.mass_low[i_m]) & \
                    (child_cq.mass <= mass_bins.mass_high[i_m]) & \
                    (child_cq.gal_type != '') 
            print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]

            # indices of galaxies within mass range
            mass_bin_index = child_cq.get_index(mass_bin_bool) 
        
            mbin_ngal = len(child_cq.mass[mass_bin_bool])       # Ngal in mass bin 
            if mbin_ngal == 0:                  # if there are no galaxies within mass bin
                continue 
            
            if 'fq' in kwargs.keys():           
                mbin_qf = util.get_fq(mass_bins.mass_mid[i_m], child_cq.zsnap, 
                    lit=kwargs['fq']) 
            else:                       # make sure qf is specified
                raise TypeError("Specify quiescent fraction type") 

            print 'z = ', child_cq.zsnap, ' M* = ', mass_bins.mass_mid[i_m], ' fq = ', mbin_qf

            mbin_exp_n_q = int( np.rint(mbin_qf * np.float(mbin_ngal)) )    # Ngal,Q_expected
            
            # number of SF galaxies that need to be quenched 
            child_gal_type = child_cq.gal_type[mass_bin_bool] 

            ngal_2quench = mbin_exp_n_q - len(child_gal_type[child_gal_type == 'quiescent'])  
            print ngal_2quench, ' SF galaxies will be quenched'
            if ngal_2quench <= 0: 
                print ngal_2quench, ' is negative'
                continue 
            
            # randomly sample mass_bin_n_q galaxies from the mass bin 
            mbin_sf_index = child_cq.get_index( 
                    [ mass_bin_bool & (child_cq.gal_type == 'star-forming')] 
                    ) 

            if len(child_gal_type[child_gal_type == 'star-forming']) < ngal_2quench: 
                quench_index = mbin_sf_index
            else: 
                quench_index = random.sample(mbin_sf_index, ngal_2quench) 
            
            child_cq.gal_type[quench_index] = 'quiescent'  # boom quenched 

            if 'tau' in kwargs.keys(): 
                child_cq.tau[quench_index] = util.get_quenching_efold(
                        child_cq.parent_mass[quench_index], type=kwargs['tau'])
            else: 
                raise TypeError('specify quenching e-fold: tau = instant, constant, linear') 

            tau_quench = np.log10( np.exp( 
                -(child_cq.t_cosmic - parent_cq.t_cosmic) / child_cq.tau[quench_index]
                ))      # SFR quenching amount  

            child_cq.sfr[quench_index] = child_cq.parent_sfr[quench_index] + tau_quench 
            child_cq.ssfr[quench_index] = child_cq.sfr[quench_index] - child_cq.mass[quench_index]
            child_cq.q_ssfr[quench_index] = 0.2 * np.random.randn(ngal_2quench) - 12.0
        
        # deal with orphans
        print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
        child_cq.AssignSFR(child_cq.nsnap, **kwargs) 
        print len(child_cq.gal_type[child_cq.gal_type == '']), ' child galaxies are orphans' 
        print child_cq.mass[child_cq.gal_type == '']
        
        # deal with galaxies that are being quenched beyond the designated final quenching  
        over_quenched = child_cq.ssfr < child_cq.q_ssfr 
        child_cq.ssfr[over_quenched] = child_cq.q_ssfr[over_quenched]
        child_cq.tau[over_quenched] = -999.0            # done quenching 

        child_cq.writeout(nsnap=child_cq.nsnap, file_type='evol from '+str(origin_nsnap), 
                columns = ['mass', 'sfr', 'ssfr', 'gal_type', 'tau', 'q_ssfr', 
                    'parent_sfr', 'parent_mass', 'parent', 'child', 'ilk', 'snap_index'], 
                **kwargs)  

        parent_cq = child_cq
        print 'Quiescent Fraction = ', np.float(len(parent_cq.gal_type[parent_cq.gal_type == 'quiescent']))/np.float(len(parent_cq.gal_type)) 

def build_cenque_importsnap(**kwargs): 
    ''' 
    '''
    for i_snap in [13]: 
        snap = CenQue() 
        snap.readin(nsnap=i_snap)
        snap.AssignSFR(i_snap, **kwargs) 
        snap.writeout(nsnap=i_snap, file_type='sf assign', **kwargs)

if __name__=='__main__': 
    #build_cenque_importsnap(fq='cosmosinterp') 
    #EvolveCenQue(13, 1, fq='cosmosinterp', tau='instant') 
    #EvolveCenQue(13, 1, fq='cosmosinterp', tau='constant') 
    #EvolveCenQue(13, 1, fq='cosmosinterp', tau='linear') 

    #build_cenque_importsnap(fq='wetzel') 
    EvolveCenQue(13, 1, fq='wetzel', tau='instant') 
    EvolveCenQue(13, 1, fq='wetzel', tau='constant') 
    EvolveCenQue(13, 1, fq='wetzel', tau='linear') 
