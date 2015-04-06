'''

Documentation here 

Author(s): ChangHoon Hahn 

'''

import numpy as np
import conque_utility as util

class CenQue: 
    '''
    Central quenching (cenque) galaxy catalogue
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

    def ImportSnapshot(self, nsnap): 
        ''' Import mass values from TreePM --> SHAM of snapshot 
        '''
        snapshot_dir = '/data1/hahn/wetzel_tree/'
        snapshot_file = ''.join([ snapshot_dir, 
            'subhalo_sham_centrals_snapshot', str(nsnap), '.fits' 
            ]) 
        snapshot = mrdfits(snapshot_file)   # read snapshot file
        
        # import snapshot values
        self.mass = snapshot.mass 
        self.parent = snapshot.parent
        self.child = snapshot.child
        self.parent = snapshot.parent
        self.snap_index = snapshot.index
        self.pos = snap.index 

    def AssignSFR(self, mass_bin=None): 
        ''' Assign SFR based on galaxy mass 
        '''

        # get mass bins
        if mass_bin = None:     # if mass bin isn't specified 
            mass_bins = util.simple_mass_bin()  # use simplest mass bins
        else: 
            raise NameError("not yet coded") 
        
        # loop through mass bins and assign SFRs
        for i_m in range(mass_bins.nbins): 
            mass_range_bool = (self.mass > mass_bins.mass_low[i_m]) & \
                    (self.mass <= mass_bins.mass_high[i_m])     # boolean list for mass range

            mass_bin_qf = 
