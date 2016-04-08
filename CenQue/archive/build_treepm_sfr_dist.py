

class treepm_sfr: 
    def __init__(self): 
        pass
    def get_massbin(self, n_massbin, minmass='9.5', maxmass='11.5', delta_mass='0.01') :
        pass
    def get_nsnap_z(self, z): 
        '''
        Given redshift, get closest n_snap 
        '''
        nsnap_table_file = 'snapshot_table.dat'
        nsnap_table = np.loadtxt(nsnap_table_file)
        return  



