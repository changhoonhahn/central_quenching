
class ssfr: 
    def __init__(self, nsnap, catalog='mock', nsnap_i=13, mbin=0.01, tau='constant'): 
        '''
        Read SSFR in distribution 
        '''
        # for evoled treepm mocks
        if catalog == 'mock': 
            # quenching timescale 
            if tau == 'constant':       # constant tau 
                tau_str = '_constanttau'
            elif tau == 'linear':       # tau linear with mass 
                tau_str = '_lineartau'
            elif tau == 'instant': 
                tau_str = '_instanttau'
            
            # get file name
            file_dir = '/data1/hahn/wetzel_tree/'
            if nsnap == nsnap_i: 
                file_name = file_dir+'subhalo_sham_centrals_snapshot'+str(nsnap_i)+'_ssfr_mbin'+str(mbin)+'.fits'
            else: 
                file_name = file_dir+'subhalo_sham_centrals_snapshot'+str(nsnap)+'_ssfr_mbin'+str(mbin)+\
                        '_evolfrom_snapshot'+str(nsnap_i)+tau_str+'.fits'
            self.file_name = file_name          # store file name 

            # read file and extract relevant data ---------------------------------------
            ssfr_data = mrdfits(file_name) 

            # columns we want to extract
            data_columns = ['mass', 'ssfr']
            self.columns = data_columns 

            for column in data_columns: 
                setattr(self, column, getattr(ssfr_data, column))

        # for Jeremy's SDSS group catalog 
        elif catalog == 'data': 
            # nsnap determines absmag and mass cuts, which are hardcoded
            if nsnap == 0: 
                absmag_cut = 18
                mass_cut = 9.4
            elif nsnap == 1: 
                absmag_cut = 19
                mass_cut = 9.8
            elif nsnap == 2: 
                absmag_cut = 20
                mass_cut = 10.2
            else:
                raise NameError("nsnap = 0, 1, or 2")

            file_dir = '/data1/hahn/group_catalog/'
            file_name = file_dir+'massSFR_clf_groups_M'+str(absmag_cut)+'_'+str(mass_cut)+'_D360.galdata_corr.central.fits'
            
            self.file_name = file_name
            self.absmag_cut = absmag_cut
            self.mass_cut = mass_cut
            
            # read file and extract relevant data ---------------------------------------
            ssfr_data = mrdfits(file_name)

            # columns we want to extract
            data_columns = ['mass', 'ssfr']
            self.columns = data_columns 

            for column in data_columns: 
                setattr(self, column, getattr(ssfr_data, column))

def get_z_nsnap(nsnap): 
    '''
    redshift of snapshot 
    '''
    iz, redshift, tcosmis = np.loadtxt('/home/users/hahn/research/pro/tinker/central_quenching/snapshot_table.dat', unpack=True, 
            usecols=[0,2,3]) 
    return  redshift[iz == nsnap]

def get_fq(Mstar, redshift, scheme='cosmosinterp'): 
    '''
    calculate the quiescent/red fraction at Mstar and z
    '''
    if "cosmos" in scheme.lower(): 
        cosmos_zbins = [0.36, 0.66, 0.88]
        cosmos_dir = '/data1/hahn/wetzel_tree/'

    if scheme.lower() == 'cosmosinterp': 
        cosmos_zbins = [0.36, 0.66, 0.88]
       
        fqMstar = np.zeros(len(cosmos_zbins))
        for i_z, cosmos_zbin in enuemrate(cosmos_zbins): 
            cosmos_file = cosmos_dir+'qf_z'+str(cosmos_zbin)+'cen.dat'
            comsos_mass, cosmos_qf = np.loadtxt(cosmos_file, unpack=True, usecols=[0,1])

            fqMstar[i_z] = np.interp(Mstar, cosmos_mass, cosmos_qf)
        return np.interp(redshift, cosmos_zbin, fqMstar) 

def plot_evomock_ssfr_dist(nsnap=[1], tau='constant'): 
    ''' plot SSFR distribution for the TREEPM evolved mocks
    '''
    # configure figure 
    prettyplot() 
    fig = plt.figure(1, figsize=(15,15))
    fig.subplots_adjust(hspace=0., wspace=0.)
    sub0 = fig.add_subplot(221) 
    sub1 = fig.add_subplot(222) 
    sub2 = fig.add_subplot(223) 
    sub3 = fig.add_subplot(224) 

    # set pretty colors: 
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
            (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
            (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
            (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
            (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  

    for i in range(len(tableau20)):  
        r, g, b = tableau20[i]  
        tableau20[i] = (r / 255., g / 255., b / 255.)  

    # mass bin 
    mass_bins = [[9.5, 10.0], [10.0, 10.5], [10.5, 11.0], [11.0, 11.5]]
    
    for ii_snap, i_snap in enumerate(nsnap): 
        # read in snapshot ssfr data 
        snap_data = ssfr(i_snap, catalog='mock', tau=tau) 
        isnap_z = get_z_nsnap(i_snap) 

        # for each mass bin 
        for i_mass, mass_bin in enumerate(mass_bins): 
            mass_limit = (snap_data.mass >= mass_bin[0]) & (snap_data.mass < mass_bin[1])

            # SSFR bins
            ssfr_hist_bin = [round(-13.0+i*0.2,10) for i in range(0,31)]
            ssfr_hist_bin_mid = np.array([round(0.5*(ssfr_hist_bin[i]+ssfr_hist_bin[i+1]),10) for i in range(0,30)]) 

            # SSFR histogram 
            ssfr_hist, ssfr_bin_edges = np.histogram(snap_data.ssfr[mass_limit], bins=ssfr_hist_bin, density=True)

            ssfr_hist_label = r'$z='+str(isnap_z[0])+'$'
            plt.subplot(2,2,i_mass+1)
            plt.plot(ssfr_hist_bin_mid, ssfr_hist, color=tableau20[ii_snap], lw=2, label=ssfr_hist_label) 
            plt.text(-11.0, 1.4, r'$log \; M_{star} = ['+str(mass_bin[0])+', '+str(mass_bin[1])+']$', fontsize=24)
    
    # overplot SDSS SSFR distribution  
    linestyles = ['-', '--', '-.'] 
    for i_sdss, style in enumerate(linestyles): 
        sdss_data = ssfr(i_sdss, catalog='data')
        for i_mass, mass_bin in enumerate(mass_bins): 
            mass_limit = (sdss_data.mass >= mass_bin[0]) & (sdss_data.mass < mass_bin[1])

            # SSFR bins
            ssfr_hist_bin = [round(-13.0+i*0.2,10) for i in range(0,31)]
            ssfr_hist_bin_mid = np.array([round(0.5*(ssfr_hist_bin[i]+ssfr_hist_bin[i+1]),10) for i in range(0,30)]) 

            # SSFR histogram 
            ssfr_hist, ssfr_bin_edges = np.histogram(sdss_data.ssfr[mass_limit], bins=ssfr_hist_bin, density=True)

            ssfr_hist_label = 'SDSS Group Catalog ('+str(sdss_data.absmag_cut)+','+str(sdss_data.mass_cut)+')'
            plt.subplot(2,2,i_mass+1)
            plt.plot(ssfr_hist_bin_mid, ssfr_hist, color='black', ls=style, lw=3, label=ssfr_hist_label) 

    # set axes limits
    sub0.set_xlim([-13.0, -7.0])
    sub1.set_xlim([-13.0, -7.0])
    sub2.set_xlim([-13.0, -7.0])
    sub3.set_xlim([-13.0, -7.0])

    sub0.set_ylim([0.0, 1.6])
    sub1.set_ylim([0.0, 1.6])
    sub2.set_ylim([0.0, 1.6])
    sub3.set_ylim([0.0, 1.6])
    # set y-axes labels
    sub0.set_ylabel(r'$P(log \; SSFR)$') 
    sub1.set_yticklabels([])
    sub2.set_ylabel(r'$P(log \; SSFR)$') 
    sub3.set_yticklabels([])

    # set y-axes labels
    sub0.set_xticklabels([])
    sub1.set_xticklabels([])
    sub2.set_xlabel(r'$log \; SSFR \;[yr^{-1}]$') 
    sub3.set_xlabel(r'$log \; SSFR \;[yr^{-1}]$') 
    # legends 
    sub0.legend(loc='upper left', prop={'size':'12'}) 
   
    fig_dir = '/home/users/hahn/research/figures/tinker/'
    fig_name = fig_dir+'subhalo_sham_centrals_ssfr_dist_allsnapshots_evolfrom_snapshot13_'+tau+'tau.png'
    print fig_name
    fig.savefig(fig_name, bbox_inches="tight") 
    fig.clear() 

def plot_fq_evol(scheme='cosmosinterp'): 
   pass 

if __name__ == "__main__": 
    plot_evomock_ssfr_dist(nsnap=range(1,14), tau='constant')
    plot_evomock_ssfr_dist(nsnap=range(1,14), tau='instant')
    plot_evomock_ssfr_dist(nsnap=range(1,14), tau='linear')
