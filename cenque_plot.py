'''

Plotting code for CentralQuenching (CenQue) project

Author(s): ChangHoon Hahn 

'''

import numpy as np
import matplotlib.pyplot as plt
import mpfit
import scipy.stats as scistat

#----- Local -----
from utility.plotting import prettyplot
from utility.plotting import prettycolors 
import cenque_utility as util 
import cenque as cq 
import cenque_groupcat as cq_group
import sf_mainseq as sfms 

# ssfr distribution ------------------------------------------
def plot_cenque_ssfr_dist(cenque, fig=None, **kwargs): 
    ''' Plot CenQue data
    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
    
    if fig == None: 
        fig = plt.figure(1, figsize=(25,8))
        fig.subplots_adjust(hspace=0., wspace=0.)
        new_plot = True
    else: 
        new_plot = False

    # mass bins of the panel        ( hardcoded ) 
    panel_mass_bins = [ 
            [9.7, 10.2], [10.2, 10.7], [10.7, 11.2], [11.2, 11.7]
            ]
    
    for i_mass, panel_mass in enumerate(panel_mass_bins):       # loop through each panel 

        mass_limit = (cenque.mass >= panel_mass[0]) & (cenque.mass < panel_mass[1])

        # SSFR histogram 
        ssfr_hist, ssfr_bin_edges = np.histogram(cenque.ssfr[mass_limit], 
                range=[-13.0, -7], bins=40, normed=True)
        
        # SSFR bins 
        ssfr_bin_low = ssfr_bin_edges[:-1]
        ssfr_bin_high = ssfr_bin_edges[1:]
        ssfr_bin_mid = [ 0.5*(ssfr_bin_low[i] + ssfr_bin_high[i]) 
                for i in range(len(ssfr_bin_low)) ] 
    
        if 'label' in kwargs: 
            ssfr_hist_label = kwargs['label']
        else: 
            try: 
                cenque.zsnap
            except AttributeError: 
                ssfr_hist_label = 'Centrals'     # redshift lable 
            else: 
                ssfr_hist_label = r'$z='+str(cenque.zsnap)+'$'     # redshift lable 
    
        if 'line_color' in kwargs: 
            line_color = kwargs['line_color']
        else: 
            try: 
                line_color = pretty_colors[cenque.nsnap]
            except TypeError: 
                line_color = 'black'

        if 'line_style' in kwargs:
            line_style = kwargs['line_style'] 
        else: 
            line_style = '-'

        sub = fig.add_subplot(1, 4, i_mass+1)       # panel subplot 
        sub.plot(ssfr_bin_mid, ssfr_hist, 
                color=line_color, lw=4, ls=line_style, label=ssfr_hist_label) 
        
        if new_plot == True:        # put mass labels in panels 
            plt.text(-10.75, 1.4, 
                    r'$\mathtt{log \; M_{*} = ['+str(panel_mass[0])+', '+str(panel_mass[1])+']}$', 
                    fontsize=24)

        # set axes limits
        sub.set_xlim([-13.0, -7.0])
        sub.set_ylim([0.0, 1.6])
        # set y-axes labels
        if i_mass == 0: 
            sub.set_ylabel(r'$P(log \; SSFR)$') 
        elif i_mass in [1, 2]:
            sub.set_yticklabels([])
            sub.set_xlabel(r'$log \; SSFR \;[yr^{-1}]$') 
        else: 
            sub.set_yticklabels([])
                
            try: 
                fig_leg.remove() 
            except UnboundLocalError: 
                pass
            
            fig_leg = sub.legend(loc='lower right', prop={'size':'12'})        # legends 

    return fig   

def plot_cenque_ssfr_dist_evolution(Mrcut=18, **kwargs): 
    ''' Plot evolution of the CenQue SSFR distribution 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 
    nsnaps : List of snapshot #s to plot  


    '''
    snap = cq.CenQue() 
    snap.readin(nsnap=13, file_type='sf assign', **kwargs)  # starting CenQue 
    ssfr_fig = plot_cenque_ssfr_dist(snap)      # plot!
    
    # determine which snapshots to plot 
    if 'nsnaps' in kwargs.keys():   
        nsnaps = kwargs['nsnaps']
    else: 
        nsnaps = [1] 

    for i_nsnap in nsnaps:  
        next_snap = cq.CenQue() 
        next_snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
        
        ssfr_fig = plot_cenque_ssfr_dist(next_snap, fig=ssfr_fig) 
    
    # overplot SDSS group catalog sSFR dist
    central_ssfr = cq_group.central_catalog(Mrcut=Mrcut, clobber=True) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, label= 'Mrcut = '+str(Mrcut)) 
    
    if 'sfms_slope' in kwargs.keys(): 
        slope_str = str("%.2f" % kwargs['sfms_slope']) 
        yint_str = str("%.2f" % kwargs['sfms_yint']) 
        sfms_str = '_sfms_slope'+slope_str+'_yint'+yint_str
    else: 
        sfms_str = ''

    if kwargs['tau'] == 'discrete': 
        fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
            'cenque_ssfr_evol_', '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] ), 
            'tau_', kwargs['fq'], 'fq_Mrcut', str(Mrcut),'.png']) 

    elif kwargs['tau'] == 'linefit':
        fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
            'cenque_ssfr_evol', sfms_str, 
            '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] ) , 
            'tau_', kwargs['fq'], 'fq_Mrcut', str(Mrcut),'.png']) 
    else: 
        fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
            'cenque_ssfr_evol_', '_'.join([str(t) for t in kwargs['tau']]), 'tau_', 
            kwargs['fq'], 'fq_Mrcut', str(Mrcut), sfms_str, '.png']) 
    ssfr_fig.savefig(fig_file, bbox_inches='tight') 
    ssfr_fig.clear() 

def plot_cenque_ssfr_dist_evolution_match2isedfit(Mrcut=18, **kwargs): 
    snap = cq.CenQue() 
    snap.readin(nsnap=13, file_type='sf assign', **kwargs) 
    ssfr_fig = plot_cenque_ssfr_dist(snap)

    for i_nsnap in reversed(range(1,13)):
        next_snap = cq.CenQue() 
        next_snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
        
        ssfr_fig = plot_cenque_ssfr_dist(next_snap, fig=ssfr_fig) 

    central_ssfr = cq_group.central_catalog_match2isedfit(Mrcut=Mrcut, clobber=True) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, label= 'Mrcut = '+str(Mrcut)) 
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'cenque_ssfr_evol_', kwargs['tau'], 'tau_', kwargs['fq'], 'fq_Mrcut', str(Mrcut),
        '_match2isedfit.png']) 
    ssfr_fig.savefig(fig_file, bbox_inches='tight') 
    ssfr_fig.clear() 

def plot_sdss_group_cat(): 
    ''' plot ssfr distribution for SDSS group catalogs
    '''
    pretty_colors = prettycolors() 
    central_ssfr = cq_group.central_catalog(Mrcut=18) 
    print np.min(central_ssfr.mass), np.max(central_ssfr.mass) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, 
            label='Mrcut=18', line_color=pretty_colors[1], line_style='-') 

    central_ssfr = cq_group.central_catalog(Mrcut=19) 
    print np.min(central_ssfr.mass), np.max(central_ssfr.mass) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, 
            label='Mrcut=19', line_color=pretty_colors[3], line_style='--') 

    central_ssfr = cq_group.central_catalog(Mrcut=20) 
    print np.min(central_ssfr.mass), np.max(central_ssfr.mass) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, 
            label='Mrcut=20', line_color=pretty_colors[5], line_style='-.') 

    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'sdss_groupcat_ssfr.png']) 
    ssfr_fig.savefig(fig_file, bbox_inches='tight')

def plot_sdss_group_cat_bestfit(Mrcut=19): 
    ''' plot ssfr distribution for SDSS group catalogs
    '''
    central_ssfr = cq_group.central_catalog(Mrcut=Mrcut, clobber=True) 

    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, label='Mrcut='+str(Mrcut), 
            line_color='black', line_style='--') 

    ssfr_hist, ssfr_bin_edges = np.histogram(central_ssfr.ssfr, 
            range=[-13.0, -7], bins=20, normed=True)

    ssfr_bin_low = ssfr_bin_edges[:-1]
    ssfr_bin_high = ssfr_bin_edges[1:]
    ssfr_bin_mid = [ 0.5*(ssfr_bin_low[i] + ssfr_bin_high[i]) 
            for i in range(len(ssfr_bin_low)) ] 
    
    output = cq_group.double_gaussian_fit(central_ssfr) 
    
    panel_mass_bins = [
            [10.0, 10.5], [10.5, 11.0], [11.0, 11.5]
            ]

    for i_mass, panel_mass in enumerate(panel_mass_bins):       # loop through each panel 

        sub = ssfr_fig.add_subplot(1, 3, i_mass+1)       # panel subplot 

        #sub.plot(ssfr_bin_mid, cq_group.double_gaussian(ssfr_bin_mid, (output[i_mass])[2]), 
        #        color='blue', lw=4, ls='-', label='Best Fit') 
        
        if i_mass == 0: 
            byeye_param = np.array([0.2, -11.8, -10.318, 0.19, 0.25]) 
        elif i_mass == 1: 
            byeye_param = np.array([0.3, -12.15, -10.57, 0.19, 0.25]) 
        elif i_mass == 2: 
            byeye_param = np.array([0.6, -12.38, -10.75, 0.175, 0.25]) 


        sub.plot(ssfr_bin_mid, cq_group.double_gaussian(ssfr_bin_mid, byeye_param), 
                color='blue', lw=4, ls='-', label='Best Fit') 
        sub.vlines(((output[i_mass])[2])[1], 0.0, 1.5, lw=4)
        sub.vlines(((output[i_mass])[2])[2], 0.0, 1.5, lw=4)

        sub.text(-11, 0.7, 'Best fit fQ = '+str(((output[i_mass])[2])[0]))
        sub.text(-11, 0.65, 'Best fit SSFR Q = '+str(((output[i_mass])[2])[1]))
        sub.text(-11, 0.6, 'Best fit SSFR SF = '+str(((output[i_mass])[2])[2]))
        
        sub.text(-11, 0.55, 'Best fit sig Q = '+str(((output[i_mass])[2])[3]))
        sub.text(-11, 0.5, 'Best fit sig SF = '+str(((output[i_mass])[2])[4]))

    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'sdss_groupcat_ssfr_bestift_mrcut', str(Mrcut), '.png']) 
    ssfr_fig.savefig(fig_file, bbox_inches='tight')
    ssfr_fig.clear()

def plot_cenque_sf_mainseq(): 
    ''' Given CenQue data, plot the Star-Forming Main Sequence
    '''
    masses = [10.0 + np.float(i)*0.1 for i in range(15)]        # masses
    redshift = [ 0.05, 0.2, 0.5, 0.9 ] 

    fig = plt.figure(1, figsize=(20,5)) 

    for i_z, z_in in enumerate(redshift): 
        sub = fig.add_subplot(1, len(redshift), i_z+1) 

        sfr_m_z = [] 
        for mass in masses: 
            sfr, sig_sfr = util.get_sfr_mstar_z(mass, z_in, lit='primusfit')
            sfr_m_z.append(sfr) 

        print masses 
        print sfr_m_z
        sub.scatter( masses, sfr_m_z, c='b' ) 
        if i_z == 0.0:
            sub.scatter([10.25, 10.75, 11.25], [-0.0777, 0.248, 0.483], c='k', s=6) 
        sub.text(10.0, 1.0, 'z = '+str(z_in))
        sub.set_xlim([9.5, 12.0]) 
        sub.set_ylim([-2.0, 2.0]) 
        sub.set_xlabel('log M') 
        if i_z == 0: 
            sub.set_ylabel('log SFR') 
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'sf_ms_evol.png']) 
    fig.savefig(fig_file) 

# quiescent fraction -----------------------------------
def plot_fq_evol(): 
    ''' plot fq evolution 
    '''

    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[zbin < 1.0]
    #zbin = [0.1*np.float(i) for i in range(1,10)]   # zbin 

    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    fig, subs = plt.subplots(1,4, figsize=[20, 5]) 
    subs = subs.ravel() 

    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    fq_types = ['cosmosinterp', 'cosmosfit', 'wetzel'] 
    for i_fq, fq_type in enumerate(fq_types): 
        for i_z, z in enumerate(zbin): 

            # plot fq(Mass) 
            fq_mass = [util.get_fq(mass_bin.mass_mid[i], z, lit=fq_type) 
                    for i in range(len(mass_bin.mass_mid))]
            subs[i_fq].plot(mass_bin.mass_mid, fq_mass, 
                    color=pretty_colors[i_z], lw=4, label='z = '+str(z) ) 

        subs[i_fq].set_title(fq_type) 

        subs[i_fq].set_xlim([9.0, 12.0])
        subs[i_fq].set_ylim([0.0, 1.0])

        subs[i_fq].set_xlabel('Mass') 

    for i_z, z in enumerate(zbin): 

        # plot fq(Mass) 
        fq_mass = [util.get_fq_alpha(mass_bin.mass_mid[i], z, -1.5) 
                for i in range(len(mass_bin.mass_mid))]
        subs[3].plot(mass_bin.mass_mid, fq_mass, 
                color=pretty_colors[i_z], lw=4, label='z = '+str(z) ) 

    subs[3].set_title('Evolved from z = 0.88') 

    subs[3].set_xlim([9.0, 12.0])
    subs[3].set_ylim([0.0, 1.0])

    subs[3].set_xlabel('Mass') 

    subs[0].set_ylabel('Quiescent Fraction') 
    subs[0].legend(loc='upper left') 
    return fig

def plot_fq_evol_w_geha(): 
    ''' plot fq evolution with literature results overplotted
    '''

    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[ zbin < 1.0 ] 
    #zbin = [0.1*np.float(i) for i in range(1,10)]   # zbin 
    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    fig, subs = plt.subplots(1,3, figsize=[15, 5]) 
    subs = subs.ravel() 

    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    ## load literature data 
    # modified tinker
    mod_tink_file = ''.join(['/data1/hahn/central_quenching/literature/modified_tinker_fq.dat']) 
    mod_tink_mass, mod_tink_fq = np.loadtxt(mod_tink_file, unpack=True, usecols=[0,1])   

    fq_types = ['cosmosinterp', 'cosmosfit', 'wetzel'] 
    for i_fq, fq_type in enumerate(fq_types): 
        for i_z, z in enumerate(zbin): 

            # plot fq(Mass) 
            fq_mass = [util.get_fq(mass_bin.mass_mid[i], z, lit=fq_type) 
                    for i in range(len(mass_bin.mass_mid))]
            subs[i_fq].plot(mass_bin.mass_mid, fq_mass, 
                    color=pretty_colors[i_z], lw=4, label='z = '+str(z) ) 

        subs[i_fq].plot(mod_tink_mass, mod_tink_fq, 
                color='black', lw=6, label='Modified Tinker Group' ) 

        subs[i_fq].set_title(fq_type) 

        subs[i_fq].set_xlim([9.0, 12.0])
        subs[i_fq].set_ylim([0.0, 1.0])

        subs[i_fq].set_xlabel('Mass') 
        subs[i_fq].set_xlabel('Mass') 

    subs[0].set_ylabel('Quiescent Fraction') 
    subs[0].legend(loc='upper left') 
    return fig

def plot_quenching_efold(): 
    ''' plot fq evolution 
    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    tau_types = ['instant', 'constant', 'linear'] 
    fig, subs = plt.subplots(1, len(tau_types), figsize=[15, 5]) 
    subs = subs.ravel() 
    for i_tau, tau in enumerate(tau_types): 
        # plot fq(Mass) 
        tau_mass = util.get_quenching_efold(np.array(mass_bin.mass_mid), type=tau) 
        subs[i_tau].plot(mass_bin.mass_mid, tau_mass, 
                color=pretty_colors[5], lw=4) 

        subs[i_tau].set_title(tau) 

        subs[i_tau].set_xlim([9.0, 12.0])
        subs[i_tau].set_ylim([0.0, 1.0])

        subs[i_tau].set_xlabel('Mass') 
        subs[i_tau].set_xlabel('Mass') 

    subs[0].set_ylabel('Quenching e-fold') 
    return fig

def plot_group_cat_bigauss_bestfit(): 
    plot_sdss_group_cat_bestfit(Mrcut=18)
    plot_sdss_group_cat_bestfit(Mrcut=19)
    plot_sdss_group_cat_bestfit(Mrcut=20)

# SF-MS ---------------------------
def plot_sfms_data(lowz_slope, lowz_yint, Mrcut=18): 
    ''' Plot StarForming Main Sequence from iSEDfit data and flexible SFMS fit function 

    Parameters
    ----------
    lowz_slope : designated low-z SDSS SF-MS slope
    lowz_yint : designated low-z SDSS SF-MS y-intercept

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 
    zbins = [ (0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8)]  # zbin

    fig, subs = plt.subplots(1, len(zbins), figsize=[20, 5]) 
    subs = subs.ravel() 

    # SDSS group catalog best fit 
    #groupcat_slope, groupcat_yint = sfms.sdss_groupcat_sfms_bestfit()
    #subs[0].plot(
    #        np.array(mass_bin.mass_mid), 
    #        (groupcat_slope * (np.array(mass_bin.mass_mid)-10.5)) + groupcat_yint, 
    #        c='k', lw=6, ls='--') 

    zmids, slopes, yints = sfms.get_sfmsfit_sfr(lowz_slope, lowz_yint, clobber=True) 

    fits = [] 
    for i_z, zbin in enumerate(zbins): 
        
        mass, avg_sfrs, var_sfrs = [], [], [] 
        for i_mass in range(len(mass_bin.mass_low)): 

            avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z_envcount(mass_bin.mass_mid[i_mass], 
                    0.5*(zbin[0] + zbin[1]))

            if ngal < 10: 
                continue    # skip mass bin if there aren't many galaxies
            if mass_bin.mass_mid[i_mass] < 9.5: 
                continue    # skip low mass

            mass.append(mass_bin.mass_mid[i_mass]) 
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)
        
            
        if i_z == 0:    # for SDSS panel 
            #centrals = cq_group.central_catalog(Mrcut=Mrcut) 
            #subs[i_z].hist2d(centrals.mass, centrals.sfr, bins=[30, 1000])
            centrals = sfms.sf_centrals(Mrcut=Mrcut) 
            subs[i_z].hist2d(centrals.mass, centrals.sfr, bins=[30, 100])

        subs[i_z].errorbar(mass, avg_sfrs, yerr=var_sfrs, 
                lw=4, c=pretty_colors[1], label='Avg SFR (EnvCount)')

        if i_z == 0: 
            param = sfms.get_bestfit_envcount_sfms()
            subs[i_z].plot(mass, param[0]*(np.array(mass)-10.5) + param[1], 
                    c=pretty_colors[2])

            # Plot average SFR for SF Group Catalog 
            mass, avg_sfrs, var_sfrs = [], [], [] 
            for i_mass in range(len(mass_bin.mass_low)): 

                avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z_groupcat(
                        mass_bin.mass_mid[i_mass], Mrcut=Mrcut)

                if ngal < 10: 
                    continue    # skip mass bin if there aren't many galaxies
                if mass_bin.mass_mid[i_mass] < 9.5: 
                    continue    # skip low mass

                mass.append(mass_bin.mass_mid[i_mass]) 
                avg_sfrs.append(avg_sfr)
                var_sfrs.append(var_sfr)

            subs[i_z].errorbar(mass, avg_sfrs, yerr=var_sfrs, 
                    lw=4, c=pretty_colors[3], label='Avg SFR (GroupCat)')

        #print str(zbin[0]) + ' < ' + str(zmids[i_z]) + ' < ' + str(zbin[1])
        subs[i_z].plot(mass, slopes[i_z]*(np.array(mass)-10.5) + yints[i_z], c=pretty_colors[2])
    
        bestfit_sfr = [] 
        for m in mass: 
            blah = util.get_sfr_mstar_z_bestfit(m, 0.5*(zbin[0]+zbin[1]))
            bestfit_sfr.append(blah[0]) 

        subs[i_z].plot(mass, bestfit_sfr, c=pretty_colors[3])
        
        subs[i_z].text(10.0, 1.25, '$\mathtt{z \sim '+str(0.5*(zbin[0]+zbin[1]))+'}$') 
        subs[i_z].text(10.75, -0.75,'slope= '+str('%.2f' % slopes[i_z])) 
        subs[i_z].text(10.75, -1.0, 'y-int= '+str('%.2f' % yints[i_z])) 

        subs[i_z].set_xlim([9.5, 12.0]) 
        subs[i_z].set_ylim([-1.5, 1.5]) 
        if i_z in (1, 2):
            subs[i_z].set_xlabel('log(M)') 
        if i_z == 0: 
            subs[i_z].set_ylabel('log(SFR)') 
        #if i_z == 3: 
        #    subs[i_z].legend()
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'sf_ms_data_', str("%.2f" % lowz_slope), '_', str("%.2f" % lowz_yint), '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear() 

def plot_sfms_groupcat(Mrcut=18):
    ''' Plot StarForming Main Sequence (mass vs SFR) from Star-Forming SDSS group catalog  

    Parameters
    ----------
    None

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 

    fig = plt.figure(1)#, figsize=[7,5])
    subs = fig.add_subplot(111) 

    # SDSS group catalog best fit 
    #groupcat_slope, groupcat_yint = sfms.sdss_groupcat_sfms_bestfit()
    #subs[0].plot(
    #        np.array(mass_bin.mass_mid), 
    #        (groupcat_slope * (np.array(mass_bin.mass_mid)-10.5)) + groupcat_yint, 
    #        c='k', lw=6, ls='--') 

    #zmids, slopes, yints = sfms.get_sfmsfit_sfr(lowz_slope, lowz_yint, clobber=True) 

    mass, avg_sfrs, var_sfrs = [], [], [] 
    for i_mass in range(len(mass_bin.mass_low)): 
        avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z_groupcat(mass_bin.mass_mid[i_mass], 
                Mrcut=Mrcut) 

        if ngal < 10: 
            continue    # skip mass bin if there aren't many galaxies
        if mass_bin.mass_mid[i_mass] < 9.5: 
            continue    # skip low mass

        mass.append(mass_bin.mass_mid[i_mass]) 
        avg_sfrs.append(avg_sfr)
        var_sfrs.append(var_sfr)
        
    sf_cen = sfms.sf_centrals(Mrcut=Mrcut) 
    subs.hist2d(sf_cen.mass, sf_cen.sfr, bins=[30, 100])
    #subs[i_z].scatter(centrals.mass, centrals.sfr, s=2, color=pretty_colors[3]) 
    #centrals = cq_group.central_catalog(Mrcut=19) 
    #subs[i_z].scatter(centrals.mass, centrals.sfr, s=2, color=pretty_colors[4]) 
    #centrals = cq_group.central_catalog(Mrcut=20) 
    #subs[i_z].scatter(centrals.mass, centrals.sfr, s=2, color=pretty_colors[5]) 

    subs.errorbar(mass, avg_sfrs, yerr=var_sfrs, 
            lw=4, c=pretty_colors[1], label='Average SFR')

    bestfit_params = sfms.get_bestfit_groupcat_sfms(Mrcut=Mrcut) 
    subs.plot(mass, util.line(np.array(mass)-10.5, bestfit_params), 
            lw=4, ls='--', c=pretty_colors[2], label='MPfit line') 
    subs.text(11.25, 0.0, 'Slope = '+str("%.2f" % bestfit_params[0]))
    subs.text(11.25, -0.5, 'Y-int = '+str("%.2f" % bestfit_params[1]))

    subs.set_xlim([9.5, 12.0]) 
    subs.set_ylim([-1.5, 1.5]) 
    subs.set_xlabel('log(M)') 
    subs.set_ylabel('log(SFR)') 
    subs.legend(loc='upper right') 
    
    fig_name = ''.join(['/home/users/hahn/research/figures/tinker/',
        'sf_ms_groupcat_', str(Mrcut), '.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear()

def plot_q_groupcat(Mrcut=18): 
    ''' Plot mass vs sSFR for the Quiescet SDSS group catalog  

    Parameters
    ----------
    Mrcut : Absolute magnitude cut-off that specifies the group catalog

    Notes
    -----
    * Mainly used for determing the final sSFR of the quiescent population 
        * Because the SFR for the quiescent populations are just upper limits, this is just to make sure that the offset with the quiescent mode in the SSFR mode does not hinder the fitting procedure. 

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 

    fig = plt.figure(1)
    subs = fig.add_subplot(111) 
        
    q_cen = sfms.q_centrals(Mrcut=Mrcut) # Quiescent centrals  
    subs.hist2d(q_cen.mass, q_cen.ssfr, bins=[30, 1000])  # 2D histogram 
    
    med_ssfrs = [] 
    masses = np.arange(9.75, 11.5, 0.25) 
    for mass in masses: 
        med_ssfr, var, ngal = sfms.get_ssfr_mstar_qgroupcat(mass) 
        med_ssfrs.append(med_ssfr) 

    subs.plot(masses, med_ssfrs, 
            lw=4, c=pretty_colors[1]) 

    fit_line_param = sfms.get_bestfit_qgroupcat_ssfr(Mrcut=Mrcut, clobber=True) 
    subs.plot(masses, fit_line_param[0].item() * (masses - 10.5) + fit_line_param[1].item(), 
            lw=4, ls='--', c=pretty_colors[3]) 

    subs.set_xlim([9.5, 12.0]) 
    subs.set_ylim([-15.0, -10.0]) 
    subs.set_xlabel('log(M)') 
    subs.set_ylabel('log(sSFR)') 
    
    fig_name = ''.join(['/home/users/hahn/research/figures/tinker/',
        'ssfr_mass_q_groupcat_mrcut', str(Mrcut), '.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear()

if __name__=='__main__': 
    #plot_group_cat_bigauss_bestfit()
    #plot_sdss_group_cat() 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau='instant') 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau='linear') 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau='constant') 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau=[0.6, 0.4, 0.1, 0.0], 
    #        sfms_slope=0.3, sfms_yint=-0.1) 
    #plot_cenque_ssfr_dist_evolution(Mrcut=18, fq='wetzel', tau='linefit', tau_param=[-0.5, 0.4], 
    #        sfms_slope=0.7, sfms_yint=0.125) 
    #plot_cenque_ssfr_dist_evolution(Mrcut=19, fq='wetzel', tau='linefit', tau_param=[-0.5, 0.4], 
    #        sfms_slope=0.7, sfms_yint=0.125) 
    #plot_cenque_ssfr_dist_evolution(Mrcut=20, fq='wetzel', tau='linefit', tau_param=[-0.5, 0.4], 
    #        sfms_slope=0.7, sfms_yint=0.125) 
    #plot_q_groupcat(Mrcut=18)
    plot_cenque_ssfr_dist_evolution(nsnaps=[1], fq='wetzel', tau='instant') #tau='linefit', tau_param=[-0.5, 0.4])

    #plot_sfms_data(0.51, -0.22)

    #plot_sfms_groupcat(Mrcut=18)
    #plot_sfms_groupcat(Mrcut=19)
    #plot_sfms_groupcat(Mrcut=20)
    #plot_cenque_sf_mainseq()

    #fq_fig = plot_fq_evol_w_geha() 
    #fq_fig = plot_fq_evol()
    #fq_fig.savefig('/home/users/hahn/research/figures/tinker/fq_evol_fig_lit.png', bbox_inches='tight')
    #
    #tau_fig = plot_quenching_efold() 
    #tau_fig.savefig('quenching_efold_fig.png', bbox_inches='tight')
