'''

Plotting code for CentralQuenching (CenQue) project

Author(s): ChangHoon Hahn 

'''

import numpy as np
import matplotlib.pyplot as plt
import mpfit

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
        fig = plt.figure(1, figsize=(20,8))
        fig.subplots_adjust(hspace=0., wspace=0.)
        new_plot = True
    else: 
        new_plot = False

    # mass bins of the panel        ( hardcoded ) 
    panel_mass_bins = [
            [10.0, 10.5], [10.5, 11.0], [11.0, 11.5]
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

        sub = fig.add_subplot(1, 3, i_mass+1)       # panel subplot 
        sub.plot(ssfr_bin_mid, ssfr_hist, 
                color=line_color, lw=4, ls=line_style, label=ssfr_hist_label) 
        
        if new_plot == True:        # put mass labels in panels 
            plt.text(-10.5, 1.4, 
                    r'$\mathtt{log \; M_{*} = ['+str(panel_mass[0])+', '+str(panel_mass[1])+']}$', 
                    fontsize=24)

        # set axes limits
        sub.set_xlim([-13.0, -7.0])
        sub.set_ylim([0.0, 1.6])
        # set y-axes labels
        if i_mass == 0: 
            sub.set_ylabel(r'$P(log \; SSFR)$') 
        elif i_mass == 1:
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
    snap = cq.CenQue() 
    snap.readin(nsnap=13, file_type='sf assign', **kwargs) 
    ssfr_fig = plot_cenque_ssfr_dist(snap)

    for i_nsnap in [1]:#reversed(range(1,13)):
        next_snap = cq.CenQue() 
        next_snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
        
        ssfr_fig = plot_cenque_ssfr_dist(next_snap, fig=ssfr_fig) 

    central_ssfr = cq_group.central_catalog(Mrcut=Mrcut, clobber=True) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, label= 'Mrcut = '+str(Mrcut)) 
    
    if isinstance(kwargs['tau'], list) == False: 
        fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
            'cenque_ssfr_evol_', kwargs['tau'], 'tau_', 
            kwargs['fq'], 'fq_Mrcut', str(Mrcut),'.png']) 
    else: 
        fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
            'cenque_ssfr_evol_', '_'.join([str(t) for t in kwargs['tau']]), 'tau_', 
            kwargs['fq'], 'fq_Mrcut', str(Mrcut),'.png']) 
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

        sub.plot(ssfr_bin_mid, cq_group.double_gaussian(ssfr_bin_mid, (output[i_mass])[2]), 
                color='blue', lw=4, ls='-', label='Best Fit') 
        sub.text(-11, 0.7, 'Best fit fQ = '+str(((output[i_mass])[2])[0]))
        sub.text(-11, 0.65, 'Best fit SSFR Q = '+str(((output[i_mass])[2])[1]))
        sub.text(-11, 0.6, 'Best fit SSFR SF = '+str(((output[i_mass])[2])[2]))

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
def plot_sfms_data(): 
    ''' plot SF MS from data
    '''

    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 
    zbins = [ (0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8)]  # zbin

    fig, subs = plt.subplots(1, len(zbins), figsize=[20, 5]) 
    subs = subs.ravel() 

    # SDSS group catalog best fit 
    groupcat_slope, groupcat_yint = sfms.sdss_groupcat_sfms_bestfit()
    subs[0].plot(
            np.array(mass_bin.mass_mid), 
            (groupcat_slope * (np.array(mass_bin.mass_mid)-10.5)) + groupcat_yint, 
            c='k', lw=6, ls='--') 

    fits = [] 
    for i_z, zbin in enumerate(zbins): 
        
        mass = [] 
        avg_sfrs = []
        var_sfrs = [] 

        for i_mass in range(len(mass_bin.mass_low)): 

            avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z(mass_bin.mass_mid[i_mass], 
                    0.5*(zbin[0] + zbin[1]))

            if ngal < 10: 
                continue    # skip mass bin if there aren't many galaxies
            if mass_bin.mass_mid[i_mass] < 9.5: 
                continue    # skip low mass

            mass.append(mass_bin.mass_mid[i_mass]) 
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)
            
        subs[i_z].errorbar(mass, avg_sfrs, yerr=var_sfrs, c=pretty_colors[1])

        p0 = [0.0+i_z*0.2] 
        fa = {'x': np.array(mass)-10.5, 'y': np.array(avg_sfrs), 'err': np.array(var_sfrs)} 

        bestfit = mpfit.mpfit(sfms.mpfit_line_fixedslope, p0, functkw=fa, nprint=0)
                
        print bestfit.params
        subs[i_z].plot(mass, sfms.line_fixedslope(np.array(mass)-10.5, bestfit.params), c=pretty_colors[2]) 
        
        subs[i_z].text(10.0, 1.0, '$\mathtt{z \sim '+str(0.5*(zbin[0]+zbin[1]))+'}$') 
        subs[i_z].set_xlim([9.5, 12.0]) 
        subs[i_z].set_ylim([-0.5, 1.5]) 
        if i_z in (1, 2):
            subs[i_z].set_xlabel('log(M)') 
        if i_z == 0: 
            subs[i_z].set_ylabel('log(SFR)') 

    return fig 

if __name__=='__main__': 
    #plot_group_cat_bigauss_bestfit()
    #plot_sdss_group_cat() 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau='instant') 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau='linear') 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau='constant') 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau=[0.9, 0.8, 0.4, 0.3]) 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau=[0.9, 0.8, 0.4, 0.0]) 
    #plot_cenque_ssfr_dist_evolution(fq='wetzel', tau=[0.7, 0.4, 0.4, 0.1]) 
    fig = plot_sfms_data()
    fig.savefig('/home/users/hahn/research/figures/tinker/sf_ms_data.png', bbox_inches='tight')
    #plot_cenque_sf_mainseq()

    #fq_fig = plot_fq_evol_w_geha() 
    #fq_fig = plot_fq_evol()
    #fq_fig.savefig('/home/users/hahn/research/figures/tinker/fq_evol_fig_lit.png', bbox_inches='tight')
    #
    #tau_fig = plot_quenching_efold() 
    #tau_fig.savefig('quenching_efold_fig.png', bbox_inches='tight')
