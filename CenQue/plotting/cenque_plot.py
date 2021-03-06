'''

Plotting code for CentralQuenching (CenQue) project

Author(s): ChangHoon Hahn 

'''

import numpy as np
import matplotlib.pyplot as plt
import mpfit
import scipy.stats as scistat
import pylab as pyl

#----- Local -----
from utility.plotting import prettyplot
from utility.plotting import prettycolors 
import cenque_utility as util 
import cenque as cq 
import cenque_groupcat as cq_group
import sf_mainseq as sfms 
import bovy_plot as bovy

# ssfr distribution ------------------------------------------
def plot_cenque_ssfr_dist(cenque, fig=None, **kwargs): 
    ''' Plot CenQue data
    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
    
    if fig == None: 
        fig = plt.figure(1, figsize=(16,16))
        #fig.subplots_adjust(hspace=0., wspace=0.)
        new_plot = True
    else: 
        new_plot = False

    # mass bins of the panel        ( hardcoded ) 
    panel_mass_bins = [ 
            [9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]
            ]
    
    for i_mass, panel_mass in enumerate(panel_mass_bins):       # loop through each panel 

        mass_limit = (cenque.mass >= panel_mass[0]) & (cenque.mass < panel_mass[1])

        # SSFR histogram 
        ssfr_hist, ssfr_bin_edges = np.histogram(cenque.ssfr[mass_limit], 
                range=[-13.0, -7], bins=40, normed=True)
        ngal_bin = len(cenque.ssfr[mass_limit])

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
                ssfr_hist_label = r'$\mathtt{z='+str(cenque.zsnap)+'}$'     # redshift lable 
    
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

        if 'lw' in kwargs: 
            line_width = kwargs['lw'] 
        else:
            line_width = 4

        #sub = fig.add_subplot(1, 5, i_mass+1)       # panel subplot 
        sub = fig.add_subplot(2, 2, i_mass+1)       # panel subplot 
        sub.plot(ssfr_bin_mid, ssfr_hist, 
                color=line_color, lw=line_width, ls=line_style, label=ssfr_hist_label) 
   
        '''
        if cenque.zsnap: 
            ssfr_cut = -11.15 + 0.76*(cenque.zsnap - 0.05) - 0.35 * (panel_mass[0] - 10.5)
            sub.vlines( ssfr_cut, 0.0, 2.0, colors='black', linestyles='dashed')
            ssfr_cut = -11.15 + 0.76*(cenque.zsnap - 0.05) - 0.35 * ( 0.5*(panel_mass[0]+panel_mass[1]) - 10.5)
            sub.vlines( ssfr_cut, 0.0, 2.0, colors='black', linestyles='solid')
            ssfr_cut = -11.15 + 0.76*(cenque.zsnap - 0.05) - 0.35 * (panel_mass[1] - 10.5)
            sub.vlines( ssfr_cut, 0.0, 2.0, colors='black', linestyles='dashed')
        '''
        #sub.text(-12.0, 1.2, r'$N_{gal}='+str(ngal_bin)+'$')
        if new_plot == True:        # put mass labels in panels 
            plt.text(-11.25, 1.4, 
                    r'$\mathtt{log \; M_{*} = ['+str(panel_mass[0])+', '+str(panel_mass[1])+']}$', 
                    fontsize=24)

        # set axes limits
        sub.set_xlim([-13.0, -7.0])
        sub.set_ylim([0.0, 1.6])
        # set y-axes labels
        if i_mass == 0: 
            sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$') 
        elif i_mass == 1: 
            sub.set_yticklabels([])
        elif i_mass == 2:
            #sub.set_yticklabels([])
            sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$') 
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$') 
                
            #try: 
            #    fig_leg.remove() 
            #except UnboundLocalError: 
            #    pass
            #
            #fig_leg = sub.legend(loc='lower left', prop={'size':'28'})        # legends 
    ''' 
    if 'tau' in kwargs.keys(): 
        sub = fig.add_subplot(1, 5, 5)       # tau panel 

        mass_bin = util.simple_mass_bin()       # mass bin 
        
        # satellite quenching fraction 
        tau_mass = util.get_quenching_efold(np.array(mass_bin.mass_mid), 
                type='linear') 

        sub.plot(mass_bin.mass_mid, tau_mass, color='black', lw=4, ls='--',
                label=r'$\tau_\mathtt{satellite}$') 
        
        # central quenching fraction 
        tau_mass = util.get_quenching_efold(np.array(mass_bin.mass_mid), 
                type=kwargs['tau'], param=kwargs['tau_param']) 

        sub.plot(mass_bin.mass_mid, tau_mass, color=pretty_colors[5], lw=4) 
        
        sub.set_xlim([9.0, 12.0])
        sub.set_ylim([0.05, 2.0])
        sub.set_yscale('log') 
        sub.set_xlabel('Mass') 
        sub.set_ylabel(r'Quenching e-Fold time $\tau$') 
        sub.yaxis.set_label_position('right')
        sub.legend(loc='lower left', prop={'size':'24'}) 
    '''
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
    ssfr_fig = plot_cenque_ssfr_dist(snap, lw=2, line_style='--')      # plot!
    
    # determine which snapshots to plot 
    if 'nsnaps' in kwargs.keys():   
        nsnaps = kwargs['nsnaps']
    else: 
        nsnaps = [1] 

    for i_nsnap in nsnaps:  
        next_snap = cq.CenQue() 
        next_snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
        
        #ssfr_fig = plot_cenque_ssfr_dist(next_snap, fig=ssfr_fig) 
        ssfr_fig = plot_cenque_ssfr_dist(next_snap)#, fig=ssfr_fig) 
    
    # overplot SDSS group catalog sSFR dist
    central_ssfr = cq_group.central_catalog(Mrcut=Mrcut) 
    #ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, label= r'$M_\mathtt{r,cut} = '+str(Mrcut)+'$', **kwargs) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, label= r'SDSS Group Catalog', **kwargs) 
    
    # file name ----------------------------------------------------------------------------
    if 'sfms_slope' in kwargs.keys():       # sfms specifier
        slope_str = str("%.2f" % kwargs['sfms_slope']) 
        yint_str = str("%.2f" % kwargs['sfms_yint']) 
        sfms_str = '_sfms_slope'+slope_str+'_yint'+yint_str
    else: 
        sfms_str = ''

    # Quenching Fraction specifier 
    if 'fqing_slope' in kwargs.keys(): 
        fqing_slope_str = str(kwargs['fqing_slope'])
    else: 
        fqing_slope_str = str(0.63)

    if 'fqing_yint' in kwargs.keys(): 
        fqing_yint_str = str(kwargs['fqing_yint'])
    else: 
        fqing_yint_str = str(-6.04) 

    fqing_str = ''.join([fqing_slope_str, '_', fqing_yint_str, 'fqing']) 

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # Stellar mass specifier 
    if kwargs['stellmass'].lower() == 'integrated': 
        mass_str = '_integ'
    elif kwargs['stellmass'].lower() == 'sham': 
        mass_str = '_sham'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 

    fig_file = ''.join(['figure/', 
        'cenque_ssfr_evol_', tau_str, '_', fqing_str, '_Mrcut', str(Mrcut), 
        file_type_str, '.png'])
    ssfr_fig.savefig(fig_file, bbox_inches='tight') 
    ssfr_fig.clear() 
    plt.close(ssfr_fig)

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
    
    fig_file = ''.join(['figure/tinker/', 
        'cenque_ssfr_evol_', kwargs['tau'], 'tau_', kwargs['fq'], 'fq_Mrcut', str(Mrcut),
        '_match2isedfit.png']) 
    ssfr_fig.savefig(fig_file, bbox_inches='tight') 
    ssfr_fig.clear() 

def plot_cenque_quenching_ssfr_dist(nsnap, **kwargs): 
    ''' Plot sSFR distribution of snapshot nsnap, with 'quenching' population highlighted  

    Parameters
    ----------
    nsnap : snapshot #

    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
    
    fig = plt.figure(figsize=(25,8))
    fig.subplots_adjust(hspace=0., wspace=0.)

    # mass bins of the panel        ( hardcoded ) 
    panel_mass_bins = [ 
            [9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]
            ]
    
    # import snapshot 
    cenque = cq.CenQue() 
    cenque.readin(nsnap=nsnap, file_type='evol from 13', **kwargs) 

    for i_mass, panel_mass in enumerate(panel_mass_bins):       # loop through each panel 

        mass_limit = (cenque.mass >= panel_mass[0]) & (cenque.mass < panel_mass[1])
    
        quenching_tau = (cenque.tau > 0.0) 
        else_tau = (cenque.tau < 0.0) 
        
        sub = fig.add_subplot(1,4, i_mass+1)

        if np.sum(mass_limit & quenching_tau) == 0: 
            # stacked histogram 
            sub.hist(cenque.ssfr[mass_limit & else_tau], 25, normed=True) 
        else: 
            # stacked histogram 
            sub.hist([cenque.ssfr[mass_limit & quenching_tau], cenque.ssfr[mass_limit & else_tau]], 
                    25, stacked=True, normed=True) 

        sub.text(-10.75, 1.2, r'$N_{quenching}/N_{bin} = '+\
                str('%.2f' % (np.float(len(cenque.ssfr[mass_limit & quenching_tau]))/np.float(len(cenque.ssfr[mass_limit]))))+'$') 

        if 'label' in kwargs: 
            ssfr_hist_label = kwargs['label']
        else: 
            try: 
                cenque.zsnap
            except AttributeError: 
                ssfr_hist_label = 'Centrals'     # redshift lable 
            else: 
                ssfr_hist_label = r'$z='+str(cenque.zsnap)+'$'     # redshift lable 
    
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
                
            fig_leg = sub.legend(loc='lower right', prop={'size':'12'})        # legends 

    # sfms specifier
    if 'sfms_slope' in kwargs.keys(): 
        slope_str = str("%.2f" % kwargs['sfms_slope']) 
        yint_str = str("%.2f" % kwargs['sfms_yint']) 
        sfms_str = '_sfms_slope'+slope_str+'_yint'+yint_str
    else: 
        sfms_str = ''

    # Quenching Fraction specifier 
    if 'fqing_slope' in kwargs.keys(): 
        fqing_slope_str = str(kwargs['fqing_slope'])
    else: 
        fqing_slope_str = str(0.63)

    if 'fqing_yint' in kwargs.keys(): 
        fqing_yint_str = str(kwargs['fqing_yint'])
    else: 
        fqing_yint_str = str(-6.04) 

    fqing_str = ''.join([fqing_slope_str, '_', fqing_yint_str, 'fqing']) 

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
   
    # Stellar mass specifier 
    if kwargs['stellmass'].lower() == 'integrated': 
        mass_str = '_integ'
    elif kwargs['stellmass'].lower() == 'sham': 
        mass_str = '_sham'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    fig_file = ''.join(['figure/', 
        'cenque_ssfr_evol_', tau_str, '_', fqing_str, file_type_str, '_nsnap', str(nsnap), 
        '_quenching_component.png'])
    fig.savefig(fig_file, bbox_inches='tight') 
    fig.clear() 
    plt.close(fig)

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

    fig_file = ''.join(['figure/tinker/', 
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
            [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]
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

    fig_file = ''.join(['figure/tinker/', 
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
    
    fig_file = ''.join(['figure/tinker/', 
        'sf_ms_evol.png']) 
    fig.savefig(fig_file) 

# quiescent fraction -----------------------------------
def plot_fq_evol(): 
    ''' Plot fq evolution comparison
    '''

    # snapshot redshifts
    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[zbin < 1.0]
    #zbin = [0.1*np.float(i) for i in range(1,10)]   # zbin 

    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 
    
    # load literature data 
    # modified tinker
    mod_tink_file = ''.join(['dat/central_quenching/literature/modified_tinker_fq.dat']) 
    mod_tink_mass, mod_tink_fq = np.loadtxt(mod_tink_file, unpack=True, usecols=[0,1])   

    fq_types = ['cosmosinterp', 'wetzel', 'wetzelsmooth'] 
    
    fig, subs = plt.subplots(1, len(fq_types), figsize=[len(fq_types)*5, 5]) 
    subs = subs.ravel() 

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

    subs[0].set_ylabel('Quiescent Fraction') 
    subs[0].legend(loc='upper left') 
    fig_name = ''.join(['figure/fq_evol_comp_', '_'.join(fq_types), '.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear() 

def plot_fq_geha_groupcat(Mrcut=18): 
    ''' Plot comparison of Modified Tinker Catalog fq from Geha and SDSS group catalog 

    Notes
    -----
    * Mrcut should only be = 18 because the modified tinker catalog is constructed from Mrcut=-18 with z <= 0.06. 
    * SDSS Group Catalog fQ is calculated based on SFR-M* cut. 
    * Geha Modified Tinker catalog fQ is calculated based on EW Halpha and Dn4000. 
    * Should they match?  

    '''
    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    if Mrcut != 18: 
        print 'Mrcut should only be = 18 because the modified tinker catalog is constructed from Mrcut=-18 with z <= 0.06.'
    
    # load literature data 
    # modified tinker
    mod_tink_file = ''.join(['dat/central_quenching/literature/modified_tinker_fq.dat']) 
    mod_tink_mass, mod_tink_fq = np.loadtxt(mod_tink_file, unpack=True, usecols=[0,1])   
    
    # read group catalog 
    central = cq_group.central_catalog(Mrcut=Mrcut) 
    central_z = central.z    # redshifts 
    median_z = np.median(central_z)
   
    fq_mass, fq = [], [] 

    mass_bins = util.simple_mass_bin() 
    for i_m in range(mass_bins.nbins): 

        mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
        mass_bin_mid = round(mass_bins.mass_mid[i_m], 2)
        mass_bin_high = round(mass_bins.mass_high[i_m], 2) 

        # boolean list for mass range
        mass_bin_bool = (central.mass > mass_bin_low) & (central.mass <= mass_bin_high)
        
        if np.sum(mass_bin_bool) < 10: 
            continue 

        sfqs = util.sfq_classify( central.mass[mass_bin_bool], central.sfr[mass_bin_bool], median_z ) 
        ngal = np.float(len(sfqs))
        ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
        fq_mass.append(mass_bin_mid) 
        fq.append(ngal_q/ngal) 

    fig = plt.figure(figsize=[8, 8]) 
    sub = fig.add_subplot(111)

    sub.plot(mod_tink_mass, mod_tink_fq, 
        color='black', lw=6, ls='--', label='Modified Tinker Group' ) 
    
    sub.plot(fq_mass, fq, 
            color='black', lw=6, label=r'SDSS Group Catalog; z = '+str(median_z))

    sub.set_xlim([9.0, 12.0])
    sub.set_ylim([0.0, 1.0])

    sub.set_xlabel('Mass') 
    sub.set_ylabel('Quiescent Fraction') 
    sub.legend(loc='upper left') 

    fig_name = ''.join(['figure/fq_evol_groupcat', str(Mrcut), '_geha_comp.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear() 

def plot_snapshot_fqobs_evol(nsnaps=[1,2,3,4,5,6,7,8,9,10,11,12],fq_type='wetzelsmooth', **kwargs): 
    ''' Plot the observed quiescent fraction of snapshots

    Parameters
    ----------
    nsnaps : list of snapshots
    fqtype : type of queiscent fraction  
    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    # snapshot redshifts
    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[nsnaps]
    #zbin = [0.1*np.float(i) for i in range(1,10)]   # zbin 

    mass_bins = util.simple_mass_bin()                    # mass bin 
       
    fig, subs = plt.subplots(1,2, figsize=[10, 5]) 
    subs = subs.ravel() 

    snap = cq.CenQue() 
    snap.readin(nsnap=13, file_type='sf assign', **kwargs)  # starting CenQue 

    fq_mass, fq = [], [] 
    for i_m in range(mass_bins.nbins): 
        mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
        mass_bin_mid = round(mass_bins.mass_mid[i_m], 2)
        mass_bin_high = round(mass_bins.mass_high[i_m], 2) 

        # boolean list for mass range
        mass_bin_bool = (snap.mass > mass_bin_low) & (snap.mass <= mass_bin_high)
        
        if np.sum(mass_bin_bool) == 0: 
            continue 

        sfqs = util.sfq_classify( snap.mass[mass_bin_bool], snap.sfr[mass_bin_bool], snap.zsnap ) 
        ngal = np.float(len(sfqs))
        ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
        fq_mass.append(mass_bin_mid) 
        fq.append(ngal_q/ngal) 
    
    subs[0].plot(fq_mass, fq, color='black', ls='--', lw=4, label='z= '+str(snap.zsnap)) 
    
    for i_nsnap in nsnaps: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 

        fq_mass, fq = [], [] 

        for i_m in range(mass_bins.nbins): 
            mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
            mass_bin_mid = round(mass_bins.mass_mid[i_m], 2)
            mass_bin_high = round(mass_bins.mass_high[i_m], 2) 

            # boolean list for mass range
            mass_bin_bool = (snap.mass > mass_bin_low) & (snap.mass <= mass_bin_high)
            
            if np.sum(mass_bin_bool) == 0: 
                continue 

            sfqs = util.sfq_classify( snap.mass[mass_bin_bool], snap.sfr[mass_bin_bool], snap.zsnap ) 
            ngal = np.float(len(sfqs))
            ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
            fq_mass.append(mass_bin_mid) 
            fq.append(ngal_q/ngal) 
        
        subs[0].plot(fq_mass, fq, color=pretty_colors[i_nsnap-1], lw=4, label='z= '+str(snap.zsnap)) 
    
    subs[0].set_title('Snapshots')
    subs[0].set_xlim([9.0, 12.0])
    subs[0].set_ylim([0.0, 1.0])
    subs[0].set_xlabel('Mass') 
    for i_z, z in enumerate(zbin): 

        # plot fq(Mass) 
        fq_mass = [util.get_fq(mass_bins.mass_mid[i], z, lit=fq_type) 
                for i in range(len(mass_bins.mass_mid))]
        subs[1].plot(mass_bins.mass_mid, fq_mass, 
                color=pretty_colors[i_z], lw=4, label='z = '+str(z) ) 
    
    subs[1].set_title(fq_type) 

    subs[1].set_xlim([9.0, 12.0])
    subs[1].set_ylim([0.0, 1.0])

    subs[1].set_xlabel('Mass') 

    subs[0].set_ylabel('Quiescent Fraction') 
    
    # file name ----------------------------------------------------------------
    # Quenching Fraction specifier 
    if 'fqing_slope' in kwargs.keys(): 
        fqing_slope_str = str(kwargs['fqing_slope'])
    else: 
        fqing_slope_str = str(0.63)

    if 'fqing_yint' in kwargs.keys(): 
        fqing_yint_str = str(kwargs['fqing_yint'])
    else: 
        fqing_yint_str = str(-6.04) 

    fqing_str = ''.join([fqing_slope_str, '_', fqing_yint_str, 'fqing']) 

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'

    # Stellar mass specifier 
    if kwargs['stellmass'].lower() == 'integrated': 
        mass_str = '_integ'
    elif kwargs['stellmass'].lower() == 'sham': 
        mass_str = '_sham'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 

    fig_file = ''.join(['figure/', 
        'fq_obs_snapshots_', tau_str, '_', fqing_str, file_type_str, '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear()
    plt.close(fig)

def plot_snapshot_fqobs(nsnap, fq_type='wetzel', **kwargs): 
    ''' Plot the observed quiescent fraction of a snapshot with parameterized fQ

    Parameters
    ----------
    nsnap : (int) snapshot #
    fqtype : type of queiscent fraction  

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    # snapshot redshifts
    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[nsnap]

    mass_bins = util.simple_mass_bin()                    # mass bin 
       
    fig = plt.figure(figsize=[5, 5]) 
    subs = fig.add_subplot(111)

    snap = cq.CenQue() 
    snap.readin(nsnap=nsnap, file_type='evol from 13', **kwargs) 

    fq_mass, fq = [], [] 

    for i_m in range(mass_bins.nbins): 
        mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
        mass_bin_mid = round(mass_bins.mass_mid[i_m], 2)
        mass_bin_high = round(mass_bins.mass_high[i_m], 2) 

        # boolean list for mass range
        mass_bin_bool = (snap.mass > mass_bin_low) & (snap.mass <= mass_bin_high)
        
        if np.sum(mass_bin_bool) == 0: 
            continue 

        sfqs = util.sfq_classify( snap.mass[mass_bin_bool], snap.sfr[mass_bin_bool], snap.zsnap ) 
        ngal = np.float(len(sfqs))
        ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
        fq_mass.append(mass_bin_mid) 
        fq.append(ngal_q/ngal) 
    
    subs.plot(fq_mass, fq, color=pretty_colors[3], lw=4, label='Snapshot '+str(snap.nsnap)) 
    
    # parameterized fq 
    fq_mass = [util.get_fq(mass_bins.mass_mid[i], zbin, lit=fq_type) for i in range(len(mass_bins.mass_mid))]
    subs.plot(mass_bins.mass_mid, fq_mass, 
            color='black', lw=4, ls='--', label='Wetzel; z = '+str(zbin) ) 

    subs.set_xlim([9.0, 12.0])
    subs.set_ylim([0.0, 1.0])
    subs.set_xlabel('Mass') 
    subs.set_ylabel('Quiescent Fraction') 
    subs.legend(loc='upper left') 

    # Quenching Fraction specifier 
    if 'fqing_slope' in kwargs.keys(): 
        fqing_slope_str = str(kwargs['fqing_slope'])
    else: 
        fqing_slope_str = str(0.63)

    if 'fqing_yint' in kwargs.keys(): 
        fqing_yint_str = str(kwargs['fqing_yint'])
    else: 
        fqing_yint_str = str(-6.04) 

    fqing_str = ''.join([fqing_slope_str, '_', fqing_yint_str, 'fqing']) 

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # Stellar mass specifier 
    if kwargs['stellmass'].lower() == 'integrated': 
        mass_str = '_integ'
    elif kwargs['stellmass'].lower() == 'sham': 
        mass_str = '_sham'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 

    fig_file = ''.join(['figure/', 
        'fq_obs_snapshot', str(nsnap), '_', tau_str, '_', fqing_str, file_type_str, '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear()
    plt.close(fig)

def plot_quenching_efold(taus, tau_params): 
    ''' Plot quenching e-fold (tau) as a function of M*

    Parameters
    ----------
    taus : list of tau types 
    tau_params : list of tau paramters

    Notes
    -----
    * taus and tau_params have to be the same length 

    '''
    if len(taus) != len(tau_params): 
        return 'asdfasdaf'

    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 
       
    #fig, subs = plt.subplots(1, len(taus), figsize=[5*len(taus), 5]) 
    #subs = subs.ravel() 
    fig = plt.figure(1, figsize=[8,8])
    subs = fig.add_subplot(111)
    
    tau_str = ''
    for i_tau, tau in enumerate(taus): 
        

        tau_mass = util.get_quenching_efold(np.array(mass_bin.mass_mid), 
                type=tau, param=tau_params[i_tau]) 

        subs.plot(10**np.array(mass_bin.mass_mid), tau_mass, 
                color=pretty_colors[3], lw=4, label='Central') 
        
        # satellite quenching fraction 
        sat_mass = np.array([9.75 + 0.2*np.float(i) for i in range(9)])
        tau_mass = util.get_quenching_efold(sat_mass, 
                type='linear') 
    
        subs.plot(10**sat_mass, tau_mass, color='black', lw=4, ls='--',
                label=r'Satellite (Wetzel et al. 2013)') 

        #subs.set_title(tau) 
        subs.set_xlim([10**9.5, 10**11.75])
        subs.set_xscale('log')
        subs.set_ylim([0.0, 2.5])
        subs.set_xlabel(r'Stellar Mass $[M_\odot]}$',fontsize=28) 
        subs.legend(loc='upper right') 

        # tau specifier 
        if tau == 'discrete': 
            tau_str += '_'+'_'.join( [str("%.1f" % t) for t in tau_params[i_tau]] )
        elif tau == 'linefit':
            tau_str += '_line'+'_'.join( [str("%.2f" % t) for t in tau_params[i_tau]] )
        else: 
            tau_str += '_'+tau

        tau_str += 'tau'

    subs.set_ylabel(r'Quenching Timescales $(\tau_\mathtt{Q})$ [Gyr]',fontsize=28) 
                
    fig.savefig('figure/quenching_efold'+tau_str+'.png', bbox_inches='tight')
    fig.clear() 

def plot_group_cat_bigauss_bestfit(): 
    plot_sdss_group_cat_bestfit(Mrcut=18)
    plot_sdss_group_cat_bestfit(Mrcut=19)
    plot_sdss_group_cat_bestfit(Mrcut=20)

# Group Catalog ------------------------------
def plot_groupcat_zdist(): 
    ''' Plot redshift distribution of all three group catalogs 
    
    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 

    fig = plt.figure(1)
    sub = fig.add_subplot(111)
    
    # loop through redshift distribution
    for i_mr, mrcut in enumerate([18, 19, 20]): 

        # read group catalog 
        central = cq_group.central_catalog(Mrcut=mrcut) 
        central_z = central.z    # redshifts 
        median_z = np.median(central_z)
        
        z_hist, z_bin_edges = np.histogram(central_z, range=[0.0, 0.2], bins=40, normed=True)
        
        z_bin_low = z_bin_edges[:-1]
        z_bin_high = z_bin_edges[1:]
        z_bin_mid = [ 0.5*(z_bin_low[i] + z_bin_high[i]) for i in range(len(z_bin_low)) ] 

        sub.plot(z_bin_mid, z_hist, 
                lw=4, c= pretty_colors[i_mr], label=r'$\mathtt{M_r = -'+str(mrcut)+'}$')
        sub.text(median_z, 50-i_mr*5, 'Median z ='+str('%.2f' % median_z)) 

    sub.set_xlim([0.0, 0.15]) 
    sub.set_xlabel(r'$z$ (Redshift)') 
    sub.legend(loc='upper right') 

    fig_file = 'figure/groupcat_zdist.png'
    fig.savefig(fig_file, bbox_inches='tight')

def plot_groupcat_obs_fq(Mrcut=18):
    ''' Plot observed f_q using the SDSS Group Catalog and sSFR SF/Q classification 

    Parameters
    ----------
    Mrcut : Absolute magnitude cutoff specifier for group catalog 

    Notes
    ----- 
    * This is to make sure that the SF/Q classification is consistent with the Andrew Wetzel's f_Q parameterization 
    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
    
    # read group catalog 
    central = cq_group.central_catalog(Mrcut=Mrcut) 
    central_z = central.z    # redshifts 
    median_z = np.median(central_z)
   
    fq_mass, fq = [], [] 

    mass_bins = util.simple_mass_bin() 
    for i_m in range(mass_bins.nbins): 

        mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
        mass_bin_mid = round(mass_bins.mass_mid[i_m], 2)
        mass_bin_high = round(mass_bins.mass_high[i_m], 2) 

        # boolean list for mass range
        mass_bin_bool = (central.mass > mass_bin_low) & (central.mass <= mass_bin_high)
        
        if np.sum(mass_bin_bool) < 10: 
            continue 

        sfqs = util.sfq_classify( central.mass[mass_bin_bool], central.sfr[mass_bin_bool], median_z ) 
        ngal = np.float(len(sfqs))
        ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
        fq_mass.append(mass_bin_mid) 
        fq.append(ngal_q/ngal) 
    
    fig = plt.figure(1, figsize=[8,8])
    sub = fig.add_subplot(111)

    # snapshot redshifts
    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    close_z = min(range(len(zbin)), key=lambda i: abs(zbin[i] - median_z))
    zbin = [ zbin[close_z - 1], zbin[close_z], zbin[close_z + 1] ]

    for i_z, z in enumerate(zbin): 

        # plot fq(Mass) 
        fq_fit = [util.get_fq(mass_bins.mass_mid[i], z, lit='wetzelsmooth') 
                for i in range(len(mass_bins.mass_mid))]
        sub.plot(mass_bins.mass_mid, fq_fit, 
                color=pretty_colors[i_z], lw=4, ls='--', label='z = '+str(z) ) 
    
    sub.plot(fq_mass, fq, color='black', lw=6, label=r'SDSS Group Catalog; z = '+str(median_z))
    sub.set_xlim([9.0, 12.0])
    sub.set_ylim([0.0, 1.0])
    sub.set_xlabel('Mass') 
    sub.set_ylabel('Quiescent Fraction') 
    sub.legend(loc='upper left') 

    fig_name = 'figure/fq_sdss_groupcat_mrcut'+str(Mrcut)+'_sfrclass.png'
    fig.savefig(fig_name, bbox_inches='tight') 
    fig.clear() 

# Mhalo-M* ------------------------------------
def plot_mhalo_mstar(i_nsnap = 1, **kwargs): 
    ''' Plot Mhalo versus M* of galaxies for CenQue file 

    '''
    prettyplot() 
    pretty_colors = prettycolors() 

    # import evolved snapshot 
    if i_nsnap == 13: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='sf assign', **kwargs) 
    else: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 

    mhalo = snap.halo_mass 
    mstar = snap.mass
    
    bovy.scatterplot(mstar, mhalo, scatter=True, color=pretty_colors[1], s=3, 
            xrange=[9.0, 11.5], 
            xlabel='\mathtt{M_*}', ylabel='\mathtt{M_{Halo}}')

    # figure file name ------------------------------------------------------------
    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # Stellar mass specifier 
    if kwargs['stellmass'].lower() == 'integrated': 
        mass_str = '_integ'
    elif kwargs['stellmass'].lower() == 'sham': 
        mass_str = '_sham'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 
    
    fig_name1 = ''.join(['figure/', 
        'cenque_mstar_mhalo_snapshot', str(i_nsnap), tau_str, file_type_str, 
        '_scatter.png'])
    plt.savefig(fig_name1, bbox_inches='tight')

def plot_mhalo_mstar_sham_integrated(i_nsnap = 1, **kwargs): 
    ''' Plot Mhalo versus M* of galaxies for CenQue file 

    '''
    prettyplot() 
    pretty_colors = prettycolors()  

    fig1 = plt.figure(1, figsize=[8,8])
    sub1 = fig1.add_subplot(111)
    
    mass_bins = util.simple_mass_bin()  # use simplest mass bins

    for i_mstar, mstar_type in enumerate(['sham', 'integrated']): 

        # import evolved snapshot 
        if i_nsnap == 13: 
            snap = cq.CenQue() 
            snap.readin(nsnap=i_nsnap, file_type='sf assign', **kwargs) 
        else: 
            kwargs['stellmass'] = mstar_type
            snap = cq.CenQue() 
            snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 

        mhalo = snap.halo_mass 
        mstar = snap.mass
    
        #sub1.scatter(mstar, mhalo, c=pretty_colors[i_mstar+1], s=3) 

        avg_mhalo, var_mhalo = [], [] 
        for i_m in range(mass_bins.nbins): 
            mass_limit = (mstar > mass_bins.mass_low[i_m]) & \
                    (mstar <= mass_bins.mass_high[i_m]) 
            avg_mhalo.append( np.mean(mhalo[mass_limit]) ) 
            var_mhalo.append( np.std(mhalo[mass_limit]) ) 
            
            #print mass_bins.mass_low[i_m], ' - ', mass_bins.mass_high[i_m]
            #print min(mhalo[mass_limit]), max(mhalo[mass_limit])
            #print np.mean(mhalo[mass_limit]), np.std(mhalo[mass_limit])

        if mstar_type == 'sham':  
            sub1.errorbar(mass_bins.mass_mid, avg_mhalo, yerr=var_mhalo, 
                        lw=4, c=pretty_colors[1], label='SHAM')
        else: 
            sub1.errorbar(mass_bins.mass_mid, avg_mhalo, yerr=var_mhalo, 
                        lw=4, ls='--', c=pretty_colors[3], label='Integrated Mass')

    sub1.set_xlim([9.0, 11.5]) 
    sub1.set_xlabel(r'$\mathtt{M_*}$', fontsize=24)
    sub1.set_ylabel(r'$\mathtt{M_{Halo}}$', fontsize=24)
    sub1.legend(loc='upper left')

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # Stellar mass specifier 
    mass_str = '_integ_sham_comp'

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 
    
    fig_name1 = ''.join(['figure/', 
        'cenque_mstar_mhalo_snapshot', str(i_nsnap), tau_str, file_type_str, '.png'])
    fig1.savefig(fig_name1, bbox_inches='tight')
    fig1.clear()

def plot_mhalo_mstar_snapshotSHAM(i_nsnap=1, scatter=0.0): 
    ''' Plot M_halo versus SHAM stellar mass for TreePM --> SHAM snapshot

    Parameters
    ----------
    i_nsnap : snapshot # 
    scatter : set scatter of M_halo-M* in SHAM code 

    Notes
    -----
    * Essentially for testing TreePM and SHAM 


    '''
    prettyplot() 
    pretty_colors = prettycolors() 
    
    # import TreePM-->SHAM output 
    if scatter == 0.0: 
        snapshot_file = ''.join(['dat/wetzel_tree/', 
            'subhalo_sham_centrals_snapshot', str(i_nsnap), '.hdf5' 
            ]) 
    else: 
        snapshot_file = ''.join(['dat/wetzel_tree/', 
            'subhalo_sham_centrals_snapshot', str(i_nsnap), '_scatter', str(scatter), '.hdf5' 
            ]) 
    print snapshot_file 
    f = h5py.File(snapshot_file, 'r') # read snapshot file
    grp = f['cenque_data']

    mhalo = grp['halo.m.max'][:]
    mstar = grp['mass'][:]
    
    bovy.scatterplot(mstar, mhalo, scatter=True, color=pretty_colors[1], s=3, 
        xrange = [9.0, 11.5], 
        xlabel=r'$\mathtt{M_{SHAM}}$', ylabel=r'$\mathtt{M_{Halo}}$') 

    fig_name = ''.join(['figure/', 
        'subhalo_sham_centrals_snapshot', str(i_nsnap), '_scatter', str(scatter), '.png'])
    plt.savefig(fig_name, bbox_inches='tight')

def plot_mstar_msham_snapshot(i_nsnap=1, **kwargs): 
    ''' Plot M* versus Msham of galaxies for CenQue file 

    '''
    prettyplot() 
    pretty_colors = prettycolors()  

    # import evolved snapshot 
    if i_nsnap == 13: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='sf assign', **kwargs) 
    else: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
    
    gal_type = snap.gal_type 
    mass = (snap.mass)[gal_type == 'star-forming']
    sham = (snap.sham_mass)[gal_type == 'star-forming']

    bovy.scatterplot(mass, sham, scatter=True, color=pretty_colors[1], 
            xlabel=r'$\mathtt{M_*}$', ylabel=r'$\mathtt{M_{SHAM}}$', 
            xrange=[9.0, 11.5])

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # SFR specifier
    file_type_str = kwargs['sfr']
    
    fig_name = ''.join(['figure/', 
        'cenque_mstar_msham_snapshot', str(i_nsnap), tau_str, file_type_str, '.png'])
    plt.savefig(fig_name, bbox_inches='tight')

def plot_d_stellmass_snapshot(i_nsnap=1, **kwargs): 
    ''' Plot the difference betwen integrated mass and SHAM mass versus SHAM mass for 
    CenQue model  

    '''
    prettyplot() 
    pretty_colors = prettycolors()  

    # import evolved snapshot 
    if i_nsnap == 13: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='sf assign', **kwargs) 
    else: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
    
    gal_type = snap.gal_type 
    mass = (snap.mass)[gal_type == 'star-forming']
    sham = (snap.sham_mass)[gal_type == 'star-forming']
    bovy.scatterplot(sham, sham-mass, scatter=True, color=pretty_colors[1], 
            xlabel=r'$\mathtt{M_{SHAM}}$', ylabel=r'$\mathtt{\Delta M (M_{SHAM}- M_*)}$', 
            xrange=[9.0, 11.5])

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # SFR specifier
    file_type_str = kwargs['sfr']
    
    fig_name = ''.join(['figure/', 
        'cenque_d_mstellmass_snapshot', str(i_nsnap), tau_str, file_type_str, '.png'])
    plt.savefig(fig_name, bbox_inches='tight')

# CenQue SF-MS ----------------------
def plot_cenque_sfms(i_nsnap, **kwargs): 
    ''' Plot SF-MS for CenQue 

    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 

    if i_nsnap < 13: 
        snap = cq.CenQue() 
        snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
    else: 
        snap = cq.CenQue() 
        snap.readin(nsnap=13, file_type='sf assign', **kwargs) 

    sf_index = np.where(snap.gal_type == 'star-forming')    # only keep star forming galaxies
    mass = (snap.mass)[sf_index]
    sfr = (snap.sfr)[sf_index]
    
    # scatter plot with contours!
    bovy.scatterplot(mass, sfr, scatter=True, color=pretty_colors[1], s=3, 
            xlabel=r'$\mathtt{M_*}$', ylabel=r'$\mathtt{SFR}$')

    # figure file name ------------------------------------------------------------
    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # Stellar mass specifier 
    if kwargs['stellmass'].lower() == 'integrated': 
        mass_str = '_integ'
    elif kwargs['stellmass'].lower() == 'sham': 
        mass_str = '_sham'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 
    
    fig_name = ''.join(['figure/cenque_sfms', tau_str, file_type_str, '.png'])
    plt.savefig(fig_name, bbox_inches='tight')
    plt.show() 

if __name__=='__main__': 
    #plot_cenque_quenching_ssfr_dist(10, fqing_yint=-5.84, tau='linear')
    #plot_cenque_quenching_ssfr_dist(1, fqing_yint=-5.84, tau='constant')
    #plot_cenque_quenching_ssfr_dist(1, fqing_yint=-5.84, tau='linefit', tau_param=[-0.15, 0.17])
    #plot_groupcat_obs_fq(Mrcut=18)
    #plot_groupcat_obs_fq(Mrcut=19)
    #plot_groupcat_obs_fq(Mrcut=20)

    #plot_fq_geha_groupcat(Mrcut=18) 
   
    plot_quenching_efold(['linefit'], [[-0.7, 0.4]]) 
    #for i_snap in [12, 6, 1]: 
    #    plot_mhalo_mstar_snapshotSHAM(i_nsnap = i_snap, scatter=0.0)
    #    plot_mhalo_mstar_snapshotSHAM(i_nsnap = i_snap, scatter=0.2)
    
    tau_str = 'linefit'
    tau_param_str = [-0.7, 0.4]
    #sfr_str = 'sfr_func'
    sfr_str = 'sfr_avg'
    stellmass_str = 'sham'
    #stellmass_str = 'integrated'
    #for stellmass_str  in ['sham', 'integrated']: 
    cenque_params = {'tau': tau_str, 'tau_param': tau_param_str, 
            'sfr': sfr_str, 'stellmass': stellmass_str} 
    #plot_cenque_sfms(1, **cenque_params)
    #for i_snap in [12, 6, 1]: 
    #    #plot_mstar_msham_snapshot(i_nsnap=i_snap, **cenque_params)
    #    #plot_mhalo_mstar_sham_integrated(i_nsnap = i_snap, **cenque_params)
    #    plot_mhalo_mstar(i_nsnap=i_snap, **cenque_params)

    #plot_cenque_ssfr_dist_evolution(nsnaps=[2], Mrcut=20, **cenque_params)
    #plot_cenque_ssfr_dist_evolution(nsnaps=[1], Mrcut=19, **cenque_params)
    #plot_cenque_ssfr_dist_evolution(nsnaps=[1], Mrcut=18, **cenque_params)
    '''
    for i in range(1,13): 
        plot_cenque_quenching_ssfr_dist(i, **cenque_params) 

    plot_snapshot_fqobs_evol(nsnaps=[1,2,3,4,5,6,7,8,9,10,11,12], 
            fq_type='wetzelsmooth', **cenque_params)
    for i in range(1,13): 
        plot_snapshot_fqobs(i, fq_type='wetzelsmooth', **cenque_params)
    '''
