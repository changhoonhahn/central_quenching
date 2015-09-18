'''

Plot CenQue class object


'''
import numpy as np
import matplotlib.pyplot as plt

# --- Local --- 
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

def plot_cenque_ssfr_dist(cenque, fig=None, **kwargs): 
    ''' Plot sSFR distribution for CenQue data
    '''
    #prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
    
    if fig == None: 
        fig = plt.figure(1, figsize=(16,16))
        fig.subplots_adjust(hspace=0., wspace=0.)
        new_plot = True
    else: 
        new_plot = False

    # mass bins of the panel        ( hardcoded ) 
    panel_mass_bins = [ 
            [9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]
            ]
    
    for i_mass, panel_mass in enumerate(panel_mass_bins):       # loop through each panel 

        mass_limit = np.where(
                (cenque.mass >= panel_mass[0]) & 
                (cenque.mass < panel_mass[1])
                )
        ngal_bin = len(mass_limit[0])

        # SSFR histogram 
        ssfr_hist, ssfr_bin_edges = np.histogram(
                cenque.ssfr[mass_limit], 
                range = [-13.0, -7], 
                bins = 40, 
                normed = True
                )
        ssfr_bin_low = ssfr_bin_edges[:-1]
        ssfr_bin_high = ssfr_bin_edges[1:]
        ssfr_bin_mid = 0.5 * (ssfr_bin_low + ssfr_bin_high) 
    
        if 'label' in kwargs: 
            ssfr_hist_label = kwargs['label']
        else: 
            try: 
                ssfr_hist_label = r'$\mathtt{z='+str(cenque.zsnap)+'}$' 
            except AttributeError: 
                ssfr_hist_label = 'Centrals'
    
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

        sub = fig.add_subplot(2, 2, i_mass+1)       # panel subplot 
        sub.plot(
                ssfr_bin_mid, 
                ssfr_hist, 
                color = line_color, 
                lw = line_width, 
                ls = line_style, 
                label = ssfr_hist_label) 

        if new_plot:        # put mass labels for the panels 
            plt.text(-11.25, 1.4, 
                    r'$\mathtt{log \; M_{*} = ['+str(panel_mass[0])+',\;'+str(panel_mass[1])+']}$', 
                    fontsize=24
                    )

        sub.set_xlim([-13.0, -7.0])
        sub.set_ylim([0.0, 1.6])

        if i_mass == 0: 
            sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=20) 
            sub.set_xticklabels([])
        elif i_mass == 1: 
            sub.set_xticklabels([])
            sub.set_yticklabels([])
        elif i_mass == 2:
            #sub.set_yticklabels([])
            sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=20) 
            sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=20) 
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=20) 
                
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
    snap = cq.CenQue(n_snap=13, cenque_type = 'sf_assgined') 
    snap.readin()  
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
