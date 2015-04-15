'''

Plotting code for CentralQuenching (CenQue) project

Author(s): ChangHoon Hahn 

'''

import numpy as np
import matplotlib.pyplot as plt

#----- Local -----
import cenque_utility as util 
import cenque as cq 
import cenque_groupcat as cq_group

def plot_cenque_ssfr_dist(cenque, fig=None): 
    ''' Plot CenQue data
    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
    
    if fig == None: 
        fig = plt.figure(1, figsize=(15,15))
        fig.subplots_adjust(hspace=0., wspace=0.)
        new_plot = True
    else: 
        new_plot = False

    # mass bins of the panel        ( hardcoded ) 
    panel_mass_bins = [
            [9.5, 10.0], [10.0, 10.5], 
            [10.5, 11.0], [11.0, 11.5]
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

        try: 
            cenque.zsnap
        except AttributeError: 
            ssfr_hist_label = 'Centrals'     # redshift lable 
        else: 
            ssfr_hist_label = r'$z='+str(cenque.zsnap)+'$'     # redshift lable 

        try: 
            line_color = pretty_colors[cenque.nsnap]
        except TypeError: 
            line_color = 'black'

        sub = fig.add_subplot(2, 2, i_mass+1)       # panel subplot 
        sub.plot(ssfr_bin_mid, ssfr_hist, 
                color=line_color, lw=4, label=ssfr_hist_label) 
        
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
            sub.set_xticklabels([])
        elif i_mass == 1:
            sub.set_yticklabels([])
            sub.set_xticklabels([])
        elif i_mass == 2: 
            sub.set_ylabel(r'$P(log \; SSFR)$') 
            sub.set_xlabel(r'$log \; SSFR \;[yr^{-1}]$') 
        else: 
            sub.set_yticklabels([])
            sub.set_xlabel(r'$log \; SSFR \;[yr^{-1}]$') 
                
            try: 
                fig_leg.remove() 
            except UnboundLocalError: 
                pass
            
            fig_leg = sub.legend(loc='lower right', prop={'size':'24'})        # legends 

    return fig   

def plot_cenque_sf_mainseq(cenque): 
    ''' Given CenQue data, plot the Star-Forming Main Sequence
    '''
    pass

def plot_cenque_ssfr_dist_evolution(**kwargs): 
    snap = cq.CenQue() 
    snap.readin(nsnap=13, file_type='sf assign', **kwargs) 
    ssfr_fig = plot_cenque_ssfr_dist(snap)

    for i_nsnap in reversed(range(1,13)):
        next_snap = cq.CenQue() 
        next_snap.readin(nsnap=i_nsnap, file_type='evol from 13', **kwargs) 
        
        ssfr_fig = plot_cenque_ssfr_dist(next_snap, fig=ssfr_fig) 

    central_ssfr = cq_group.central_catalog(Mrcut=18) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig) 
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'cenque_ssfr_evol_', kwargs['tau'], 'tau_', kwargs['fq'], 'fq.png']) 
    ssfr_fig.savefig(fig_file, bbox_inches='tight') 
    ssfr_fig.clear() 

def plot_fq_evol(): 
    ''' plot fq evolution 
    '''

    zbin = [0.1*np.float(i) for i in range(1,10)]   # zbin 
    ###
    ###     Get z_snapshots 
    ###
    ###
    ###
    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    fig, subs = plt.subplots(1,3, figsize=[15, 5]) 
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

if __name__=='__main__': 
    plot_cenque_ssfr_dist_evolution(fq='cosmosinterp', tau='instant') 
    plot_cenque_ssfr_dist_evolution(fq='cosmosinterp', tau='constant') 
    plot_cenque_ssfr_dist_evolution(fq='cosmosinterp', tau='linear') 

    #fq_fig = plot_fq_evol() 
    #fq_fig.savefig('fq_evol_fig.png', bbox_inches='tight')
    #
    #tau_fig = plot_quenching_efold() 
    #tau_fig.savefig('quenching_efold_fig.png', bbox_inches='tight')

