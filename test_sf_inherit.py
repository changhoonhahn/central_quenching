import numpy as np
import time

from lineage import Lineage
from sf_inherit import sf_inherit 
from util.cenque_utility import get_t_nsnap

# --- plotting --- 
from plotting.plot_fq import PlotFq
from plotting.plot_sfms import PlotSFMS
from plotting.plot_ssfr import PlotSSFR
from plotting.plot_mstar_mhalo import PlotMstarMhalo

from defutility.plotting import prettyplot
from defutility.plotting import prettycolors

def qaplot_sf_inherit(
        n_step, 
        nsnap_ancestor = 20, 
        nsnap_descendant = 1, 
        scatter = 0.0, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        sfrevol_massdep = False,
        massevol_prop = {'name': 'sham'},
        ssfr=True, 
        mass_scatter=True, 
        sfms=True
        ):
    '''
    '''
    
    if sfrevol_massdep: 
        sfrevol_massdep_str = 'SFRMt_'
    else: 
        sfrevol_massdep_str = 'SFRM0t_'

    file_str = ''.join([ 
        '_', 
        str(nsnap_ancestor), 'ancestor_', 
        str(nsnap_descendant), 'descendant_', 
        str(round(scatter)), 'Mscatter_', 
        sfrevol_prop['name'], '_sfrevol_',
        sfrevol_massdep_str, 
        massevol_prop['name'], '_massevol'
        ])

    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta = np.loadtxt(
            posterior_file, 
            unpack = True
            ) 

    med_theta = [] 
    for i_param in xrange(len(theta)):
        med_theta.append( 
                np.median(theta[i_param])
                )

    if nsnap_descendant < nsnap_ancestor: 
        bloodline = sf_inherit(
                [nsnap_descendant], 
                nsnap_ancestor = nsnap_ancestor, 
                pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
                tau_prop = {
                    'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                    }, 
                sfrevol_prop = sfrevol_prop, 
                sfrevol_massdep = sfrevol_massdep,
                massevol_prop = massevol_prop, 
                quiet = True, 
                scatter = scatter
                )

        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
    else: 
        bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
        bloodline.readin([1], scatter = scatter)
        descendant = bloodline.ancestor_cq

    sfrevol_str = sfrevol_prop['name']
    
    if ssfr: 
        # SSFR distribution
        ssfr_plot = PlotSSFR()
        ssfr_plot.cenque(
                descendant, 
                line_color='red', 
                line_width=4, 
                label=r'Median $\mathtt{\tau_{Q}}$ of Posterior'
                )
        ssfr_plot.groupcat(Mrcut=18)
        ssfr_plot.set_axes()
        ssfr_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_ssfr', file_str, '.png'
            ])
        ssfr_plot.fig.savefig(ssfr_fig_file, bbox_inches="tight")
        
    # Quiescent Fraction 
    #fq_plot = PlotFq()
    #fq_plot.cenque_fq(descendant, line_color='r', label=None)
    #fq_plot.param_fq(line_color='k', label=None)
    #fq_plot.set_axes()
    #fq_fig_file = ''.join([
    #    'figure/', 
    #    'fq_abc_posterior_mass_scatter', str(scatter), '_', sfrevol_sfr, '_sfrevol.png'
    #    ])
    #fq_plot.fig.savefig(fq_fig_file, bbox_inches="tight")
    
    # Quenching Timescale
    #tau_plot = plot_tau([{'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}])
    #tau_fig_file = ''.join([
    #    'figure/', 
    #    'tau_abc_posterior_mass_scatter', str(scatter), '.png'
    #    ])
    #tau_plot.savefig(tau_fig_file, bbox_inches="tight")
    #plt.close()
    
    if sfms:
        # Star Forming Main Sequence 
        sfms_plot = PlotSFMS()
        sfms_plot.cenque(descendant, justsf=True)
        sfms_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_sfms', file_str, '.png'
            ])
        plt.savefig(sfms_fig_file, bbox_inches="tight")
        plt.close()
   
    if mass_scatter: 
        # Stellar Mass - Halo Mass 
        mass_scatter_plot = PlotMstarMhalo()
        mass_scatter_plot.cenque(descendant)
            
        mass_scatter_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_mass_scatter', file_str, '.png'
            ])
        plt.savefig(mass_scatter_fig_file, bbox_inches="tight")
        plt.close()

def qaplot_sf_inherit_nosfr_scatter(
        n_step, 
        nsnap_descendants,
        nsnap_ancestor = 20, 
        scatter = 0.0, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        sfrevol_massdep = False,
        massevol_prop = {'name': 'sham'}
        ):
    '''
    Test Lineage SF Inherit function for the case where there's no scatter in 
    the SFMS
    '''
    
    if sfrevol_massdep: 
        sfrevol_massdep_str = 'SFRMt_'
    else: 
        sfrevol_massdep_str = 'SFRM0t_'

    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta = np.loadtxt(
            posterior_file, 
            unpack = True
            ) 

    med_theta = [] 
    for i_param in xrange(len(theta)):
        med_theta.append( 
                np.median(theta[i_param])
                )

    bloodline = sf_inherit(
            nsnap_descendants, 
            nsnap_ancestor = nsnap_ancestor, 
            ancestor_sf_prop = {'name': 'average_noscatter'}, 
            pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
            tau_prop = {
                'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                }, 
            sfrevol_prop = sfrevol_prop, 
            sfrevol_massdep = sfrevol_massdep,
            massevol_prop = massevol_prop, 
            quiet = True, 
            scatter = scatter
            )
    
    prettyplot()
    pretty_colors = prettycolors()
    
    sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

    for nsnap_descendant in nsnap_descendants: 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))

        fig = plt.figure()
        sub = fig.add_subplot(111)

        sf_gal = np.where(descendant.gal_type == 'star-forming')

        sub.scatter(
                descendant.mass[sf_gal], 
                descendant.sfr[sf_gal], 
                color=pretty_colors[nsnap_descendant],
                s=10
                )
        sub.plot(np.arange(8.5, 12.0, 0.1), 
                sfr_mstar_z(np.arange(8.5, 12.0, 0.1), get_z_nsnap(nsnap_descendant)), 
                c = 'k', 
                ls = '--', 
                lw = 3 
                )
        sub.set_xlim([9.0, 12.0])
        sub.set_ylim([-5.0, 2.0])
        sub.set_xlabel(r'$\mathtt{M_*}$')
        sub.set_ylabel(r'$\mathtt{log\;SFR}$')
    
        file_str = ''.join([ 
            '_', 
            str(nsnap_ancestor), 'ancestor_', 
            str(nsnap_descendant), 'descendant_', 
            str(round(scatter)), 'Mscatter_', 
            sfrevol_prop['name'], '_sfrevol_',
            sfrevol_massdep_str, 
            massevol_prop['name'], '_massevol'
            ])

        sfms_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_sfms', file_str, '_no_sfr_scatter_twoslope_sfms.png'
            ])
        fig.savefig(sfms_fig_file, bbox_inches='tight')
        plt.close()

def qaplot_sf_inherit_average_scatter(
        n_step, 
        nsnap_descendants,
        nsnap_ancestor = 20, 
        scatter = 0.0, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        sfrevol_massdep = False,
        massevol_prop = {'name': 'sham'}
        ):
    '''
    Test Lineage SF Inherit function
    '''

    # Get media theta from ABC posterior
    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta = np.loadtxt(
            posterior_file, 
            unpack = True
            ) 

    med_theta = [] 
    for i_param in xrange(len(theta)):
        med_theta.append( 
                np.median(theta[i_param])
                )
    # ------------------
    # SF inheritance 
    start_time = time.time()
    bloodline = sf_inherit(
            nsnap_descendants, 
            nsnap_ancestor = nsnap_ancestor, 
            ancestor_sf_prop = {'name': 'average'}, 
            pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
            tau_prop = {
                'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                }, 
            sfrevol_prop = sfrevol_prop, 
            sfrevol_massdep = sfrevol_massdep,
            massevol_prop = massevol_prop, 
            quiet = True, 
            scatter = scatter
            )
    print 'SF inherit took ', time.time() - start_time
    
    if sfrevol_massdep: 
        sfrevol_massdep_str = 'SFRMt_'
    else: 
        sfrevol_massdep_str = 'SFRM0t_'
    
    for nsnap_descendant in nsnap_descendants: 

        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
        
        start_time = time.time()
        # SFMS plot
        sfms_plot = descendant.plotSFMS(justsf=True, bovyplot=True)
        sfms_plot.param_sfms()      # plot parameterized SFMS

        file_str = ''.join([ 
            '_', 
            str(nsnap_ancestor), 'ancestor_', 
            str(nsnap_descendant), 'descendant_', 
            str(round(scatter)), 'Mscatter_', 
            sfrevol_prop['name'], '_sfrevol_',
            sfrevol_massdep_str, 
            massevol_prop['name'], '_massevol'
            ])

        sfms_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_sfms_justsf', file_str, '_twoslope_sfms.png'
            ])
        print sfms_fig_file 
        sfms_plot.save_fig(sfms_fig_file)
        print 'SFMS plotting took ', time.time() - start_time

        plt.close()

        sfms_plot = descendant.plotSFMS(justsf=False, bovyplot=True)
        sfms_plot.param_sfms()      # plot parameterized SFMS

        sfms_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_sfms', file_str, '_twoslope_sfms.png'
            ])
        print sfms_fig_file 
        sfms_plot.save_fig(sfms_fig_file)
        print 'SFMS plotting took ', time.time() - start_time

        plt.close()

def sf_inherited_lineage(
        n_step, 
        nsnap_ancestor = 20, 
        scatter = 0.0, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        sfrevol_massdep = False,
        massevol_prop = {'name': 'sham'}
        ): 
    '''
    Write out SF Inherited Lineage to file in order for faster access in the future. 
    '''
    if sfrevol_massdep: 
        sfrevol_massdep_str = 'SFRMt_'
    else: 
        sfrevol_massdep_str = 'SFRM0t_'

    file_str = ''.join([ 
        '_', 
        str(nsnap_ancestor), 'ancestor_', 
        str(round(scatter)), 'Mscatter_', 
        sfrevol_prop['name'], '_sfrevol_',
        sfrevol_massdep_str, 
        massevol_prop['name'], '_massevol_',
        str(n_step), '_abcpost'
        ])

    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta = np.loadtxt(
            posterior_file, 
            unpack = True
            ) 

    med_theta = [] 
    for i_param in xrange(len(theta)):
        med_theta.append( 
                np.median(theta[i_param])
                )

    nsnap_range = np.arange(1, nsnap_ancestor)

    bloodline = sf_inherit(
            list(nsnap_range),
            nsnap_ancestor = nsnap_ancestor, 
            pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
            tau_prop = {
                'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                }, 
            sfrevol_prop = sfrevol_prop, 
            sfrevol_massdep = sfrevol_massdep,
            massevol_prop = massevol_prop, 
            quiet = True, 
            scatter = scatter
            )

    bloodline.file_name = ''.join([
            'dat/lineage/', 
            'lineage_ancestor_', 
            str(bloodline.nsnap_ancestor), 
            '_descendants',
            '_mass_scatter', 
            str(round(bloodline.mass_scatter, 1)), 
            file_str, 
            '.hdf5'
            ]) 
    bloodline.writeout()

    return None

def track_sf_pop_sfms_evol(
        n_step, 
        n_gal, 
        nsnap_ancestor = 20, 
        scatter = 0.0, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        sfrevol_massdep = False,
        massevol_prop = {'name': 'sham'}
        ):
    '''
    Directly track the evolution of n_gal galaxies in SFR versus M* paradigm. 
    Makes cool plots, 10 out of 10 would do again.
    '''

    if sfrevol_massdep: 
        sfrevol_massdep_str = 'SFRMt_'
    else: 
        sfrevol_massdep_str = 'SFRM0t_'

    file_str = ''.join([ 
        '_', 
        str(nsnap_ancestor), 'ancestor_', 
        str(round(scatter)), 'Mscatter_', 
        sfrevol_prop['name'], '_sfrevol_',
        sfrevol_massdep_str, 
        massevol_prop['name'], '_massevol_',
        str(n_step), '_abcpost'
        ])
    
    nsnap_range = np.arange(1, nsnap_ancestor)
    
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.file_name = ''.join([
            'dat/lineage/', 
            'lineage_ancestor_', 
            str(nsnap_ancestor), 
            '_descendants',
            '_mass_scatter', 
            str(round(scatter, 1)), 
            file_str, 
            '.hdf5'
            ]) 
    bloodline.readin(nsnap_range)

    for i_snap in nsnap_range: 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(i_snap))
        
        if i_snap == nsnap_range[0]: 
            gal_assigned = np.where(descendant.gal_type != 'quiescent')

        if i_snap == nsnap_range[0]: 
            sf_masses = descendant.mass[gal_assigned]
            sf_sfrs = descendant.sfr[gal_assigned]
        else: 
            sf_masses = np.vstack((sf_masses, descendant.mass[gal_assigned]))
            sf_sfrs = np.vstack((sf_sfrs, descendant.sfr[gal_assigned]))
    
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    for i in xrange(n_gal): #sf_masses.shape[1]):
        sub.plot(
                sf_masses[:,i], 
                sf_sfrs[:,i],
                color=pretty_colors[i % 20],
                lw=2
                )
        sub.scatter(
                sf_masses[:,i], 
                sf_sfrs[:,i],
                color=pretty_colors[i % 20],
                s=16
                )
    sub.set_xlim([9.0, 12.0])
    sub.set_ylim([-5.0, 2.0])
    sub.set_xlabel(r'$\mathtt{M_*}$')
    sub.set_ylabel(r'$\mathtt{log\;SFR}$')
    
    sfms_fig_file = ''.join([
        'figure/', 
        'qaplot_sf_inherit_sfms', file_str, '_sf_galaxytrack.png'
        ])
    fig.savefig(sfms_fig_file, bbox_inches="tight")
    plt.close()
    
    for attr in ['sfr', 'mass', 'ssfr']: 

        fig = plt.figure(figsize=(15,8))
        sub = fig.add_subplot(111)
    
        for i in xrange(n_gal): #sf_masses.shape[1]):
            
            if attr == 'sfr': 
                attr_val = sf_sfrs[:,i]
            elif attr == 'mass': 
                attr_val = sf_masses[:,i]
            elif attr == 'ssfr': 
                attr_val = sf_sfrs[:,i] - sf_masses[:,i]

            sub.plot( get_t_nsnap(nsnap_range), 
                    attr_val, 
                    color=pretty_colors[i % 20],
                    lw=2
                    )

        sub.set_xlim([3, 14.0])
        sub.set_ylabel(r'$\mathtt{log\;'+attr.upper()+'}$')
        sub.set_xlabel(r'$\mathtt{t_{cosmic}}$')
    
        attr_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_', attr, '_', file_str, '_galaxytrack.png'
            ])
        fig.savefig(attr_fig_file, bbox_inches="tight")
        #plt.show()
        plt.close()

if __name__=="__main__":
    #qaplot_sf_inherit_average_scatter(
    #        29, 
    #        [1],
    #        sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
    #        sfrevol_massdep = True, 
    #        massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.1}
    #        )

    sf_inherited_lineage(
            29, 
            nsnap_ancestor = 20, 
            scatter = 0.0, 
            sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
            sfrevol_massdep = True, 
            massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.01}
            )

    #track_sf_pop_sfms_evol(
    #        29, 
    #        5,
    #        sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
    #        sfrevol_massdep = True, 
    #        massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.01}
    #        )
    
'''
    qaplot_sf_inherit_average_scatter(
            29, 
            [1], #[1,3,5,7,9,11,13,15,17,19],
            sfrevol_prop = {'name': 'amp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1]},
            sfrevol_massdep = True, 
            massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.1}
            )

    for sfrevoprop in [{'name': 'notperiodic'}, {'name': 'squarewave', 'freq_range': [100.0, 1000.0], 'phase_range': [0, 1]}]:
        track_sf_pop_sfms_evol(
                29, 
                10,
                sfrevol_prop = sfrevoprop,
                sfrevol_massdep = False, 
                massevol_prop = {'name': 'sham'}
                )
        track_sf_pop_sfms_evol(
                29, 
                10, 
                sfrevol_prop = sfrevoprop, 
                sfrevol_massdep = False, 
                massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.05}
                )
        track_sf_pop_sfms_evol(
                29, 
                10,
                sfrevol_prop = sfrevoprop,
                sfrevol_massdep = True, 
                massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.05}
                )

    for i in np.arange(1,20)[::-1]: #[1, 5]:#, 13, 18, 19]:#, 11, 15, 19]: 
        for sfrevoprop in [{'name': 'notperiodic'}, {'name': 'squarewave', 'freq_range': [100.0, 1000.0], 'phase_range': [0, 1]}]:
            print i 
            qaplot_sf_inherit(
                    29, 
                    nsnap_descendant = i, 
                    sfrevol_prop = sfrevoprop,
                    sfrevol_massdep = False, 
                    massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.01}, 
                    ssfr=False, 
                    sfms=True, 
                    mass_scatter=False
                    )
            
            qaplot_sf_inherit(
                    29, 
                    nsnap_descendant = i, 
                    sfrevol_prop = sfrevoprop,
                    sfrevol_massdep = True, 
                    massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.01}, 
                    ssfr=False, 
                    sfms=True, 
                    mass_scatter=False
                    )


        #qaplot_sf_inherit(
        #        29, 
        #        nsnap_descendant = i, 
        #        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        #        massevol_prop = {'name': 'integrated', 'f_retain': 1.0, 't_step': 0.1}
        #        )

        #qaplot_sf_inherit(
        #        29, 
        #        nsnap_descendant = i, 
        #        sfrevol_prop = {'name': 'notperiodic'}
        #        )
'''
