import numpy as np
import time

from lineage import Lineage
from sf_inherit import sf_inherit 
from util.cenque_utility import get_t_nsnap

# --- plotting --- 
from plotting.plots import PlotFq
from plotting.plots import PlotSFMS
from plotting.plots import PlotSSFR
from plotting.plots import PlotMstarMhalo

from defutility.plotting import prettyplot
from defutility.plotting import prettycolors

def qaplot_sf_inherit(
        nsnap_ancestor = 20, nsnap_descendant = 1, n_step = 29, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
        evol_prop = {
            'sfr': {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
            'mass': {'name': 'sham'}},
        ssfr=True, fq=False, tau=False, mass_scatter=False, sfms=False, smf=False
        ):
    '''
    QAPlots for SF inherit module. 
    '''
    med_theta = abc_posterior_median(n_step)
    # quenching probabilyt and tau properties from ABC posterior
    pq_dict = {'slope': med_theta[0], 'yint': med_theta[1]}
    tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}
    evol_prop['pq'] = pq_dict
    evol_prop['tau'] = tau_dict

    bloodline = sf_inherit(
            nsnap_descendant, 
            nsnap_ancestor = nsnap_ancestor, 
            subhalo_prop = subhalo_prop, 
            sfr_prop = sfr_prop,
            evol_prop = evol_prop)

    descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))

    if evol_prop['mass']['name'] == 'sham': 
        massevol_str = 'sham'
    else: 
        massevol_str = ''.join([evol_prop['mass']['type'], evol_prop['mass']['name']])

    lineage_str = ''.join([ 
        '_', 
        str(nsnap_ancestor), 'ancestor_', 
        str(nsnap_descendant), 'descendant_', 
        str(round(subhalo_prop['scatter'], 1)), 'Mscatter_', 
        evol_prop['sfr']['name'], '_sfrevol_SFRMt_',
        massevol_str, '_massevol'
        ])

    attr_list = []
    if ssfr: 
        attr_list.append('ssfr')
    if fq: 
        attr_list.append('fq')
    if tau: 
        attr_list.append('tau')
    if sfms: 
        attr_list.append('sfms')
    if mass_scatter: 
        attr_list.append('mass_scatter')
    if smf: 
        attr_list.append('smf')

    for attr in attr_list: 
        fig_name = ''.join([
            'figure/', 
            'qaplot_sf_inherit_', attr, lineage_str, '.png'
            ])

        if attr == 'ssfr':                          # SSFR
            descendant.plotSsfr(
                    line_color='red', 
                    line_width=4, 
                    label=r'Median $\mathtt{\tau_{Q}}$ of Posterior', 
                    groupcat = True, 
                    savefig = fig_name                    
                    )
            # quenching = True,
        elif attr == 'fq':                          # Quiescent Fraction 
            descendant.plotFq(
                    param=True, 
                    line_color='r', 
                    label = None,
                    savefig= fig_name                    
                    )
        elif attr == 'tau':                         # Quenching Timescale
            descendant.plotTau(tau_dict, savefig = fig_name)
        elif attr == 'smf':                         # Stellar Mass Function
            descendant.plotSMF(savefig = fig_name)
        elif attr == 'sfms':                        # Star Forming Main Sequence 
            descendant.plotSFMS(justsf=True, bovyplot=True, savefig = fig_name)
        elif attr == 'mass_scatter':                # Stellar Mass - Halo Mass 
            descendant.plotMstarMhalo(savefig = fig_name)
    
    return None

def qaplot_sf_inherit_average_scatter(
        nsnap_descendants, nsnap_ancestor = 20, scatter = 0.0, n_step = 29,
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        massevol_prop = {'name': 'sham'}
        ):
    '''
    Test Lineage SF Inherit function
    '''
    med_theta = abc_posterior_median(n_step)         # Get media theta from ABC posterior

    # quenching probabilyt properties
    pq_dict = {'slope': med_theta[0], 'yint': med_theta[1]}
    # tau properties
    tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}

    start_time = time.time()
    bloodline = sf_inherit(
            nsnap_descendants, 
            nsnap_ancestor = nsnap_ancestor, 
            ancestor_sf_prop = {'name': 'average'}, 
            pq_prop = pq_dict, 
            tau_prop = tau_dict,            
            sfrevol_prop = sfrevol_prop, 
            massevol_prop = massevol_prop, 
            quiet = True, 
            scatter = scatter
            )
    print 'SF inherit took ', (time.time() - start_time)/60.0, ' minutes' 
    
    for nsnap_descendant in nsnap_descendants: 

        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
        
        if massevol_prop['name'] == 'sham': 
            massevol_str = massevol_prop['name']
        else: 
            massevol_str = massevol_prop['type'] + massevol_prop['name']

        lineage_str = ''.join([ 
            '_', 
            str(nsnap_ancestor), 'ancestor_', 
            str(nsnap_descendant), 'descendant_', 
            str(round(scatter)), 'Mscatter_', 
            sfrevol_prop['name'], '_sfrevol_SFRMt_',
            massevol_str, '_massevol'
            ])

        for justsf in [True, False]: 
            justsf_str = ''
            if justsf: 
                justsf_str = '_justsf'

            fig_name = ''.join([
                'figure/', 
                'qaplot_sf_inherit_sfms', lineage_str, justsf_str, '_twoslope_sfms.png'
                ])

            sfms_plot = descendant.plotSFMS(justsf=justsf, bovyplot=True)
            sfms_plot.param_sfms()
            
            print fig_name
            sfms_plot.save_fig(fig_name)

        plt.close()

def sf_inherited_lineage(
        n_step, 
        nsnap_ancestor = 20, 
        scatter = 0.0, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        massevol_prop = {'name': 'sham'}
        ): 
    '''
    Write out SF Inherited Lineage to file for record keeping purposes and faster access for plotting 
    '''

    med_theta = abc_posterior_median(n_step)

    nsnap_range = np.arange(1, nsnap_ancestor)
    
    bloodline = sf_inherit(
            list(nsnap_range),
            nsnap_ancestor = nsnap_ancestor, 
            pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
            tau_prop = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}, 
            sfrevol_prop = sfrevol_prop, 
            massevol_prop = massevol_prop, 
            quiet = True, 
            scatter = scatter
            )

    if massevol_prop['name'] == 'sham': 
        massevol_str = ''.join([massevol_prop['name'], '_Mevol_'])
    else: 
        massevol_str = ''.join([
            massevol_prop['type'], massevol_prop['name'], '_Mevol_', 
            str(round(massevol_prop['t_step'], 3)), '_tstep_'
            ])

    lineage_str = ''.join([ 
        '_', 
        str(round(scatter)), 'Mscatter_', 
        sfrevol_prop['name'], '_sfrevol_SFRMt_', 
        massevol_str, str(n_step), '_abcpost'
        ])

    bloodline.file_name = ''.join([
            'dat/lineage/', 
            'lineage_ancestor_', 
            str(bloodline.nsnap_ancestor), 
            '_descendants',
            lineage_str, 
            '.hdf5'
            ]) 
    bloodline.writeout()

    return None

def track_sf_pop_sfms_evol(
        n_gal, nsnap_ancestor = 20, n_step=29, scatter = 0.0, justsf=True, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        massevol_prop = {'name': 'sham'}
        ):
    '''
    Directly track the evolution of n_gal galaxies in SFR versus M* paradigm. 
    Makes cool plots, 10 out of 10 would do again.
    '''
    if massevol_prop['name'] == 'sham': 
        massevol_str = ''.join([massevol_prop['name'], '_Mevol_'])
    else: 
        massevol_str = ''.join([
            massevol_prop['type'], massevol_prop['name'], '_Mevol_', 
            str(round(massevol_prop['t_step'], 3)), '_tstep_'
            ])

    lineage_str = ''.join([ 
        '_', 
        str(round(scatter)), 'Mscatter_', 
        sfrevol_prop['name'], '_sfrevol_SFRMt_', 
        massevol_str, str(n_step), '_abcpost'
        ])

    nsnap_range = np.arange(1, nsnap_ancestor)
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.file_name = ''.join([
            'dat/lineage/', 
            'lineage_ancestor_', 
            str(bloodline.nsnap_ancestor), 
            '_descendants',
            lineage_str, 
            '.hdf5'
            ]) 
    bloodline.readin(nsnap_range)

    for i_snap in nsnap_range: 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(i_snap))
        
        if i_snap == nsnap_range[0]: 
            if justsf: 
                gal_assigned = np.where(descendant.gal_type != 'quiescent')
            else: 
                gal_assigned = np.where(descendant.gal_type != '')

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
    
    justsf_str = ''
    if justsf: 
        justsf_str = '_justsf'
    sfms_fig_file = ''.join([
        'figure/', 
        'qaplot_sf_inherit_sfms', lineage_str, justsf_str, '_galaxytrack.png'
        ])
    fig.savefig(sfms_fig_file, bbox_inches="tight")
    plt.close()
    
    fig = plt.figure(figsize=(15,15))
    fig.subplots_adjust(hspace=0., wspace=0.)
    
    for i_attr, attr in enumerate(['sfr', 'mass', 'ssfr']): 
        sub = fig.add_subplot(3, 1, i_attr+1) 

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

        if i_attr == 2:
            sub.set_ylim([-11.0, -8.0])
            sub.set_xlabel(r'$\mathtt{t_{cosmic}}$', fontsize=30)
        elif i_attr == 1: 
            sub.set_ylim([7.5, 11.0])
            sub.set_xticklabels([])
        elif i_attr == 0: 
            sub.set_ylim([-2.0, 1.5])
            sub.set_xticklabels([])
    
    attr_fig_file = ''.join([
        'figure/', 
        'qaplot_sf_inherit_ssfr', lineage_str, justsf_str, '_galaxytrack.png'
        ])
    fig.savefig(attr_fig_file, bbox_inches="tight")
    #plt.show()
    plt.close()

def qaplot_sf_inherited_lineage(
        nsnap_range=None, nsnap_ancestor = 20, scatter = 0.0, n_step = 29, 
        sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
        massevol_prop = {'name': 'sham'}, 
        ssfr=True, fq=False, tau=False, sfms=False, mass_scatter=False, smf=False
        ): 
    '''
    '''
    if massevol_prop['name'] == 'sham': 
        massevol_str = ''.join([massevol_prop['name'], '_Mevol_'])
    else: 
        massevol_str = ''.join([
            massevol_prop['type'], massevol_prop['name'], '_Mevol_', 
            str(round(massevol_prop['t_step'], 3)), '_tstep_'
            ])

    lineage_str = ''.join([ 
        '_', 
        str(round(scatter)), 'Mscatter_', 
        sfrevol_prop['name'], '_sfrevol_SFRMt_', 
        massevol_str, str(n_step), '_abcpost'
        ])
    
    if nsnap_range is None: 
        nsnap_range = np.arange(1, nsnap_ancestor)

    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.file_name = ''.join([
            'dat/lineage/', 
            'lineage_ancestor_', 
            str(bloodline.nsnap_ancestor), 
            '_descendants',
            lineage_str, 
            '.hdf5'
            ]) 
    bloodline.readin(nsnap_range)
    
    attr_list = []
    if ssfr: 
        attr_list.append('ssfr')
    if fq: 
        attr_list.append('fq')
    if tau: 
        attr_list.append('tau')
    if sfms: 
        attr_list.append('sfms')
    if smf: 
        attr_list.append('smf')
    if mass_scatter: 
        attr_list.append('mass_scatter')
    
    for i_snap in nsnap_range: 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(i_snap))

        for attr in attr_list: 
            fig_name = ''.join([
                'figure/', 
                'qaplot_sf_inherit_lineage_', attr, lineage_str, '_', str(i_snap), 'of', str(nsnap_ancestor),'.png'
                ])

            if attr == 'ssfr':                          # SSFR
                descendant.plotSsfr(
                        line_color='red', 
                        line_width=4, 
                        label=r'Median $\mathtt{\tau_{Q}}$ of Posterior', 
                        groupcat = True, 
                        savefig = fig_name                    
                        )
            elif attr == 'fq':                          # Quiescent Fraction 
                descendant.plotFq(
                        param=True, 
                        line_color='r', 
                        label = None,
                        savefig= fig_name                    
                        )
            elif attr == 'tau':                         # Quenching Timescale
                descendant.plotTau(savefig = fig_name)
            elif attr == 'sfms':                        # Star Forming Main Sequence 
                descendant.plotSFMS(justsf=True, bovyplot=True, savefig = fig_name)
            elif attr == 'smf': 
                descendant.plotSMF(savefig = fig_name)
            elif attr == 'mass_scatter':                # Stellar Mass - Halo Mass 
                descendant.plotMstarMhalo(savefig = fig_name)
    return None

def abc_posterior_median(n_step): 
    '''
    Median theat from ABC posterior 
    '''
    # read in ABC posterior file 
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

    return med_theta 

"""
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
"""

if __name__=="__main__":

    for scat in [0.0, 0.2]: 
        start_time = time.time()
        bloodline = Lineage(nsnap_ancestor = 20)
        bloodline.descend(subhalo_prop = {'scatter': scat, 'source': 'li-march'}, clobber=True) 
        bloodline.assign_sfr_ancestor(sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}})
        bloodline.writeout()
        print 'lineage construction and write out takes ', (time.time() - start_time)/60.0

        for id in [1, 5, 9, 13, 17]:#, 3, 5, 7, 9, 11, 13, 15, 17, 19]:
            qaplot_sf_inherit(
                nsnap_ancestor = 20, nsnap_descendant = id, 
                subhalo_prop = {'scatter': scat, 'source': 'li-march'}, 
                sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
                evol_prop = {
                    'sfr': {'name': 'newamp_squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
                    'mass': {'name': 'sham'} 
                    },
                ssfr=True, fq=True, tau=False, mass_scatter=True, sfms=True, smf=True
                )
            qaplot_sf_inherit(
                nsnap_ancestor = 20, nsnap_descendant = id, 
                subhalo_prop = {'scatter': scat, 'source': 'li-march'}, 
                sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
                evol_prop = {
                    'sfr': {'name': 'newamp_squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
                    'mass': {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.05} 
                    },
                ssfr=True, fq=True, tau=False, mass_scatter=True, sfms=True, smf=True
                )

    # {'name': 'squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1]}
    # {'name': 'newamp_squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},

    # {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.05}
    #'name': 'squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1]
    #qaplot_sf_inherit_average_scatter(
    #        [1],
    #        sfrevol_prop = 
    #        massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.01}
    #        )   # {'name': 'notperiodic'}

    #for integ in ['euler', 'rk4']:
    #start_time = time.time()
    #sf_inherited_lineage(
    #        29, 
    #        nsnap_ancestor = 20, 
    #        scatter = 0.0, 
    #        sfrevol_prop = {
    #            'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3
    #            },
    #        massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.0025}
    #        )
    #print (time.time() - start_time)/60.0, ' minutes'
    #for tstep in [0.01, 0.0025]: 
    #    for integ in ['euler', 'rk4']: 
    #        qaplot_sf_inherited_lineage(
    #                nsnap_range=[1],
    #                sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3}, 
    #                massevol_prop = {'name': 'integrated', 'type': integ, 'f_retain': 0.6, 't_step': tstep}, 
    #                ssfr=False, fq=False, tau=False, sfms=False, mass_scatter=False, smf=True
    #                )
    #track_sf_pop_sfms_evol(
    #        10,
    #        sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
    #        massevol_prop = {'name': 'integrated', 'type': 'rk4', 'f_retain': 0.6, 't_step': 0.0025}, 
    #        justsf=False
    #        )
