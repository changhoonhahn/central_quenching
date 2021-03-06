import numpy as np

from sf_inherit import sf_inherit 
from plotting.plot_fq import PlotFq
from plotting.plot_tau import plot_tau
from plotting.plot_sfms import PlotSFMS
from plotting.plot_ssfr import PlotSSFR
from plotting.plot_mstar_mhalo import PlotMstarMhalo


def qaplot_sf_inherit(n_step, sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]}):
    '''
    '''
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
            [1], 
            nsnap_ancestor = nsnap_ancestor, 
            pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
            tau_prop = {
                'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                }, 
            sfrevol_prop = sfrevol_prop, 
            quiet = True, 
            scatter = scatter
            )

    sfrevol_str = sfrevol_prop['name']

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
        'ssfr_abc_posterior_mass_scatter', str(scatter), '_', sfrevol_str, '_sfrevol.png'
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
    
    # Star Forming Main Sequence 
    sfms_plot = PlotSFMS()
    sfms_plot.cenque(descendant)
    sfms_fig_file = ''.join([
        'figure/', 
        'sfms_abc_posterior_mass_scatter', str(scatter), '_', sfrevol_str, '_sfrevol.png'
        ])
    plt.savefig(sfms_fig_file, bbox_inches="tight")
    plt.close()
    
    # Stellar Mass - Halo Mass 
    #mass_scatter_plot = PlotMstarMhalo()
    #mass_scatter_plot.bloodline(
    #        bloodline, 
    #        nsnap_descendant
    #        )
    #    
    #mass_scatter_fig_file = ''.join([
    #    'figure/', 
    #    'mass_scatter_abc_posterior_mass_scatter', str(scatter), '.png'
    #    ])
    #plt.savefig(mass_scatter_fig_file, bbox_inches="tight")
    #plt.close()

if __name__=="__main__":
    qaplot_sf_inherit(29, sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]})
