import numpy as np

from sf_inherit import sf_inherit 
from plotting.plot_fq import PlotFq
from plotting.plot_tau import plot_tau
from plotting.plot_cenque import PlotCenque

def abc_posterior(n_step, nsnap_ancestor = 20, nsnap_descendant = 1): 
    """
    Plot SSFR from PMC-ABC derived posterior of tau_Q
    """
    
    # import posterior thetas 
    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta = np.loadtxt(
            posterior_file, 
            unpack = True
            ) 

    ssfr_plot = PlotCenque()

    fq_plot = PlotFq()
    """ 
    random_thetas = np.random.choice(len(theta[0]), 50, replace=False)

    for i_theta in random_thetas: 
        print i_theta

        bloodline = sf_inherit(
                [1], 
                nsnap_ancestor = nsnap_ancestor, 
                pq_prop = {
                    'slope': theta[0][i_theta], 
                    'yint': theta[1][i_theta]
                    }, 
                tau_prop = {
                    'name': 'line', 
                    'fid_mass': 11.1, 
                    'slope': theta[2][i_theta], 
                    'yint': theta[3][i_theta]
                    }, 
                quiet = True
                )

        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))

        ssfr_plot.cenque_ssfr_dist(descendant, line_color='grey', lw=2, label=None)
        fq_plot.cenque_fq(descendant, line_color = 'grey', label = None)
    """ 
    
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
            quiet = True
            )
    descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))

    tau_plot = plot_tau([{
                'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                }])
        
    ssfr_plot.cenque_ssfr_dist(
            descendant, 
            line_color='red', 
            line_width=4, 
            label=r'Median $\mathtt{\tau_{Q}}$ of Posterior'
            )
        
    fq_plot.cenque_fq(descendant, line_color = 'r', label = None)
    
    
    fq_plot.param_fq(
            line_color = 'k', 
            label = None
            )
    fq_plot.set_axes()
    
    ssfr_plot.groupcat_ssfr_dist(Mrcut=18)
    ssfr_plot.set_axes()
    ssfr_fig_file = ''.join([
        'figure/', 
        'ssfr_abc_posterior.png'
        ])
    fq_fig_file = ''.join([
        'figure/', 
        'fq_abc_posterior.png'
        ])
    tau_fig_file = ''.join([
        'figure/', 
        'tau_abc_posterior.png'
        ])

    ssfr_plot.fig.savefig(ssfr_fig_file, bbox_inches="tight")
    fq_plot.fig.savefig(fq_fig_file, bbox_inches="tight")
    tau_plot.savefig(tau_fig_file, bbox_inches="tight")
    #plt.show()

def abc_post_fq(n_step, descendant_nsnaps, nsnap_ancestor = 20): 
    """
    Plot SSFR from PMC-ABC derived posterior of tau_Q
    """
    
    # import posterior thetas 
    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta = np.loadtxt(
            posterior_file, 
            unpack = True
            ) 

    ssfr_plot = PlotCenque()

    fq_plot = PlotFq()
    
    # median theta 
    med_theta = [] 
    for i_param in xrange(len(theta)):
        med_theta.append( 
                np.median(theta[i_param])
                )
        
    bloodline = sf_inherit(
            descendant_nsnaps, 
            nsnap_ancestor = nsnap_ancestor, 
            pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
            tau_prop = {
                'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                }, 
            quiet = True
            )

    for nsnap in descendant_nsnaps: 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap))
        
        if nsnap != 1: 
            ssfr_plot.cenque_ssfr_dist(
                    descendant, 
                    lw=2, 
                    line_color = ssfr_plot.pretty_colors[descendant.nsnap], 
                    label=None
                    )

            fq_plot.cenque_fq(
                    descendant, 
                    lw=2, 
                    label = 't=' + str(round(descendant.t_cosmic, 1))+',z='+str(round(descendant.zsnap, 2))  
                    )
            fq_plot.param_fq(lw=2, line_color = fq_plot.pretty_colors[descendant.nsnap], label = None)

        else: 
            ssfr_plot.cenque_ssfr_dist(descendant, line_color='r')
            fq_plot.cenque_fq(
                    descendant, 
                    label = 't=' + str(round(descendant.t_cosmic, 1))+',z='+str(round(descendant.zsnap, 2))  
                    )
            fq_plot.param_fq(
                    line_color = fq_plot.pretty_colors[descendant.nsnap], 
                    label = None
                    )

    tau_plot = plot_tau([{
                'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                }])
            
    ssfr_plot.groupcat_ssfr_dist(Mrcut=18)
    fq_plot.set_axes()
    ssfr_plot.set_axes()
    #ssfr_fig_file = ''.join([
    #    'figure/', 
    #    'ssfr_abc_posterior.png'
    #    ])
    #fq_fig_file = ''.join([
    #    'figure/', 
    #    'fq_abc_posterior.png'
    #    ])
    tau_fig_file = ''.join([
        'figure/', 
        'tau_abc_posterior.png'
        ])

    #ssfr_plot.fig.savefig(ssfr_fig_file, bbox_inches="tight")
    #fq_plot.fig.savefig(fq_fig_file, bbox_inches="tight")
    tau_plot.savefig(tau_fig_file, bbox_inches="tight")
    plt.show()

if __name__=="__main__": 
    #abc_post_fq(29, [1,2,3,6,9,12,15,18], nsnap_ancestor = 20)
    abc_posterior(29)
