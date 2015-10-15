import numpy as np

from evolve_lineage import evolve_lineage
from plotting.plot_cenque import PlotCenque

def abc_post_ssfr(n_step, nsnap_ancestor = 13, nsnap_descendant = 1): 
    """
    Plot SSFR from PMC-ABC derived posterior of tau_Q
    """
    
    # import posterior thetas 
    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta_slope, theta_yint = np.loadtxt(
            posterior_file, 
            unpack = True, 
            usecols = [0,1]
            ) 

    ssfr_plot = PlotCenque()

    random_thetas = np.random.choice(len(theta_slope), 50, replace=False)

    for i_theta in random_thetas: 
        print i_theta
        evol_cq = evolve_lineage(
            nsnap_ancestor = 13, 
            nsnap_descendant = 1, 
            tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': theta_slope[i_theta], 'yint': theta_yint[i_theta]}
            )
        
        ssfr_plot.cenque_ssfr_dist(evol_cq, line_color='grey', lw=2, label=None)

    med_theta_slope = np.median(theta_slope)
    med_theta_yint = np.median(theta_yint)

    evol_cq = evolve_lineage(
        nsnap_ancestor = 13, 
        nsnap_descendant = 1, 
        tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': med_theta_slope, 'yint': med_theta_yint}
        )
        
    ssfr_plot.cenque_ssfr_dist(evol_cq, line_color='red', line_width=4, label=r'Median $\mathtt{\tau_{Q}}$ of Posterior')
    
    evol_cq = evolve_lineage(
        nsnap_ancestor = 13, 
        nsnap_descendant = 1, 
        tau_prop = {'name': 'satellite'}
        )
        
    ssfr_plot.cenque_ssfr_dist(evol_cq, line_color='blue', lw=2, label=r'Satellite $\mathtt{\tau_{Q}}$')

    ssfr_plot.groupcat_ssfr_dist(Mrcut=18)
    ssfr_plot.set_axes()
    fig_file = ''.join([
        'figure/', 
        'ssfr_abc_posterior_tau.png'])

    ssfr_plot.fig.savefig(fig_file, bbox_inches="tight")
    #plt.show()

if __name__=="__main__": 
    abc_post_ssfr(18)
