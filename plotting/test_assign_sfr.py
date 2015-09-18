
from cenque import CenQue
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from plotting.plot_cenque import plot_cenque_ssfr_dist

def test_assign_sfr(n_snap = 13):
    """
    """
    snap = CenQue(
            n_snap = n_snap, 
            cenque_type = 'sf_assigned'
            ) 
    snap.readin()  
    ssfr_fig = plot_cenque_ssfr_dist(snap, lw=2, line_style='--')      # plot!
    
    fig_file = ''.join(['figure/', 
        'qaplot_assign_sfr_', str(n_snap), '.png'
        ])
    ssfr_fig.savefig(
            fig_file, bbox_inches='tight'
            ) 
    ssfr_fig.clear() 
    plt.close(ssfr_fig)

if __name__=="__main__": 
    test_assign_sfr()
