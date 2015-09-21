
from cenque import CenQue
from evolve import evolve_cq
from assign_sfr import assign_sfr 
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from plotting.plot_cenque import plot_cenque_ssfr_dist

def test_evolve(nsnap = 13, final_nsnap=1, clobber=False):
    """
    """
    if clobber: 
        snap = CenQue()
        snap.import_treepm(nsnap)
        snap.writeout()
        snap = assign_sfr(snap)
        snap.writeout()
        snap = evolve_cq(snap, quiet=True)
    else:
        snap = CenQue(n_snap = final_nsnap, cenque_type='evol_from'+str(nsnap))
        snap.readin()
    
    ssfr_fig = plot_cenque_ssfr_dist(snap, lw=2, line_style='--')      # plot!

    plt.show()
     
    #fig_file = ''.join(['figure/', 
    #    'qaplot_evolve_', str(n_snap), '.png'
    #    ])
    #ssfr_fig.savefig(
    #        fig_file, bbox_inches='tight'
    #        ) 
    #ssfr_fig.clear() 
    #plt.close(ssfr_fig)

if __name__=="__main__": 
    test_evolve()
