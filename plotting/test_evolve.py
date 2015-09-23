
from cenque import CenQue
from evolve import evolve_cq
from assign_sfr import assign_sfr 
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from plotting.plot_cenque import PlotCenque

def test_evolve(n_snaps=[12,11,10,9,8,7,6,5,4,3,2,1], Mrcut=18, **kwargs): 
    ''' Plot evolution of the CenQue SSFR distribution 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 
    nsnaps : List of snapshot #s to plot  
    '''

    # Plot original CenQue object with assigned SFR 
    snap = CenQue(n_snap=13, cenque_type = 'sf_assigned') 
    snap.readin()  

    ssfr_fig = PlotCenque(cenque=snap, lw=2, linestyle='--')
    
    # Overplot CenQue of specified Snapshots 
    for i_nsnap in n_snaps:  
        next_snap = CenQue(n_snap = i_nsnap, cenque_type = 'evol_from13') 
        next_snap.readin()
        
        ssfr_fig.cenque_ssfr_dist(next_snap)
    
    if Mrcut != None: 
        ssfr_fig.groupcat_ssfr_dist(Mrcut=Mrcut)
    
    for i_mass, panel_mass in enumerate(ssfr_fig.panel_mass_bins):       # loop through each panel 

        ssfr_cut = -11.35 + 0.76*(next_snap.zsnap-0.05) - 0.35*((0.5 * np.sum(panel_mass))-10.5)

        ssfr_fig.subs[i_mass].vlines(ssfr_cut, 0.0, 10.0, lw=4)

    ssfr_fig.set_axes()

    plt.show()

def test_quenching_population(
        n_snaps=[12,11,10,9,8,7,6,5,4,3,2,1], 
        tau_prop = {'name': 'instant'}, 
        **kwargs
        ): 
    ''' Plot evolution of the quenching population in the 
    CenQue SSFR distribution 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 
    nsnaps : List of snapshot #s to plot  
    '''

    ssfr_fig = PlotCenque()
    
    # Overplot CenQue of specified Snapshots 
    for i_nsnap in n_snaps:  
        next_snap = CenQue(n_snap = i_nsnap, cenque_type = 'evol_from13') 
        next_snap.tau_prop = tau_prop
        next_snap.readin()
        
        ssfr_fig.cenque_quenching_ssfr_dist(next_snap)
    
    for i_mass, panel_mass in enumerate(ssfr_fig.panel_mass_bins):       # loop through each panel 

        ssfr_cut = -11.35 + 0.76*(next_snap.zsnap-0.05) - 0.35*((0.5 * np.sum(panel_mass))-10.5)

        ssfr_fig.subs[i_mass].vlines(ssfr_cut, 0.0, 10.0, lw=4)

    ssfr_fig.set_axes()

    fig_file = ''.join([
        '/home/users/hahn/research/pro/tinker/central_quenching/figure/', 
        'qaplot_quenching_pop_',
        ''.join([str(i_nsnap) for i_nsnap in n_snaps]), 
        '.png'
        ]) 

    ssfr_fig.fig.savefig(fig_file) 
    ssfr_fig.fig.clear()
    plt.close()

if __name__=="__main__": 
    for i_snap in [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]: 
        test_quenching_population(n_snaps=[i_snap], tau_prop = {'name': 'satellite'})
