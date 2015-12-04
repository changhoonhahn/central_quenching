import numpy as np

from cenque import CenQue
from evolve import evolve_cq
from assign_sfr import assign_sfr 
from plotting.plot_tau import plot_tau
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from util.cenque_utility import get_z_nsnap
from plotting.plot_cenque import PlotCenque
from plotting.plot_cenque import plot_ssfr_cenque
from plotting.plot_fq import plot_fqobs_snapshot_evol

def test_quenching_population(
        n_snaps=[12,11,10,9,8,7,6,5,4,3,2,1], 
        tau_prop = {'name': 'instant'}, 
        original_nsnap = 13,
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
        if i_nsnap == original_nsnap: 
            next_snap = CenQue(n_snap = i_nsnap, cenque_type = 'sf_assigned') 
        else: 
            next_snap = CenQue(n_snap = i_nsnap, cenque_type = 'evol_from'+str(original_nsnap)) 

        next_snap.tau_prop = tau_prop
        next_snap.readin()
        
        ssfr_fig.cenque_quenching_ssfr_dist(next_snap)
    
    for i_mass, panel_mass in enumerate(ssfr_fig.panel_mass_bins):       # loop through each panel 

        ssfr_cut = -11.35 + 0.76*(next_snap.zsnap-0.05) - 0.35*((0.5 * np.sum(panel_mass))-10.5)

        ssfr_fig.subs[i_mass].vlines(ssfr_cut, 0.0, 10.0, lw=4)

    ssfr_fig.set_axes()
    
    if tau_prop['name'] == 'line': 
        tau_str = ''.join([
            tau_prop['name'], 
            '_slope', str(tau_prop['slope']), 
            '_yint', str(tau_prop['yint'])
            ])
    else: 
        tau_str = tau_prop['name']

    fig_file = ''.join([
        '/home/users/hahn/research/pro/tinker/central_quenching/figure/', 
        'qaplot_ssfr_quenching_population_nsnap',
        ''.join([str(i_nsnap) for i_nsnap in n_snaps]), 
        '_', tau_str, 
        '.png'
        ]) 
        
    ssfr_fig.fig.savefig(fig_file) 
    ssfr_fig.fig.clear()
    plt.close()

def test_quenching_fraction(
        tau_prop = {'name': 'instant'}, 
        n_snap0 = 13,
        **kwargs
        ): 
    ''' Plot evolution of the quenching population in the 
    CenQue SSFR distribution 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 
    nsnaps : List of snapshot #s to plot  
    '''

    if tau_prop['name'] in ('instant', 'constant', 'satellite', 'long'): 
        tau_str = ''.join(['_', tau_prop['name'], 'tau'])
    elif tau_prop['name'] in ('line'):
        tau_str = ''.join([
            '_', tau_prop['name'], 'tau', 
            '_Mfid', str(tau_prop['fid_mass']), 
            '_slope', str(round(tau_prop['slope'], 4)), 
            '_yint', str(round(tau_prop['yint'],4))
            ])
    prettyplot()
    pretty_color = prettycolors() 

    fig1 = plt.figure(1, figsize=(15, 7))
    sub1 = fig1.add_subplot(111)
    
    fig2 = plt.figure(2, figsize=(12, 7))
    sub2 = fig2.add_subplot(111)

    # Overplot CenQue of specified Snapshots 
    snaps = [] 
    for i_nsnap in reversed(xrange(1, n_snap0)):  
        fqing_file = ''.join([
            '/data1/hahn/f_quenching/', 
            'quenching_fraction', 
            tau_str, 
            '_nsnap', 
            str(i_nsnap), 
            '.dat'
            ]) 

        mass_bin, fqing = np.loadtxt(
                fqing_file, 
                skiprows=1, 
                unpack=True, 
                usecols=[0,1]
                ) 

        #print i_nsnap, 'slope = ', (fqing[-1] - fqing[0])/(mass_bin[-1] - mass_bin[0])
        #print i_nsnap, 'slope = ', (fqing[-2] - fqing[0])/(mass_bin[-2] - mass_bin[0])
        #print i_nsnap, 'slope = ', (fqing[-3] - fqing[0])/(mass_bin[-3] - mass_bin[0])
        #print ((fqing[-1] - fqing[0])/(mass_bin[-1] - mass_bin[0]) + (fqing[-2] - fqing[0])/(mass_bin[-2] - mass_bin[0]) + (fqing[-3] - fqing[0])/(mass_bin[-3] - mass_bin[0]))/3.0

        sub1.plot(mass_bin, fqing, c=pretty_color[i_nsnap+1], lw='4', label='Snapshot'+str(i_nsnap)) 
        #sub1.plot(mass_bin, (0.03 * (np.array(mass_bin) - 9.5))*(1.8 - get_z_nsnap(i_nsnap))**5.0, ls='--', lw='3', c=pretty_color[i_nsnap+1])

        if i_nsnap == n_snap0-1: 
            fqing_m = [] 
            for i_m, mass in enumerate(mass_bin): 
                fqing_m.append([])
    
        for i_m, mass in enumerate(mass_bin): 
            fqing_m[i_m].append(fqing[i_m])

        snaps.append(get_z_nsnap(i_nsnap))
    
    #print fqing_m
    for i_m, mass in enumerate(mass_bin): 
        #print mass, np.around(fqing_m[i_m],4)
        #print 'slope = ', (fqing_m[i_m][-1] - fqing_m[i_m][0])/(snaps[0] - snaps[-1])
        #print 'slope = ', (fqing_m[i_m][-2] - fqing_m[i_m][0])/(snaps[0] - snaps[-2])
        #print 'slope = ', (np.log10(fqing_m[i_m][-1]) - np.log10(fqing_m[i_m][1]))/(snaps[1] - snaps[-1])

        fqing_massbin = np.array(fqing_m[i_m])
        fqing_massbin[np.where(fqing_massbin == 0.)] = 10.**-10
        #print snaps, fqing_massbin 
        sub2.scatter(snaps, fqing_massbin, c=pretty_color[i_m+1])
        sub2.plot(snaps, fqing_massbin, lw='3', c=pretty_color[i_m+1], label=r'$\mathtt{M_* =\;}$'+str(mass))
        #sub2.plot(snaps, (0.03 * (mass - 9.5)) * (1.8 - np.array(snaps))**2.0, ls='--', lw='3', c=pretty_color[i_m+1])

    del snaps
    del fqing_m

    sub1.legend(loc='lower right')
    sub1.set_yscale('log')
    sub1.set_ylim([0.0, 1.0]) 
    sub1.set_xlim([9.5, 14.0])
    sub1.set_xlabel('Stellar Mass ($M_*$)')
    sub1.set_ylabel('Predicted Quenching Fraction')
    
    sub2.legend(loc='lower right')
    sub2.set_yscale('log')
    sub2.set_xlim([0.9, -0.3]) 
    sub2.set_ylim([0.0001, 1.0]) 
    sub2.set_xlabel('Redshift ($\mathtt{z}$)')
    sub2.set_ylabel('Predicted Quenching Fraction')

    fig1.savefig(
            ''.join([
                '/home/users/hahn/research/pro/tinker/central_quenching/figure/', 
                'f_quenching', 
                tau_str, 
                '.png']), 
            bbox_inches='tight'
            )
    fig2.savefig(
            ''.join([
                '/home/users/hahn/research/pro/tinker/central_quenching/figure/', 
                'f_quenching_evol', 
                tau_str, 
                '.png']), 
            bbox_inches='tight'
            )
    fig1.clear()
    fig2.clear()
    plt.close()

if __name__=="__main__": 
    #tau_prop = {'name': 'instant'}
    tau_prop = {'name': 'satellite'}
    #tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}
    #tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.45, 'yint': 0.7}
    #plot_ssfr_cenque(n_snaps=[1], tau_prop = tau_prop, Mrcut=18)
    # 
    #tau_params = [tau_prop]
    #plot_tau(tau_params)

    #tau_prop = {'name': 'satellite'}
    #plot_ssfr_cenque(n_snaps=[1], tau_prop = tau_prop, Mrcut=18)
    for i_snap in reversed(range(1,14)): 
        test_quenching_population(n_snaps=[i_snap], tau_prop = tau_prop, Mrcut=18)
    
    plot_fqobs_snapshot_evol(
            nsnaps = [1,2,3,4,5,6,7,8,9,10,11,12], 
            cenque_type = 'evol_from13', 
            fq_prop = {'name': 'wetzelsmooth'}, 
            sf_prop = {'name': 'average'}, 
            tau_prop= tau_prop
            )
    
    test_quenching_fraction(
        tau_prop = tau_prop, 
        n_snap0 = 13)

