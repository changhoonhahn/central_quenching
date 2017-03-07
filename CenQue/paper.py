'''

Calculations and figures pertaining to the 
central quenching timescale paper

'''
import os 
import pickle
import numpy as np 
import astropy.cosmology as astrocosmo
from scipy.interpolate import interp1d

import util.util as Util 
import sfr_evol
from inherit import Inherit
from abcee import PlotABC
from abcee import ReadABCrun

from gal_prop import Fq 
from gal_prop import SMF
from observations import GroupCat 
from observations import FqCen_bestfit

import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def fig_SMFevol(): 
    ''' Plot the evolution of the SMF used for SHAM and also plot the 
    extreme cases.
    '''
    zrange = [0.9, 0.5, 0.05]
    source_list = ['li-march', 'li-march-extreme']

    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(7, 7))
    sub = fig.add_subplot(111)

    smf = SMF()
    
    for zz in np.arange(0.9, 1.7, 0.1): 
        m, phi1 = smf.analytic(zz, source='li-march-extreme')
        m, phi2 = smf.analytic(zz, source='li-march')
        print zz, phi1/phi2

    for i_s, source in enumerate(source_list): 
        for iz, z in enumerate(zrange): 
            mass, phi = smf.analytic(z, source=source)  # SMF phi(M*) 
                
            # labels for the SMF sources 
            if source == 'li-march': 
                source_lbl = 'Marchesini+2009 \&\nLi+2009 (Fiducial)'
            elif source == 'li-march-extreme': 
                #source_lbl = 'Li$+$2009 $\&$ \nextreme evolution'
                source_lbl = 'extreme SMF evolution'

            if source == 'li-march':  
                if z == 0.05: 
                    sub.plot(mass, phi, lw=4, color='#b35a00', label=source_lbl)
                    sub.text(10.675, 0.01, "z", color='#b35a00', fontsize=25)
                    sub.text(10.8, 0.01, r"$\mathtt{= 0.05}$", color='#b35a00', fontsize=25)
                elif z == 0.5: 
                    sub.plot(mass, phi, lw=4, color='#ff8100')
                    sub.text(10.8, 0.75*0.01, r"$\mathtt{= 0.5}$", color='#ff8100', fontsize=25)
                elif z == 0.9: 
                    sub.plot(mass, phi, lw=4, color='#ffa74d') 
                    sub.text(10.8, 0.75*0.75*0.01, r"$\mathtt{= 0.9}$", color='#ffa74d', fontsize=25)
            else: 
                if z == 0.5:  
                    sub.plot(mass, phi, lw=2, ls='-.', c='k', label=source_lbl)
                elif z == 0.9:  
                    sub.plot(mass, phi, lw=2, ls='-.', c='k')
    
    sub.set_ylim([10**-3.75, 10**-1.75])
    sub.set_xlim([8.8, 11.5])
    sub.set_xticks([9.0, 10.0, 11.0])
    sub.minorticks_on()
    sub.set_yscale('log')
    sub.set_xlabel(r'log($\mathtt{M_*} /\mathtt{M_\odot}$)', fontsize=25) 
    sub.set_ylabel(r'log($\mathtt{\Phi / Mpc^{-3}\;dex^{-1}}$)', fontsize=25) 
    sub.legend(loc='lower left', scatteryoffsets=[0.6])
    
    fig_file = ''.join(['figure/paper/', 'SMF_evol.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    Util.png2pdf(fig_file)
    return None


def fig_SFMSevol(): 
    ''' Plot the evolution of the SMF used for SHAM and also plot the 
    extreme cases.
    '''
    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(7, 6))
    sub = fig.add_subplot(111)
    
    z_arr = np.array([0.1, 0.3, 0.5, 0.7, 0.9]) 

    sfms_dict = {'name': 'linear', 'zslope': 1.14}
    mu_sfr_z = sfr_evol.AverageLogSFR_sfms(np.array([10.5]), z_arr, sfms_prop=sfms_dict) 
    print mu_sfr_z 
    sig_sfr_z = sfr_evol.ScatterLogSFR_sfms(np.array([10.5]), z_arr, sfms_prop=sfms_dict)

    # Lee et al. (2015) 
    # log SFR(M*, z) = S0(z) - log( 1 + (M*/M0)^-gamma)
    z_mid = np.array([0.36, 0.55, 0.70, 0.85, 0.99, 1.19])
    S0 = np.array([0.80, 0.99, 1.23, 1.35, 1.53, 1.72])
    S0_err = np.array([0.019, 0.015, 0.016, 0.014, 0.017, 0.024])
    M0 = np.array([10.03, 9.82, 9.93, 9.96, 10.10, 10.31]) 
    gamma = np.array([0.92, 1.13, 1.11, 1.28, 1.26, 1.07]) 

    sfr_lee = lambda mm: S0 - np.log10(1.+np.power(10, mm - M0)**(-1. * gamma))

    sub.errorbar(z_arr, mu_sfr_z, yerr=sig_sfr_z, color='#b35a00')
    sub.scatter(z_mid, sfr_lee(10.5), color='k')
            
    #sub.set_ylim([10**-5, 10**-1.5])
    sub.set_xlim([0., 1.0])
    #sub.set_yscale('log')
    #sub.set_xlabel(r'log($\mathtt{M_*} /\mathtt{M_\odot}$)', fontsize=25) 
    #sub.set_ylabel(r'log($\mathtt{\Phi / Mpc^{-3}\;dex^{-1}}$)', fontsize=25) 
    sub.legend(loc='lower left', scatteryoffsets=[0.6])
    plt.show()
    
    #fig_file = ''.join(['figure/paper/', 'SMF_evol.png'])
    #fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    #Util.png2pdf(fig_file)
    return None


def fig_fQcen_evol(): 
    ''' Plot the evolution of the central galaxy quiescent fraction
    '''
    # best fit alpha(M*) values 
    m_mid, alpha_m = FqCen_bestfit(clobber=True) 

    # mass binnning we impose 
    m_bin = np.arange(9.5, 12.5, 0.5)
    if not np.array_equal(m_mid, 0.5 * (m_bin[1:] + m_bin[:-1])):
        raise ValueError

    # read in SDSS fQ^cen data from Jeremy and evaulate at M_mid
    fq_file = ''.join(['dat/observations/cosmos_fq/', 'fcen_red_sdss_scatter.dat']) 
    m_sdss, fqcen_sdss, N_sdss = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2])
    
    fqcen_sdss_rebin = np.zeros(len(m_bin)-1)  
    for im in xrange(len(m_bin)-1): 
        sdss_mbin = np.where((m_sdss >= m_bin[im]) & (m_sdss < m_bin[im+1])) 
    
        fqcen_sdss_rebin[im] = np.sum(fqcen_sdss[sdss_mbin] * N_sdss[sdss_mbin].astype('float'))/np.sum(N_sdss[sdss_mbin].astype('float'))
    
    # Read in COSMOS fQ^cen data from Jeremy
    fqcen_cosmos_rebin = [] 
    fqcen_low_cosmos_rebin = [] 
    fqcen_high_cosmos_rebin = [] 
    for iz, z in enumerate([0.36, 0.66, 0.88]): 
        fq_file = ''.join(['dat/observations/cosmos_fq/', 
            'stats_z', str(iz+1), '.fq_cen']) 
        
        m_cosmos, fqcen_cosmos, fqcen_cosmos_low, fqcen_cosmos_high = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2,3])
        m_cosmos = np.log10(m_cosmos)

        fqcen_interp = interp1d(m_cosmos, fqcen_cosmos) 
        fqcen_low_interp = interp1d(m_cosmos, fqcen_cosmos_low) 
        fqcen_high_interp = interp1d(m_cosmos, fqcen_cosmos_high) 
    
        fqcen_cosmos_rebin.append( fqcen_interp(m_mid) ) 
        fqcen_low_cosmos_rebin.append( fqcen_low_interp(m_mid) ) 
        fqcen_high_cosmos_rebin.append( fqcen_high_interp(m_mid) ) 

    fqcen_cosmos_rebin = np.array(fqcen_cosmos_rebin) 
    fqcen_low_cosmos_rebin = np.array(fqcen_low_cosmos_rebin)
    fqcen_high_cosmos_rebin = np.array(fqcen_high_cosmos_rebin) 

    z_arr = np.array([0.025, 0.36, 0.66, 0.88])
    z_fin = np.arange(0.0, 1.1, 0.1)
    
    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(7, 7))
    sub = fig.add_subplot(111)#, len(m_mid), im+1) 
    lines = []
    for im in range(len(m_mid)-1): 
        fqcen_z = np.array([fqcen_sdss_rebin[im]]+list(fqcen_cosmos_rebin[:,im]))
        fqcen_z_err_low = np.array(
                [0.0] + list(fqcen_cosmos_rebin[:,im] - fqcen_low_cosmos_rebin[:,im])
                )
        fqcen_z_err_high = np.array(
                [0.0] + list(fqcen_high_cosmos_rebin[:,im] - fqcen_cosmos_rebin[:,im])
                )
        if im == 0: 
            err_label = 'Tinker+(2013)'
            tink = sub.errorbar(z_arr, fqcen_z, yerr=[fqcen_z_err_low, fqcen_z_err_high], 
                    fmt='o', markersize=8, color=pretty_colors[2*im+1], lw=0, 
                    ecolor='k', elinewidth=2, capthick=2, label=err_label)
            tink_handles, tink_labels = sub.get_legend_handles_labels()
            # remove the errorbars
            tink_handles = [h[0] for h in tink_handles]
        else: 
            err_label = None
            sub.errorbar(z_arr, fqcen_z, yerr=[fqcen_z_err_low, fqcen_z_err_high], 
                    fmt='o', markersize=8, color=pretty_colors[2*im+1], lw=0,
                    ecolor='k', elinewidth=2, capthick=2)

        line, = sub.plot(z_fin, fqcen_sdss_rebin[im] * (1. + z_fin)**alpha_m[im], 
                c=pretty_colors[2*im+1], ls='--', lw=2, 
                label=r'log$\mathtt{\;M_* = '+str(m_mid[im])+'}$')
        lines.append(line) 

    first_legend = sub.legend(handles=lines, loc='upper right', handletextpad=0.1, 
            prop={'size': 20})

    ax = plt.gca().add_artist(first_legend)

    # use them in the legend
    sub.legend(tink_handles, tink_labels, loc='upper left', numpoints=1, markerscale=2, 
            handletextpad=-0.15) 
    
    sub.set_xlim([0.0, 1.0])
    sub.set_ylim([0.0, 1.0])
    #sub.set_ylim([0.01, 1.0]) 
    #sub.set_yscale('log') 
    sub.set_xlabel(r'Redshift $(\mathtt{z})$', fontsize=25) 
    sub.set_ylabel(r'Quiescent Fraction of Central Galaxies', fontsize=22) 

    fig_file = ''.join(['figure/paper/', 'fqcen.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    Util.png2pdf(fig_file)
    plt.close()
    return None


def fig_SFRassign(t, abcrun, prior_name='try0'): 
    ''' Demonstrate the SFR assignment scheme based through the 
    SSFR distribution function P(SSFR) with different GV prescriptions.
    '''
    abc_plot = PlotABC(t, abcrun=abcrun, prior_name=prior_name) 
    # median theta values 
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = abc_plot.med_theta

    no_gv = [0., 0.]
    med_gv = [gv_slope, gv_offset]
    if gv_offset < 0: 
        lot_gv = [2. * gv_slope, np.abs(gv_offset)]
    else: 
        lot_gv = [2. * gv_slope, 2. * gv_offset]

    # other parameters
    sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
    sim_kwargs = sfinherit_kwargs.copy()
    sim_kwargs['evol_prop']['fudge'] = {
            'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
    sim_kwargs['evol_prop']['tau'] = {
            'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

    figdata_file = lambda gv_str: ''.join(['/data1/hahn/paper/', 'fig_SFRassign.', gv_str, '.data_file.p']) 
    
    gv_dict = {} 
    # loop through different types of GV
    for i_gv, gv_list in enumerate([no_gv, med_gv, lot_gv]): 
        # green valley perscription
        if i_gv == 0: 
            gvstr = 'no_gv'
        elif i_gv == 1: 
            gvstr = 'med_gv'
        elif i_gv == 2:  
            gvstr = 'lot_gv'

        if not os.path.isfile(figdata_file(gvstr)): 
            sim_kwargs['sfr_prop']['gv'] = {
                    'slope': gv_list[0], 'fidmass': 10.5, 'offset': gv_list[1]}
            inh = Inherit([1], 
                    nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                    subhalo_prop=sim_kwargs['subhalo_prop'], 
                    sfr_prop=sim_kwargs['sfr_prop'], 
                    evol_prop=sim_kwargs['evol_prop'])
            anc = inh.ancestor
            gv_dict[gvstr] = anc

            pickle.dump(anc, open(figdata_file(gvstr), 'wb'))
        else:
            anc = pickle.load(open(figdata_file(gvstr), 'rb'))
            gv_dict[gvstr] = anc

    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(7, 7))
    sub = fig.add_subplot(111)
    for gv_key in gv_dict.keys():
        if gv_key == 'no_gv': 
            lstyle = '--'
            label = 'No green valley'
            colors = pretty_colors[1]
        elif gv_key == 'med_gv': 
            lstyle = '-' 
            if med_gv[1] > 0: 
                label = r'$\mathtt{f_{GV} = '+str(round(med_gv[0], 1))+' (log\;M_* - 10.5) + '+str(round(med_gv[1], 1))+'}$'
            else:
                label = r'$\mathtt{f_{GV} = '+str(round(med_gv[0], 1))+' (log\;M_* - 10.5) - '+str(round(np.abs(med_gv[1]), 1))+'}$'
            colors = 'k'#pretty_colors[4]
        elif gv_key == 'lot_gv': 
            lstyle = '--'
            label = r'$\mathtt{f_{GV} = '+str(round(lot_gv[0], 1))+' (log\;M_* - 10.5) + '+str(round(lot_gv[1], 1))+'}$'
            colors = pretty_colors[7]

        bin_mid, ssfr_dist = gv_dict[gv_key].Ssfr()
        
        if gv_key == 'no_gv': 
            dashes = [5,2,10,5]
            ssfr_line, = sub.plot(bin_mid[2], ssfr_dist[2], c=colors, lw=3, ls=lstyle, label=label)
            ssfr_line.set_dashes(dashes)
        else: 
            sub.plot(bin_mid[2], ssfr_dist[2], c=colors, lw=3, ls=lstyle, label=label)

    sub.set_xlim([-13.2, -8.0])
    sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
    sub.set_ylim([0.0, 1.6])
    sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
    sub.legend(loc='upper left', prop={'size': 17}, borderpad=1.25)

    fig_file = ''.join(['figure/paper/', 'assignSFR_SSFR.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    return None 


def figSFH_demo(t, abcrun, prior_name='try0'): 
    ''' Figure that demonstrates the Star-formation history 
    of Quiescent, SF, and GV galaxies. These are not actually 
    from the real model but demonstrations
    '''
    abc_plot = PlotABC(t, abcrun=abcrun, prior_name=prior_name) 
    # median theta values 
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = abc_plot.med_theta
    # other parameters
    sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
    sim_kwargs = sfinherit_kwargs.copy()
    sim_kwargs['sfr_prop']['gv'] = {
            'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
    sim_kwargs['evol_prop']['fudge'] = {
            'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
    sim_kwargs['evol_prop']['tau'] = {
            'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

    figdata_file = ''.join(['/data1/hahn/paper/', 'figSFH_demo.data_file.p']) 
    if not os.path.isfile(figdata_file): 
        # run model for median theta values 
        inh = Inherit([1], 
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        des_dict = inh()
        anc = inh.ancestor
        pickle.dump([anc, des_dict], open(figdata_file, 'wb'))
    else:
        anc, des_dict = pickle.load(open(figdata_file, 'rb'))

    z_snap, t_snap = np.loadtxt(Util.snapshottable(), unpack=True, usecols=[2, 3]) 
    z_cosmics = z_snap[1:sim_kwargs['nsnap_ancestor']+1][::-1]
    t_cosmics = t_snap[1:sim_kwargs['nsnap_ancestor']+1][::-1]

    des1 = des_dict['1']    # descendant at nsnap = 1
    q_a_index1 = des1.will[des1.q_ancestor]   
    sf_a_index1 = des1.will[des1.sf_ancestor]
   
    # quiescent 
    onlyone = 0 
    while onlyone == 0: 
        i_ran = np.random.choice(q_a_index1)
        q_Msham_i = (anc.Msham_evol[0])[i_ran][::-1]
        q_Msham_i, q_indices = keep_non_descreasing(list(q_Msham_i))
        q_Msham_i = np.array(q_Msham_i) 

        if q_Msham_i.min() > 9. and len(q_Msham_i) > 9 and q_indices[-1] == len(z_cosmics)-1:  
            q_SFRs_i = anc.ssfr[i_ran] + q_Msham_i
            onlyone+= 1

    # star-forming 
    onlyone = 0 
    while onlyone == 0: 
        i_ran = np.random.choice(sf_a_index1)
        sf_Msham_i = (anc.Msham_evol[0])[i_ran][::-1]
        sf_Msham_i, sf_indices = keep_non_descreasing(list(sf_Msham_i))
        sf_Msham_i = np.array(sf_Msham_i)

        dutycycle_prop = sfr_evol.dutycycle_param(1, 
                dutycycle_prop={'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3})
        dutycycle_prop['delta_sfr'] = 0.3

        if sf_Msham_i.min() > 9.5 and len(sf_Msham_i) > 9 and sf_indices[-1] == len(z_cosmics)-1:  
            sf_SFRs_i = np.zeros(len(z_cosmics)) 
            for i_z, z_f in enumerate(z_cosmics): 
                # log(SFR)_SFMS evolution from t0
                logsfr_sfms = \
                        sfr_evol.AverageLogSFR_sfms(sf_Msham_i[0], z_cosmics.max(), sfms_prop={'name': 'linear', 'zslope': 1.5}) + \
                        sfr_evol.DeltaLogSFR_sfms(z_cosmics.max(), z_f, sfms_prop={'name': 'linear', 'zslope': 1.5}) 

                # log(SFR)_duty cycle evolution from t0 to tQ
                logsfr_sfduty = sfr_evol.DeltaLogSFR_dutycycle(
                        t_cosmics.min(), 
                        t_cosmics[i_z], 
                        t_q=[999.], 
                        dutycycle_prop=dutycycle_prop
                        )
                logsfr_tot = logsfr_sfms + logsfr_sfduty
                sf_SFRs_i[i_z] = logsfr_tot[0]

            sf_SFRs_i = sf_SFRs_i[sf_indices]
            onlyone += 1
    # quenching  
    onlyone = 0 
    tqq = np.array([7.])
    while onlyone == 0.: 
        i_ran = np.random.choice(sf_a_index1)
        qing_Msham_i = (anc.Msham_evol[0])[i_ran][::-1]
        qing_Msham_i, qing_indices = keep_non_descreasing(list(qing_Msham_i))
        qing_Msham_i = np.array(qing_Msham_i)

        dutycycle_prop = sfr_evol.dutycycle_param(1, 
                dutycycle_prop={'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3})
        dutycycle_prop['delta_sfr'] = -0.1

        if qing_Msham_i.min() > 9.0 and qing_Msham_i[0] < 10. and len(qing_Msham_i) > 7 and qing_indices[-1] == len(z_cosmics)-1:  
            qing_SFRs_i = np.zeros(len(z_cosmics))
            for i_z, z_f in enumerate(z_cosmics): 
                # log(SFR)_SFMS evolution from t0
                logsfr_sfms = \
                        sfr_evol.AverageLogSFR_sfms(qing_Msham_i[0], z_cosmics.max(), sfms_prop={'name': 'linear', 'zslope': 1.5}) + \
                        sfr_evol.DeltaLogSFR_sfms(z_cosmics.max(), z_f, sfms_prop={'name': 'linear', 'zslope': 1.5}) 

                # log(SFR)_duty cycle evolution from t0 to tQ
                logsfr_sfduty = sfr_evol.DeltaLogSFR_dutycycle(
                        t_cosmics.min(), 
                        t_cosmics[i_z], 
                        t_q=tqq, 
                        dutycycle_prop=dutycycle_prop
                        )
                closest_tQ_index = np.abs(t_cosmics - tqq[0]).argmin() - 1
                logsfr_quench = sfr_evol.DeltaLogSFR_quenching(
                        tqq, 
                        t_cosmics[i_z],
                        M_q=qing_Msham_i[closest_tQ_index],
                        tau_prop={'name': 'line', 'fid_mass': 11.1, 'slope': -0.6, 'yint': 0.5})
                logsfr_tot = logsfr_sfms + logsfr_sfduty + logsfr_quench
                qing_SFRs_i[i_z] = logsfr_tot[0]

            qing_SFRs_i = qing_SFRs_i[qing_indices]
            onlyone += 1

    # overquenching 
    avg_q_ssfr = sfr_evol.AverageLogSSFR_q_peak(qing_Msham_i[-1])
    sigma_q_ssfr = sfr_evol.ScatterLogSSFR_q_peak(qing_Msham_i[-1])
    min_q_ssfr = (sigma_q_ssfr * np.random.randn(1) + avg_q_ssfr)[0]
    qing_ssfr_i = np.array(qing_SFRs_i) - qing_Msham_i
    qing_SFRs_i[np.where(qing_ssfr_i < min_q_ssfr)] = qing_Msham_i[np.where(qing_ssfr_i < min_q_ssfr)] + min_q_ssfr

    # figure 
    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(14, 6))
    sub0 = fig.add_subplot(121)

    sub0.quiver(q_Msham_i[:-1], q_SFRs_i[:-1], q_Msham_i[1:]-q_Msham_i[:-1], q_SFRs_i[1:]-q_SFRs_i[:-1], 
            color=pretty_colors[4], scale_units='xy', angles='xy', scale=1)
    sub0.scatter([0], [0], c=pretty_colors[4], lw=0, s=1, label='Quiescent')

    sub0.quiver(sf_Msham_i[:-1], sf_SFRs_i[:-1], sf_Msham_i[1:]-sf_Msham_i[:-1], sf_SFRs_i[1:]-sf_SFRs_i[:-1], 
            color=pretty_colors[1], scale_units='xy', angles='xy', scale=1)
    sub0.scatter([0], [0], c=pretty_colors[1], lw=0, s=1, label='Star-Forming')
    
    sub0.quiver(qing_Msham_i[:-1], qing_SFRs_i[:-1], qing_Msham_i[1:]-qing_Msham_i[:-1], qing_SFRs_i[1:]-qing_SFRs_i[:-1], 
            color=pretty_colors[6], scale_units='xy', angles='xy', scale=1)
    sub0.scatter([0], [0], c=pretty_colors[6], lw=0, s=1, label='Quenching')

    plt.sca(sub0)
    plt.xticks(np.arange(9.0, 12.0, 0.5))
    sub0.set_xlim([9.0, 12.0]) 
    sub0.set_xlabel(r'$\mathtt{log}(\mathtt{M_*}\; [\mathtt{M}_\odot])$', fontsize=25) 

    sub0.set_ylim([-5., 2.0]) 
    sub0.set_ylabel(r'$\mathtt{log}(\mathtt{SFR}\;[\mathtt{M}_\odot/\mathtt{yr}])$', fontsize=25) 
    sub0.legend(loc='lower right', scatterpoints=1, prop={'size': 20}, borderpad=0.5, handletextpad=-0.25, 
            markerscale=10.) 

    sub1 = fig.add_subplot(122)
    q_t_cosmics = t_cosmics[q_indices]
    sub1.quiver(q_t_cosmics[:-1], q_SFRs_i[:-1], q_t_cosmics[1:]-q_t_cosmics[:-1], q_SFRs_i[1:]-q_SFRs_i[:-1], 
            color=pretty_colors[4], scale_units='xy', angles='xy', scale=1)
    sf_t_cosmics = t_cosmics[sf_indices]
    sub1.quiver(sf_t_cosmics[:-1], sf_SFRs_i[:-1], sf_t_cosmics[1:]-sf_t_cosmics[:-1], sf_SFRs_i[1:]-sf_SFRs_i[:-1], 
            color=pretty_colors[1], scale_units='xy', angles='xy', scale=1)
    qing_t_cosmics = t_cosmics[qing_indices]
    sub1.quiver(qing_t_cosmics[:-1], qing_SFRs_i[:-1], qing_t_cosmics[1:]-qing_t_cosmics[:-1], qing_SFRs_i[1:]-qing_SFRs_i[:-1], 
            color=pretty_colors[6], scale_units='xy', angles='xy', scale=1)
    #sub1.plot(t_cosmics, q_SFRs[i,:], c=pretty_colors[4])
    
    sub1.set_xlim([5., 14.])
    sub1.set_xlabel(r'$\mathtt{t_{cosmic}}\;[\mathtt{Gyr}]$', fontsize=25) 

    sub1.set_yticklabels([])
    sub1.set_ylim([-5., 2.0]) 

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/paper/', 'SFH_demo.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    plt.close() 
    return None 


def figSFH_SchematicDemo(t, abcrun, prior_name='try0'): 
    ''' Figure that schematically demonstrates the Star-formation history 
    of Quiescent, SF, and GV galaxies. These are not actually 
    from the real model but demonstrations
    '''
    abc_plot = PlotABC(t, abcrun=abcrun, prior_name=prior_name) 
    # median theta values 
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = abc_plot.med_theta
    # other parameters
    sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
    sim_kwargs = sfinherit_kwargs.copy()
    sim_kwargs['sfr_prop']['gv'] = {
            'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
    sim_kwargs['evol_prop']['fudge'] = {
            'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
    sim_kwargs['evol_prop']['tau'] = {
            'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

    figdata_file = ''.join(['/data1/hahn/paper/', 'figSFH_demo.data_file.', abcrun, '.p']) 
    if not os.path.isfile(figdata_file): 
        # run model for median theta values 
        inh = Inherit([1], 
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        des_dict = inh()
        anc = inh.ancestor
        pickle.dump([anc, des_dict], open(figdata_file, 'wb'))
    else:
        anc, des_dict = pickle.load(open(figdata_file, 'rb'))

    z_snap, t_snap = np.loadtxt(Util.snapshottable(), unpack=True, usecols=[2, 3]) 
    z_cosmics = z_snap[1:sim_kwargs['nsnap_ancestor']+1][::-1]
    t_cosmics = t_snap[1:sim_kwargs['nsnap_ancestor']+1][::-1]
    #z_cosmics = np.arange(z_cosmics.min(), z_cosmics.max()+0.01, 0.01) 
    #t_cosmics = Util.get_tsnap(z_cosmics)

    des1 = des_dict['1']    # descendant at nsnap = 1

    M_anc = (anc.mass)[des1.will] # ancestor mass 
    m0_lim = np.where((M_anc > 9.5) & (M_anc < 10.)) 
    
    Msham_evol = (anc.Msham_evol[0])[des1.will[m0_lim]][::-1]

    avg_Msham_evol = np.zeros(Msham_evol.shape[1]) 
    for i in range(Msham_evol.shape[1]): 
        hasmass = np.where(Msham_evol[:,i] > 0.)
        avg_Msham_evol[i] = np.mean((Msham_evol[:,i])[hasmass])

    # star-forming 
    sf_SFRs = \
            sfr_evol.AverageLogSFR_sfms(
                    avg_Msham_evol,
                    z_cosmics.max(), 
                    sfms_prop=sim_kwargs['sfr_prop']['sfms']) + \
            sfr_evol.DeltaLogSFR_sfms(
                    z_cosmics.max(), 
                    z_cosmics.min(), 
                    sfms_prop=sim_kwargs['sfr_prop']['sfms']) 

    # quenching  
    tqq = np.array([9.0])
    qing_SFRs = np.zeros(len(z_cosmics))
    qing_SFRs_sat = np.zeros(len(z_cosmics))
    for i_z, z_f in enumerate(z_cosmics): 
        # log(SFR)_SFMS evolution from t0
        logsfr_sfms = sf_SFRs[i_z]

        closest_tQ_index = np.abs(t_cosmics - tqq[0]).argmin() - 1
        logsfr_quench = sfr_evol.DeltaLogSFR_quenching(
                tqq, 
                t_cosmics[i_z],
                M_q=avg_Msham_evol[closest_tQ_index],
                tau_prop=sim_kwargs['evol_prop']['tau'])
        logsfr_tot = logsfr_sfms + logsfr_quench
        qing_SFRs[i_z] = logsfr_tot

        logsfr_quench_sat = sfr_evol.DeltaLogSFR_quenching(
                tqq, 
                t_cosmics[i_z],
                M_q=avg_Msham_evol[closest_tQ_index],
                tau_prop={'name': 'satellite'})
        logsfr_tot_sat = logsfr_sfms + logsfr_quench_sat
        qing_SFRs_sat[i_z] = logsfr_tot_sat
    '''
    # overquenching 
    avg_q_ssfr = sfr_evol.AverageLogSSFR_q_peak(avg_Msham_evol[-1])
    sigma_q_ssfr = sfr_evol.ScatterLogSSFR_q_peak(avg_Msham_evol[-1])
    min_q_ssfr = (sigma_q_ssfr + avg_q_ssfr)#[0]
    qing_ssfr = qing_SFRs - avg_Msham_evol 
    qing_SFRs[np.where(qing_ssfr < min_q_ssfr)] = \
            avg_Msham_evol[np.where(qing_ssfr < min_q_ssfr)] + min_q_ssfr
    '''
    # figure 
    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(7, 6))
    sub = fig.add_subplot(111) 

    sub.plot(t_cosmics, qing_SFRs_sat, color=pretty_colors[7], lw=2, ls='--')
    sub.plot(t_cosmics, qing_SFRs, color=pretty_colors[6], lw=5)
    #sub.text(9.8, -.7, 'Quenching ($\mathtt{t_{Q, start}=9}$ Gyr)', rotation=-47.5) 
    sub.text(10.5, -1.025, 'Central', rotation=-47.5) 
    sub.text(9.5, -1.05, 'Satellite', rotation=-60.) 
    sub.plot(t_cosmics, sf_SFRs, color=pretty_colors[1], lw=5, ls='--')
    sub.text(8.5, -0.24, 'Star Forming', rotation=-4, fontsize=18) 

    avg_q_ssfr = sfr_evol.AverageLogSSFR_q_peak(avg_Msham_evol[-1])
    sigma_q_ssfr = sfr_evol.ScatterLogSSFR_q_peak(avg_Msham_evol[-1])
    sub.fill_between(t_cosmics, 
            np.repeat(avg_q_ssfr + sigma_q_ssfr + avg_Msham_evol[-1], len(t_cosmics)), 
            np.repeat(-3., len(t_cosmics)), color=pretty_colors[4]) 
    #print np.repeat(avg_q_ssfr + sigma_q_ssfr + avg_Msham_evol[-1], len(t_cosmics))

    sub.text(8.7, -2.18, 'Quiescent', fontsize=22) 
    sub2 = sub.twiny()
    cosmo = astrocosmo.FlatLambdaCDM(H0=70, Om0=0.274)
    reds = np.array([0.9, 0.7, 0.5, 0.3, 0.1])
    ages = cosmo.age(reds).value
    sub2.set_xticks(ages)
    sub2.set_xticklabels(['{:g}'.format(age) for age in reds])
    sub2.set_xlabel('Redshift', fontsize=25) 
    sub2.tick_params(axis='x', pad=0)
    
    sub.set_xlim([t_cosmics.min(), t_cosmics.max()])
    sub2.set_xlim([t_cosmics.min(), t_cosmics.max()])
    sub.set_xlabel(r'$\mathtt{t_{cosmic}}\;[\mathtt{Gyr}]$', fontsize=25) 
    sub.set_ylim([-2.4, 0.5]) 
    sub.set_ylabel(r'$\mathtt{log}(\mathtt{SFR}\;[\mathtt{M}_\odot/\mathtt{yr}])$', fontsize=25) 
    fig_file = ''.join(['figure/paper/', 'SFH_SchematicDemo.png']) 
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    Util.png2pdf(fig_file)
    plt.close() 
    return None 
    

def fig_SSFRevol(t, abcrun, prior_name='try0', orientation='portrait'): 
    ''' Demonstrate the SFR assignment scheme through the 
    SSFR distribution function P(SSFR) evolution over multiple
    snapshots. Also include extreme versions of the quenching
    timescale for comparison. 
    '''
    abc_plot = PlotABC(t, abcrun=abcrun, prior_name=prior_name) 
    # median theta values 
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = abc_plot.med_theta
    print '###### abc median loaded'

    # other parameters
    sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
    sim_kwargs = sfinherit_kwargs.copy()
    sim_kwargs['evol_prop']['fudge'] = {
            'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
    sim_kwargs['evol_prop']['tau'] = {
            'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
    sim_kwargs['sfr_prop']['gv'] = {
            'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}

    nsnaps = [1,4,7,10,12]
    color_scheme = ['#6c1500', '#b42400', '#c34f32', '#d27b66', '#e1a799'][::-1]
    #['#cc2800', '#991e00', '#661400', '#330a00', '#000000']

    print '###### about to run inherit'
    #nsnaps = [1,4,10]
    inh = Inherit(nsnaps, 
            nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
            subhalo_prop=sim_kwargs['subhalo_prop'], 
            sfr_prop=sim_kwargs['sfr_prop'], 
            evol_prop=sim_kwargs['evol_prop'])
    print '###### finished run inherit'
    anc = inh.ancestor
    des_dict = inh()

    prettyplot() 
    pretty_colors = prettycolors()
    if orientation == 'portrait': 
        fig = plt.figure(1, figsize=(6,10))
        sub = fig.add_subplot(211)
        sub2 = fig.add_subplot(212)
    elif orientation == 'landscape': 
        fig = plt.figure(1, figsize=(14,7))
        sub = fig.add_subplot(121)
        sub2 = fig.add_subplot(122)
    bkgd = fig.add_subplot(111, frameon=False)

    bin_mid, ssfr_dist = anc.Ssfr()
    sub.plot(bin_mid[1], ssfr_dist[1], 
            c='#f0d3cc', lw=2, ls='-', label="z = 1.")
    
    for ii, i_snap in enumerate(nsnaps[::-1]): 
        bin_mid, ssfr_dist = des_dict[str(i_snap)].Ssfr()
        label = 'z = '+str(round(des_dict[str(i_snap)].zsnap,1))
        
        if i_snap == 1: 
            lwidth = 3
            lbl = label 
        elif i_snap == 7: 
            lwidth = 2
            lbl = label 
        else: 
            lwidth = 2
            lbl = None
        sub.plot(bin_mid[1], ssfr_dist[1], 
                c=color_scheme[ii], lw=lwidth, ls='-', label=lbl)
        if i_snap == 1: 
            sub2.plot(bin_mid[1], ssfr_dist[1], 
                    c='k', lw=3, ls='-')
        #sub.text(-10.75, 1.25, r"$\mathtt{log\;M_* = [10.1, 10.5]}$", fontsize=20)

    qfrac = Fq() 
    sfq = qfrac.Classify(des_dict['1'].mass, des_dict['1'].sfr, des_dict['1'].zsnap, 
            sfms_prop=sim_kwargs['sfr_prop']['sfms'])
    ngal = len(np.where((des_dict['1'].mass >= 10.1) & (des_dict['1'].mass < 10.5))[0]) 
    ngal_q = len(np.where(
        (des_dict['1'].mass >= 10.1) & 
        (des_dict['1'].mass < 10.5) & 
        (sfq == 'quiescent'))[0]) 
    print 'Standard fq = ', np.float(ngal_q)/np.float(ngal)
    
    # plot extreme quenching timescale versions 
    sat_ppp = PlotABC(10, abcrun='SatABC_TinkerFq', prior_name='satellite')
    sat_gv_slope, sat_gv_offset, sat_fudge_slope, sat_fudge_offset = sat_ppp.med_theta
    fast_kwargs = sim_kwargs.copy()     # faster quenching
    fast_kwargs['sfr_prop']['gv'] = {'slope': sat_gv_slope, 'fidmass': 10.5, 'offset': sat_gv_offset}
    fast_kwargs['evol_prop']['fudge'] = {'slope': sat_fudge_slope, 'fidmass': 10.5, 'offset': sat_fudge_offset}
    fast_kwargs['evol_prop']['tau'] = {'name': 'satellite'}
    fast_inh = Inherit([1], 
            nsnap_ancestor=fast_kwargs['nsnap_ancestor'],
            subhalo_prop=fast_kwargs['subhalo_prop'], 
            sfr_prop=fast_kwargs['sfr_prop'], 
            evol_prop=fast_kwargs['evol_prop'])
    fast_des_dict = fast_inh()
    
    fast_bin_mid, fast_ssfr_dist = fast_des_dict['1'].Ssfr()
    sub2.plot(fast_bin_mid[1], fast_ssfr_dist[1], 
            c='k', lw=3, ls='--', label=r"Shorter $\tau_Q$" )
    
    fast_sfq = qfrac.Classify(fast_des_dict['1'].mass, fast_des_dict['1'].sfr, fast_des_dict['1'].zsnap, 
            sfms_prop=fast_kwargs['sfr_prop']['sfms'])
    fast_ngal = len(np.where((fast_des_dict['1'].mass >= 10.1) & (fast_des_dict['1'].mass < 10.5))[0]) 
    fast_ngal_q = len(np.where(
        (fast_des_dict['1'].mass >= 10.1) & 
        (fast_des_dict['1'].mass < 10.5) & 
        (fast_sfq == 'quiescent'))[0]) 
    print 'Fast fq = ', np.float(fast_ngal_q)/np.float(fast_ngal)
    
    #longtau_ppp = PlotABC(13, abcrun='FixedLongTau_TinkerFq', prior_name='longtau')
    #long_gv_slope, long_gv_offset, long_fudge_slope, long_fudge_offset = longtau_ppp.med_theta
    slow_kwargs = sim_kwargs.copy()     # slower quenching
    slow_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
    slow_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
    slow_kwargs['evol_prop']['tau'] = {
            'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset+1.}
    slow_inh = Inherit([1], 
            nsnap_ancestor=slow_kwargs['nsnap_ancestor'],
            subhalo_prop=slow_kwargs['subhalo_prop'], 
            sfr_prop=slow_kwargs['sfr_prop'], 
            evol_prop=slow_kwargs['evol_prop'])
    slow_des_dict = slow_inh()
    
    slow_bin_mid, slow_ssfr_dist = slow_des_dict['1'].Ssfr()
    sub2.plot(slow_bin_mid[1], slow_ssfr_dist[1], 
            c='k', lw=3, ls=':', label=r"Longer $\tau_Q$" )
    
    slow_sfq = qfrac.Classify(slow_des_dict['1'].mass, slow_des_dict['1'].sfr, slow_des_dict['1'].zsnap, 
            sfms_prop=slow_kwargs['sfr_prop']['sfms'])
    slow_ngal = len(np.where((slow_des_dict['1'].mass >= 10.1) & (slow_des_dict['1'].mass < 10.5))[0]) 
    slow_ngal_q = len(np.where(
        (slow_des_dict['1'].mass >= 10.1) & 
        (slow_des_dict['1'].mass < 10.5) & 
        (slow_sfq == 'quiescent'))[0]) 
    print 'Slow fq = ', np.float(slow_ngal_q)/np.float(slow_ngal)
    
    sub.set_xlim([-13.0, -8.5])
    if orientation == 'portrait': 
        sub.set_xticks([-13, -12, -11, -10, -9]) 
        sub.set_xticklabels([]) 
    sub.set_ylim([0.0, 1.4])
    sub.set_yticks([0.0, 0.4, 0.8, 1.2]) 
    sub.minorticks_on()
    sub.legend(loc='upper left', prop={'size': 20}, 
            borderaxespad=1., handletextpad=0.0)

    sub2.set_xlim([-13.0, -8.5])
    sub2.set_ylim([0.0, 1.0])
    #sub2.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
    if orientation == 'portrait': 
        #sub2.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8])
        sub2.set_xticks([-13, -12, -11, -10, -9]) 
        sub2.set_yticks([0., 0.2, 0.4, 0.6, 0.8]) 
        sub2.minorticks_on() 
    elif orientation == 'landscape': 
        sub.set_xticklabels([-13, -12, -11, -10, -9, -8])
        #sub2.set_xticklabels([-13, '', -12, '', -11, '', -10, '', -9]) 
        sub2.set_xticklabels(['', -12.5, '', -11.5, '', -10.5, '', -9.5]) 
        sub2.yaxis.tick_right()
        sub2.yaxis.set_ticks_position('both')
        sub2.yaxis.set_label_position('right')
    sub2.legend(bbox_to_anchor=(1.05, 1.05),
            loc='upper right', borderaxespad=1.5, prop={'size': 20})
    
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
    bkgd.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
    
    if orientation == 'portrait': 
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
    elif orientation == 'landscape': 
        fig.subplots_adjust(wspace=0.1, hspace=0.0)
    fig_file = ''.join(['figure/paper/', 'SSFRevol_test', 
        '.', orientation, '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    plt.close() 
    Util.png2pdf(fig_file)
    return None 


def fig_ABC_posterior(tf, abcrun=None, prior_name='try0'):
    ''' The posterior corner plots from ABC given time step, abcrun name and prior name. 
    '''
    ppp = PlotABC(tf, abcrun=abcrun, prior_name=prior_name)

    fig_file = ''.join(['figure/paper/',
        'ABC_posterior',
        '.', abcrun, 
        '.', prior_name, '_prior', 
        '.png'])
    fig_file = ppp.Corner(filename=fig_file)
    plt.close() 

    Util.png2pdf(fig_file)
    return None 


def fig_SSFR_ABC_post(tf, abcrun=None, prior_name='try0'): 
    ''' The SSFR distribution from the median value of the ABC posterior
    '''
    figdata_file = ''.join(['/data1/hahn/paper/',  
        'SSFR.ABC_posterior', '.', abcrun, '.', prior_name, '_prior', '.p'])
    if not os.path.isfile(figdata_file):
        # model 
        ppp = PlotABC(tf, abcrun=abcrun, prior_name=prior_name)
        gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = ppp.med_theta

        sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
        sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
        sim_kwargs['evol_prop']['tau'] = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
        inh = Inherit([1], 
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        des_dict = inh() 
        descendant = des_dict['1'] 
        model_bin_mid, model_ssfr_dist = descendant.Ssfr()
        pickle.dump([model_bin_mid, model_ssfr_dist], open(figdata_file, 'wb'))
    else: 
        model_bin_mid, model_ssfr_dist = pickle.load(open(figdata_file, 'rb')) 

    # group catalog
    groupcat = GroupCat(Mrcut=18, position='central')
    sdss_bin_mid, sdss_ssfr_dist = groupcat.Ssfr()

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(20, 5))
    bkgd = fig.add_subplot(111, frameon=False)

    panel_mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    for i_m, mass_bin in enumerate(panel_mass_bins): 
        sub = fig.add_subplot(1, 4, i_m+1)

        if i_m == 3:  
            model_label = 'Hahn+(2016)'
            sdss_label = 'SDSS Centrals'
        else: 
            sdss_label = None 
            model_label = None

        sub.plot(model_bin_mid[i_m], model_ssfr_dist[i_m], 
                lw=3, ls='-', c=pretty_colors[3], label=model_label)

        sub.plot(sdss_bin_mid[i_m], sdss_ssfr_dist[i_m], 
                lw=2, ls='--', c='k', label=sdss_label)

        massbin_str = ''.join([ 
            r'$\mathtt{log \; M_{*} = [', 
            str(mass_bin[0]), ',\;', 
            str(mass_bin[1]), ']}$'
            ])
        sub.text(-12., 1.4, massbin_str, fontsize=20)
    
        # x-axis
        if i_m == 3:
            sub.set_xticks([-13, -12, -11, -10, -9])
        else: 
            sub.set_xticks([-13, -12, -11, -10])
        sub.set_xlim([-13., -9.])
        # y-axis 
        sub.set_ylim([0.0, 1.7])
        sub.set_yticks([0.0, 0.5, 1.0, 1.5])
        if i_m == 0: 
            sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
        else: 
            sub.set_yticklabels([])
        
        ax = plt.gca()
        leg = sub.legend(bbox_to_anchor=(-8.5, 1.55), loc='upper right', prop={'size': 20}, borderpad=2, 
                bbox_transform=ax.transData, handletextpad=0.5)
    
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
        
    fig.subplots_adjust(wspace=0.0, hspace=0.0)

    fig_file = ''.join(['figure/paper/',
        'SSFR.ABC_posterior',
        '.', abcrun, 
        '.', prior_name, '_prior', 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)

    #ppp.Ssfr(filename=fig_file)
    plt.close() 
    Util.png2pdf(fig_file) 
    return None


def fig_SSFR_tau_satellite(tf, abcrun='rhofq_tausat', prior_name='satellite'): 
    ''' The SSFR distribution from the median of the ABC posterior when using the 
    satellite quenching timescale model. 
    '''
    figdata_file = ''.join(['/data1/hahn/paper/',  
        'SSFR.ABC_posterior', '.', abcrun, '.', prior_name, '_prior', '.p'])
    if not os.path.isfile(figdata_file):
        # model 
        ppp = PlotABC(tf, abcrun=abcrun, prior_name=prior_name)
        gv_slope, gv_offset, fudge_slope, fudge_offset = ppp.med_theta

        sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
        sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
        sim_kwargs['evol_prop']['tau'] = {'name': 'satellite'}
        inh = Inherit([1], 
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        des_dict = inh() 
        descendant = des_dict['1'] 
        model_bin_mid, model_ssfr_dist = descendant.Ssfr()
        pickle.dump([model_bin_mid, model_ssfr_dist], open(figdata_file, 'wb'))
    else: 
        model_bin_mid, model_ssfr_dist = pickle.load(open(figdata_file, 'rb')) 

    # group catalog
    groupcat = GroupCat(Mrcut=18, position='central')
    sdss_bin_mid, sdss_ssfr_dist = groupcat.Ssfr()

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(20, 5))
    bkgd = fig.add_subplot(111, frameon=False)

    panel_mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    for i_m, mass_bin in enumerate(panel_mass_bins): 
        sub = fig.add_subplot(1, 4, i_m+1)

        if i_m == 3:  
            model_label = r'$\tau_\mathtt{Q}^\mathtt{sat}$'
            sdss_label = 'SDSS Centrals'
        else: 
            sdss_label = None 
            model_label = None

        sub.plot(model_bin_mid[i_m], model_ssfr_dist[i_m], 
                lw=3, ls='-', c=pretty_colors[7], label=model_label)

        sub.plot(sdss_bin_mid[i_m], sdss_ssfr_dist[i_m], 
                lw=2, ls='--', c='k', label=sdss_label)

        massbin_str = ''.join([ 
            r'$\mathtt{log \; M_{*} = [', 
            str(mass_bin[0]), ',\;', 
            str(mass_bin[1]), ']}$'
            ])
        sub.text(-12., 1.4, massbin_str, fontsize=20)
    
        # x-axis
        if i_m == 3:
            sub.set_xticks([-13, -12, -11, -10, -9])
        else: 
            sub.set_xticks([-13, -12, -11, -10])
        sub.set_xlim([-13., -9.])
        #sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
        # y-axis 
        sub.set_ylim([0.0, 1.7])
        sub.set_yticks([0.0, 0.5, 1.0, 1.5])
        if i_m == 0: 
            sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
        else: 
            sub.set_yticklabels([])
        
        ax = plt.gca()
        leg = sub.legend(bbox_to_anchor=(-8.5, 1.55), loc='upper right', prop={'size': 20}, borderpad=2, 
                bbox_transform=ax.transData, handletextpad=0.5)
    
        bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
        
    fig.subplots_adjust(wspace=0.0, hspace=0.0)

    fig_file = ''.join(['figure/paper/',
        'SSFR.ABC_posterior',
        '.', abcrun, 
        '.', prior_name, '_prior', 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)

    #ppp.Ssfr(filename=fig_file)
    plt.close() 
    Util.png2pdf(fig_file) 
    return None


def fig_SSFR_SDSS(): 
    ''' Compare P(SSFR) between the SDSS centrals and satellites.
    '''
    # group catalog CENTRALS
    groupcat = GroupCat(Mrcut=18, position='central')
    cen_bin_mid, cen_ssfr_dist = groupcat.Ssfr()
    
    # group catalog SATELLITES
    groupcat = GroupCat(Mrcut=18, position='satellite')
    sat_bin_mid, sat_ssfr_dist = groupcat.Ssfr()

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(6, 6))

    panel_mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    for i_m, mass_bin in enumerate(panel_mass_bins): 
        if i_m != 1: 
            continue 
        sub = fig.add_subplot(111)
        cen_label = 'SDSS Centrals'
        sat_label = 'SDSS Satellites'

        sub.plot(sat_bin_mid[i_m], sat_ssfr_dist[i_m], 
                lw=3, ls='--', c=pretty_colors[3], label=sat_label)
        
        sub.plot(cen_bin_mid[i_m], cen_ssfr_dist[i_m], 
                lw=3, ls='-', c=pretty_colors[1], label=cen_label)

        #sub.fill_between(np.arange(-11., -9., 0.1), 
        #        np.repeat(0., len(np.arange(-11., -9., 0.1))), 
        #        np.repeat(2., len(np.arange(-11., -9., 0.1))), lw=0, color='k', alpha=0.25) 
        #sub.fill_between(np.arange(-13., -11.75, 0.1), 
        #        np.repeat(0., len(np.arange(-13., -11.75, 0.1))), 
        #        np.repeat(2., len(np.arange(-13., -11.75, 0.1))), lw=0, color='k', alpha=0.25) 
    
        mmbin = np.where((cen_bin_mid[i_m] >= -11.75) & (cen_bin_mid[i_m] <= -11.)) 
        print np.max([cen_ssfr_dist[i_m][mmbin], sat_ssfr_dist[i_m][mmbin]], axis=0)
        sub.fill_between(cen_bin_mid[i_m][mmbin], 
                np.repeat(0., len(mmbin[0])), 
                np.max([cen_ssfr_dist[i_m][mmbin], sat_ssfr_dist[i_m][mmbin]], axis=0),
                lw=0, color='g', alpha=0.25) 

        #sub.fill_between(np.arange(-11.75, -11., 0.1), 
        #        np.repeat(0., len(np.arange(-11.75, -11., 0.1))), 
        #        np.repeat(2., len(np.arange(-11.75, -11., 0.1))), lw=0, color='g', alpha=0.25) 

        massbin_str = ''.join([ 
            r'$\mathtt{log \; M_{*} = [', 
            str(mass_bin[0]), ',\;', 
            str(mass_bin[1]), ']}$'
            ])
        sub.text(-11.65, 0.925, massbin_str, fontsize=20)
    
        # x-axis
        sub.set_xlim([-13., -9.4])
        sub.set_xticks([-13, -12, -11, -10])
        sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
        # y-axis 
        sub.set_ylim([0.0, 1.1])
        sub.set_yticks([0.0, 0.25, 0.5, 0.75, 1.])
        sub.minorticks_on()
        sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
        
        ax = plt.gca()
        leg = sub.legend(bbox_to_anchor=(-9., 1.025), loc='upper right', prop={'size': 20}, borderpad=2, 
                bbox_transform=ax.transData, handletextpad=0.0)
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    
    fig_file = ''.join(['figure/paper/', 'SSFR_SDSS.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close() 
    Util.png2pdf(fig_file) 
    return None


def fig_SSFR_Poster(): 
    ''' The SSFR distribution from the median value of the ABC posterior
    '''
    figdata_file = ''.join(['/data1/hahn/paper/',  
        'SSFR.ABC_posterior', '.RHOssfrfq_TinkerFq_Std.updated_prior', '.p'])
    model_bin_mid, model_ssfr_dist = pickle.load(open(figdata_file, 'rb')) 
   
    sat_figdata_file = ''.join(['/data1/hahn/paper/',  
        'SSFR.ABC_posterior', '.rhofq_tausat.satellite_prior', '.p'])
    sat_model_bin_mid, sat_model_ssfr_dist = pickle.load(open(sat_figdata_file, 'rb')) 

    # group catalog
    groupcat = GroupCat(Mrcut=18, position='central')
    sdss_bin_mid, sdss_ssfr_dist = groupcat.Ssfr()

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(20, 6))

    panel_mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    for i_m, mass_bin in enumerate(panel_mass_bins): 
        sub = fig.add_subplot(1, 4, i_m+1)

        if i_m == 3:  
            model_label = 'Hahn+(2016)'
            sat_model_label = r'w/ Satellite $\tau_Q$'
            sdss_label = 'SDSS Centrals'
        else: 
            sdss_label = None 
            sat_model_label = None 
            model_label = None
        
        scale_down = model_ssfr_dist[i_m].max()/sat_model_ssfr_dist[i_m].max()
        print model_ssfr_dist[i_m].max(), sat_model_ssfr_dist[i_m].max()
        print scale_down
        sub.plot(sat_model_bin_mid[i_m], scale_down*sat_model_ssfr_dist[i_m], 
                lw=2, ls=':', c=pretty_colors[1], label=sat_model_label)

        sub.plot(model_bin_mid[i_m], model_ssfr_dist[i_m], 
                lw=3, ls='-', c=pretty_colors[3], label=model_label)

        sub.plot(sdss_bin_mid[i_m], sdss_ssfr_dist[i_m], 
                lw=2, ls='--', c='k', label=sdss_label)

        massbin_str = ''.join([ 
            r'$\mathtt{log \; M_{*} = [', 
            str(mass_bin[0]), ',\;', 
            str(mass_bin[1]), ']}$'
            ])
        sub.text(-12., 1.6, massbin_str, fontsize=20)
    
        # x-axis
        if i_m == 3:
            sub.set_xticks([-13, -12, -11, -10, -9])
        else: 
            sub.set_xticks([-13, -12, -11, -10])
        sub.set_xlim([-13., -9.])
        sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
        # y-axis 
        sub.set_ylim([0.0, 1.8])
        sub.set_yticks([0.0, 0.5, 1.0, 1.5])
        if i_m == 0: 
            sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
        else: 
            sub.set_yticklabels([])
        
        ax = plt.gca()
        leg = sub.legend(bbox_to_anchor=(-8.5, 1.75), loc='upper right', prop={'size': 20}, borderpad=2, 
                bbox_transform=ax.transData, handletextpad=0.5)
        
    fig.subplots_adjust(wspace=0.0, hspace=0.0)

    fig_file = ''.join(['figure/paper/',
        'SSFR.Poster.ABC_posterior',
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)

    #ppp.Ssfr(filename=fig_file)
    plt.close() 
    Util.png2pdf(fig_file) 
    return None


def fig_tau_ABC_post(tf, abcrun=None, prior_name='try0'): 
    ''' The tau_Q^cen for the median values of the quenching timescale parameters in the ABC posterior
    '''
    # model 
    ppp = PlotABC(tf, abcrun=abcrun, prior_name=prior_name)
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = ppp.med_theta

    med_tau_dict = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}


    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(7,7))
    sub = fig.add_subplot(111)

    m_arr = np.arange(8.5, 12.5, 0.25)
    tauq_list = [] 
    for ii in range(len(ppp.theta)): 
        tau_dict_i = {'name': 'line', 'slope': (ppp.theta[ii])[-2], 'fid_mass': 11.1, 'yint': (ppp.theta[ii])[-1]}
        sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i), 
                c=pretty_colors[8], lw=0.4, alpha=0.1)
        tauq_list.append(sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i))
    tauq_list = np.array(tauq_list)
    #a, b, c, d, e = np.percentile(tauq_list, [2.5, 16, 50, 84, 97.5], axis=0)
    #b, c, d = np.percentile(tauq_list, [16, 50, 84], axis=0)
    yerr = np.std(tauq_list, axis=0)
    
    #for ii in np.random.choice(range(len(ppp.theta)), 100):
    #    tau_dict_i = {'name': 'line', 'slope': (ppp.theta[ii])[-2], 'fid_mass': 11.1, 'yint': (ppp.theta[ii])[-1]}
    #    sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i), 
    #            c=pretty_colors[8], lw=0.5, alpha=0.2)

    sub.errorbar(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=med_tau_dict), 
            yerr=yerr, c=pretty_colors[7], fmt='o', lw=2, label='Centrals (Hahn+2016)')

    #sub.errorbar(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=med_tau_dict), 
    #        c='k', lw=3, label='Centrals (Hahn+2016)')
    m_arr = np.arange(8.5, 12.5, 0.01)
    #satplot, = sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite'}), 
    #        c='k', ls='--', lw=3, label='Satellites (Wetzel+2013)')
    satplot = sub.fill_between(10**m_arr, 
            sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite_lower'}), 
            sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite_upper'}), 
            color='none', linewidth=1, edgecolor='k', hatch='X',
            label='Satellites (Wetzel+2013)')

    # data thiefed error bar values 
    #sat_m_arr = np.array([6.34372e+9, 1.01225e+10, 1.65650e+10, 2.59064e+10, 3.93062e+10, 1.58160e+11])
    #sat_tau_top = np.array([9.50877e-1, 8.45614e-1, 7.82456e-1, 6.84211e-1, 5.29825e-1, 3.75439e-1])
    #sat_tau_bot = np.array([6.56140e-1, 5.43860e-1, 4.94737e-1, 3.82456e-1, 1.29825e-1, -2.45614e-2])
    #satplot = sub.fill_between(sat_m_arr, sat_tau_top, sat_tau_bot, 
    #        color='none', linewidth=1, edgecolor='k', hatch='X',
    #        label='Satellites (Wetzel+2013)')
    sub.set_xlabel(r"$\mathtt{log(M_*\;[M_\odot])}$", fontsize=25)
    sub.set_xscale('log') 
    sub.set_xlim([10**9.5, 10**11.5])

    sub.set_ylim([0.0, 1.7])
    sub.set_yticks([0.0, 0.5, 1.0, 1.5])
    sub.minorticks_on()
    sub.set_ylabel(r"$\tau_\mathtt{Q}\;[\mathtt{Gyr}]$", fontsize=25)
    # get handles
    handles, labels = sub.get_legend_handles_labels()
    # remove the errorbars
    for i_h, h in enumerate(handles): 
        try:
            handles[i_h] = h[0]
        except TypeError: 
            pass
    sub.legend(handles, labels, loc='upper right', numpoints=1, prop={'size': 20}, handletextpad=0.5, markerscale=3)
    
    fig_file = ''.join(['figure/paper/',
        'tau.ABC_posterior',
        '.', abcrun, 
        '.', prior_name, '_prior', 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()
    Util.png2pdf(fig_file)
    return None 


def fig_tau_SMFevol(standard_run=None, standard_tf=7, noSMF_run=None, noSMF_tf=7, extraSMF_run=None, extraSMF_tf=7): 
    ''' tau_Q^cen comparison for different SMF evolution prescription 
    '''
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(7,7))
    sub = fig.add_subplot(111)
    m_arr = np.arange(8.5, 12.1, 0.1)   # log M* 

    # Standard model 
    std = PlotABC(standard_tf, abcrun=standard_run, prior_name='updated')
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = std.med_theta
    std_med_tau_dict = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

    tauq_list = [] 
    for ii in range(len(std.theta)): 
        tau_dict_i = {
                'name': 'line', 
                'slope': (std.theta[ii])[-2], 
                'fid_mass': 11.1, 
                'yint': (std.theta[ii])[-1]
                }
        #sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i), 
        #        c=pretty_colors[8], lw=0.4, alpha=0.1)
        tauq_list.append(sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i))
    tauq_list = np.array(tauq_list)
    #a, b, c, d, e = np.percentile(tauq_list, [2.5, 16, 50, 84, 97.5], axis=0)
    b, c, d = np.percentile(tauq_list, [16, 50, 84], axis=0)
    #yerr = np.std(tauq_list, axis=0)

    #sub.errorbar(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=std_med_tau_dict), 
    #        yerr=yerr, c=pretty_colors[3], fmt='o', lw=2, label='Centrals (Hahn+2016)')
    sub.fill_between(10**m_arr, b, d, color=pretty_colors[3], alpha=0.5, edgecolor='none', label='Centrals') 

    #sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=std_med_tau_dict), 
    #        c=pretty_colors[3], lw=2, ls='--')#, label='Centrals (Hahn+2016)')

    # Extra SMF evolution  
    extrasmf = PlotABC(extraSMF_tf, abcrun=extraSMF_run, prior_name='updated')
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = extrasmf.med_theta
    extrasmf_med_tau_dict = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
    sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=extrasmf_med_tau_dict), 
            c=pretty_colors[7], lw=2, ls='--', label='Extreme SMF Evo.')

    # No SMF evolution  
    nosmf = PlotABC(noSMF_tf, abcrun=noSMF_run, prior_name='updated')
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = nosmf.med_theta
    nosmf_med_tau_dict = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
    #sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=nosmf_med_tau_dict), 
    #        c=pretty_colors[5], lw=2, label='No SMF Evo.')
    sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=nosmf_med_tau_dict), 
            c=pretty_colors[5], lw=2, ls='--', label='No SMF Evo.')
    
    m_arr = np.arange(8.5, 12.5, 0.01)
    sub.fill_between(10**m_arr, 
            sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite_lower'}), 
            sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite_upper'}), 
            color='none', linewidth=1, edgecolor='k', hatch='X',
            label='Satellites')

    #sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite'}), 
    #        c='k', ls='--', lw=3, label='Satellites')
    sub.set_xlabel(r"$\mathtt{log(M_*\;[M_\odot])}$", fontsize=25)
    sub.set_xscale('log') 
    sub.set_xlim([10**9.5, 10**11.5])

    sub.set_ylim([0.0, 1.7])
    sub.set_yticks([0.0, 0.5, 1.0, 1.5])
    sub.minorticks_on()
    sub.set_ylabel(r"$\tau_\mathtt{Q}\;[\mathtt{Gyr}]$", fontsize=25)
    sub.legend(loc='upper right', scatterpoints=1, prop={'size': 20}, handletextpad=0.5, markerscale=3)
    
    fig_file = ''.join(['figure/paper/',
        'tau.SMFevolcomparison.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()
    Util.png2pdf(fig_file)
    return None 


def fig_gas_depletion():
    ''' Explore the gas depletion time to the quenching time 
    using the gas mass scaling relation of Stewart et al. (2009)
    '''

    m_arr = np.arange(9.0, 12.1, 0.1)   # mass array 
    #print np.log10(f_stargas(m_arr, 1.0)/f_stargas(m_arr, 0.)) 
    
    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(15, 6))
    for iz, z in enumerate([0.1, 0.4, 0.8]): 
        # central quenching time 
        std = PlotABC(7, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated')
        gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = std.med_theta
        std_med_tau_dict = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
        
        f_sfr = 10**(-1.*(
                sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope': 1.14}) + 2.0))
        #(sfr_evol.AverageLogSSFR_q_peak(m_arr) + m_arr)))

        t_quench = np.log(f_sfr) * sfr_evol.getTauQ(m_arr, tau_prop=std_med_tau_dict) * -1.

        # Stewart+2009
        f_stargas = 0.04 * ((10.**m_arr)/(4.5*10.**11.))**(-0.59 * (1. + z)**0.45)
        f_gas_Stewart = 1. - 1./(1. + f_stargas)
        M_gas_Stewart = np.log10(1./(1./f_gas_Stewart-1.) * 10**m_arr)
        t_gas_Stewart = 10**M_gas_Stewart / 10**sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14})/10**9

        # Santini+2014
        if z < 0.2: 
            dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.1.dat'
        elif (z >= 0.2) & (z < 0.6): 
            dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.4.dat'
        elif (z >= 0.6) & (z < 1.): 
            dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.8.dat'
        else: 
            raise ValueError
        m_dat, fgas_dat = np.loadtxt(dat_file, delimiter=',', unpack=True, usecols=[0,1]) 
        f_gas_Santini = fgas_dat 
        M_gas_Santini = np.log10(1./(1./f_gas_Santini-1.) * 10**m_dat)
        t_gas_Santini = 10**M_gas_Santini / 10**sfr_evol.AverageLogSFR_sfms(m_dat, z, sfms_prop={'name': 'linear', 'zslope':1.14})/10**9

        # Boselli+2014
        f_gas_Boselli = 1. - 1./(1. + 10**(-0.69 * m_arr + 6.63))
        M_gas_Boselli = np.log10(1./(1./f_gas_Boselli-1.) * 10**m_arr)
        t_gas_Boselli = 10**M_gas_Boselli / 10**sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14})/10**9

        sub = fig.add_subplot(1,3,iz+1)
        sub.plot(10**m_arr, t_gas_Stewart, lw=3, c='k', label='Stewart+2009')
        sub.plot(10**m_dat, t_gas_Santini, lw=3, c='k', ls='--', 
                label='Santini+2014')
        if z < 0.2: 
            sub.plot(10**m_arr, t_gas_Boselli, lw=3, c='k', ls='-.',
                    label='Boselli+2014')
        else: 
            sub.plot([0], [0], lw=3, c='k', ls='-.',
                    label='Boselli+2014')
        sub.plot(10**m_arr, t_quench, lw=3, c=pretty_colors[3], label='Quenching Time') 
        sub.text(10**9.9, 1., r'$\mathtt{z = '+str(z)+'}$', fontsize=20)

        sub.set_xlim([10**9.7, 10**11.5]) 
        sub.set_xscale('log') 
        if iz == 1: 
            sub.set_xlabel(r'Stellar Mass', fontsize=25) 
        sub.set_ylim([0., 10.]) 
        
        if iz == 0: 
            sub.set_ylabel(r'Gas Depletion Time [Gyr]', fontsize=25) 
        elif iz == 2: 
            sub.set_yticklabels([]) 
            sub.legend(loc='best') 
        else: 
            sub.set_yticklabels([]) 

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/paper/',
        't_gas_depletion.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()
    Util.png2pdf(fig_file)
    return None 


def fig_quenching_comparison(z = 0.5): 
    ''' Compare the quenching *time* to the gas depletion time estimated
    from gas fraction of galaxies from Stewart+(2009), Boselli+(2014), 
    Santini+(2014) at the central redshift bin of our simulation. 
    '''
    m_arr = np.arange(9.0, 12.1, 0.1)   # mass array 
        
    # estimate central quenching time from the timescale
    # tQ = - tauQ * ln( SFRq/SFRsf ) 
    std = PlotABC(7, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated')
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = std.med_theta
    std_med_tau_dict = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

    tauq_list = [] 
    for ii in range(len(std.theta)): 
        tau_dict_i = {
                'name': 'line', 
                'slope': (std.theta[ii])[-2], 
                'fid_mass': 11.1, 
                'yint': (std.theta[ii])[-1]
                }
        #sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i), 
        #        c=pretty_colors[8], lw=0.4, alpha=0.1)
        tauq_list.append(sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i))
    tauq_list = np.array(tauq_list)
    b, c, d = np.percentile(tauq_list, [16, 50, 84], axis=0)

    del_logSFR = sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope': 1.14}) + 2. # - sfr_evol.ScatterLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope': 1.14})
    f_sfr = 10**(-1. * del_logSFR)                      # fractional decrease of SFR 
    t_quench = -1. * np.log(f_sfr) * sfr_evol.getTauQ(m_arr, tau_prop=std_med_tau_dict)  # centrla quenching time 
    t_quench_low = -1. * np.log(f_sfr) * b  
    t_quench_high = -1. * np.log(f_sfr) * d  
    sat_t_mig = -1. * np.log(f_sfr) * sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite'})  # satellite quenching time 
    print sfr_evol.getTauQ(m_arr, tau_prop=std_med_tau_dict)/ sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite'})
    print sat_t_mig[np.where((m_arr > 10.4) & (m_arr < 10.6))]
    print t_quench - sat_t_mig
        
    # Stewart+2009
    f_stargas = 0.04 * ((10.**m_arr)/(4.5*10.**11.))**(-0.59 * (1. + z)**0.45)
    M_gas_Stewart = f_stargas * 10**m_arr #np.log10(1./(1./f_gas_Stewart-1.) * 10**m_arr)
    t_gas_Stewart = M_gas_Stewart / 10**(sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14}) + 9.) 

    # Santini+2014
    #if z < 0.2: 
    #    dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.1.dat'
    #elif (z >= 0.2) & (z < 0.6): 
    #    dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.4.dat'
    #elif (z >= 0.6) & (z < 1.): 
    #    dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.8.dat'
    #else: 
    #    raise ValueError
    #m_dat, fgas_dat = np.loadtxt(dat_file, delimiter=',', unpack=True, usecols=[0,1]) 
    #f_gas_Santini = fgas_dat 
    if z > 0.2 and z < 0.6: 
        m_dat = np.array([10., 10.5, 11.0, 11.5]) 
        alpha = np.array([-1.53, -1.53, -1.34, -1.58]) 
        beta = np.array([-0.52, -0.52, -0.53, -0.85]) 
    elif z <= 0.2: 
        m_dat = np.array([10.25, 10.5, 11.0, 11.5]) 
        alpha = np.array([-2.17, -2.17, -2.17, -1.53]) 
        beta = np.array([-1.04, -1.04, -1.04, -0.52]) 
    f_gas_Santini = 10.**(alpha + beta * (m_dat - 11.)) 
    M_gas_Santini = np.log10(1./(1./f_gas_Santini-1.) * 10**m_dat)
    t_gas_Santini = 10**M_gas_Santini / 10**(sfr_evol.AverageLogSFR_sfms(m_dat, z, sfms_prop={'name': 'linear', 'zslope':1.14}) + 9.) 

    # Boselli+2014
    f_gas_Boselli = 1. - 1./(1. + 10**(-0.69 * m_arr + 6.63))
    M_gas_Boselli = np.log10(1./(1./f_gas_Boselli-1.) * 10**m_arr)
    t_gas_Boselli = 10**M_gas_Boselli / 10**(sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14}) + 9.) 

    # Popping+(2015) SHAMed gas results 
    # import Popping data 
    popping_file = ''.join(['dat/observations/popping2015/', 'mstar.vs.SFE.dat'])
    z_pop, logm_pop, logsfe_pop, sig_logsfr_pop = np.loadtxt(popping_file, unpack=True, skiprows=4, usecols=[0,1,2,3])
    tdep_pop = (1./10**logsfe_pop)/1e9    # convert to Gyrs

    popping_m_lowz = 10**logm_pop[np.where(z_pop == 0.0)]
    popping_delt_lowz = tdep_pop[np.where(z_pop == 0.0)]
    popping_m_highz = 10**logm_pop[np.where(z_pop == 0.5)]
    popping_delt_highz = tdep_pop[np.where(z_pop == 0.5)]

    t_lowz, t_highz = 13.8099, 8.628
    popping_delt_tinterp_lowz = ((popping_delt_lowz - popping_delt_highz)/(t_lowz - t_highz)) *\
            (11.8271 - t_highz) + popping_delt_highz
    popping_delt_tinterp_highz = ((popping_delt_lowz - popping_delt_highz)/(t_lowz - t_highz))*\
            (10.5893 - t_highz) + popping_delt_highz
    popping_delt_tinterp = ((popping_delt_lowz - popping_delt_highz)/(t_lowz - t_highz)) *\
            (11.1980 - t_highz) + popping_delt_highz

    # make the figure  
    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(6, 6))
    sub = fig.add_subplot(111)
    #stewart_plot, = sub.plot(10**m_arr, t_gas_Stewart, lw=3, c='k', 
    #        label=r'$\mathtt{\hat{t}_{gas}}$ Stewart+2009')
    #boselli_plot, = sub.plot(10**m_arr, t_gas_Boselli, lw=3, c='k', ls='-.',
    #        label=r'$\mathtt{\hat{t}_{gas}}$ Boselli+2014')
    #santini_plot, = sub.plot(10**m_dat, t_gas_Santini, lw=3, c='k', ls='dotted', 
    #        label=r'$\mathtt{\hat{t}_{gas}}$ Santini+2014')
    popping_plot = sub.fill_between(popping_m_highz, 
            popping_delt_tinterp_highz, popping_delt_tinterp_lowz, 
            color=pretty_colors[1], alpha=0.5, edgecolor='none', 
            label=r'$\mathtt{t_{gas}}$ Popping+2015') 
    sub.plot(popping_m_lowz, popping_delt_tinterp, lw=2, alpha=0.75,
            c=pretty_colors[1], ls='-.')

    #peng_m = np.array([10.121919584954604, 10.348897535667964, 10.71206225680934, 10.890402075226978, 11.110894941634243])
    #peng_m = np.array([9.92088197146563, 10.488326848249027, 10.757457846952011, 11.033073929961091, 11.123865110246435])
    peng_m = np.array([9.429542645241039, 10.11248454882571, 10.328800988875155, 10.65945611866502, 10.869592088998765, 11.088998763906057]) 
    peng_delt = np.array([5., 5., 4., 3., 2., 1.])

    #peng_plot, = sub.plot(10**np.arange(9.5, 11.1, 0.1), 
    #        np.repeat(4., len(np.arange(9.5, 11.1, 0.1))), lw=3, c='k', ls='--',
    #        label=r'$\mathtt{\hat{t}_{Q}}$ Peng+2015') 
    peng_plot, = sub.plot(10**peng_m, 
            peng_delt, lw=3, c='k', ls='--',
            label=r'$\mathtt{t_{mig}}$ Peng+2015') 
    #q_plot, = sub.plot(10**m_arr, t_quench, lw=3, c=pretty_colors[3], 
    #        label=r'$\mathtt{\hat{t}_{Q}^{cen}}$ Hahn+2016') 
    #sub.plot(10**m_arr, t_quench_low, lw=3, c=pretty_colors[3])
    #sub.plot(10**m_arr, t_quench_high, lw=3, c=pretty_colors[3])
    q_plot = sub.fill_between(10**m_arr, t_quench_low, t_quench_high, color=pretty_colors[3], 
            alpha=0.5, edgecolor='none', label=r'$\mathtt{t_{mig}^{cen}}$ Hahn+2016') 
    

    # Martig+2009
    martig_plot = sub.scatter([2.*10**11], [2.5], c='k', s=75, marker='*',
            label=r'Martig+2009(Morph)') 
    # Haywood+2016
    del_logSFR = sfr_evol.AverageLogSFR_sfms(np.array([np.log10(6.*10**10)]), z, sfms_prop={'name': 'linear', 'zslope': 1.14}) + 2. 
    f_sfr = 10**(-1. * del_logSFR)                      # fractional decrease of SFR 
    t_q_haywood = -1. * np.log(f_sfr) * (1.5/np.log(10.))
    haywood_plot = sub.scatter([6.*10**10], t_q_haywood, c='k', s=75, marker='^', zorder=5,
            label=r'Haywood+2016(Morph)') 

    sub.set_xlim([10**9.7, 10**11.5]) 
    sub.set_xscale('log') 
    sub.set_xlabel(r'$\mathtt{M}_*$ [$\mathtt{M_\odot}$]', fontsize=25) 
    sub.set_ylim([0., 6.75]) 
    sub.set_ylabel(r'$\mathtt{t_{mig}}$ [Gyr]', fontsize=25) 
    first_legend = plt.legend(handles=[martig_plot, haywood_plot], 
            loc='upper right', scatterpoints=1, 
            handletextpad=0.1, markerscale=1.3, scatteryoffsets=[0.5,0.5]) 
    
    ax = plt.gca().add_artist(first_legend)
    plt.legend(handles=[q_plot, peng_plot, popping_plot], #santini_plot, 
            handletextpad=0.1, loc='lower left')

    #first_legend = plt.legend(handles=[stewart_plot, boselli_plot, peng_plot], #santini_plot, 
    #        handletextpad=-0.1, loc='upper right')
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/paper/',
        't_quenching_comparison', 
        '.z', str(z), 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()
    Util.png2pdf(fig_file)
    return None 


def fig_quenching_comp_proposal(z = 0.5): 
    ''' Compare the quenching *time* to the gas depletion time estimated
    from gas fraction of galaxies from Stewart+(2009), Boselli+(2014), 
    Santini+(2014) at the central redshift bin of our simulation. 
    '''
    m_arr = np.arange(9.0, 12.1, 0.1)   # mass array 
        
    # estimate central quenching time from the timescale
    # tQ = - tauQ * ln( SFRq/SFRsf ) 
    std = PlotABC(7, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated')
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = std.med_theta
    std_med_tau_dict = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

    tauq_list = [] 
    for ii in range(len(std.theta)): 
        tau_dict_i = {
                'name': 'line', 
                'slope': (std.theta[ii])[-2], 
                'fid_mass': 11.1, 
                'yint': (std.theta[ii])[-1]
                }
        #sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i), 
        #        c=pretty_colors[8], lw=0.4, alpha=0.1)
        tauq_list.append(sfr_evol.getTauQ(m_arr, tau_prop=tau_dict_i))
    tauq_list = np.array(tauq_list)
    b, c, d = np.percentile(tauq_list, [16, 50, 84], axis=0)

    del_logSFR = sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope': 1.14}) + 2. # - sfr_evol.ScatterLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope': 1.14})
    f_sfr = 10**(-1. * del_logSFR)                      # fractional decrease of SFR 
    t_quench = -1. * np.log(f_sfr) * sfr_evol.getTauQ(m_arr, tau_prop=std_med_tau_dict)  # centrla quenching time 
    t_quench_low = -1. * np.log(f_sfr) * b  
    t_quench_high = -1. * np.log(f_sfr) * d  
    sat_t_mig = -1. * np.log(f_sfr) * sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite'})  # satellite quenching time 
    print sfr_evol.getTauQ(m_arr, tau_prop=std_med_tau_dict)/ sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite'})
    print sat_t_mig[np.where((m_arr > 10.4) & (m_arr < 10.6))]
    print t_quench - sat_t_mig
        
    # Stewart+2009
    f_stargas = 0.04 * ((10.**m_arr)/(4.5*10.**11.))**(-0.59 * (1. + z)**0.45)
    M_gas_Stewart = f_stargas * 10**m_arr #np.log10(1./(1./f_gas_Stewart-1.) * 10**m_arr)
    t_gas_Stewart = M_gas_Stewart / 10**(sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14}) + 9.) 

    # Santini+2014
    #if z < 0.2: 
    #    dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.1.dat'
    #elif (z >= 0.2) & (z < 0.6): 
    #    dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.4.dat'
    #elif (z >= 0.6) & (z < 1.): 
    #    dat_file = Util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.8.dat'
    #else: 
    #    raise ValueError
    #m_dat, fgas_dat = np.loadtxt(dat_file, delimiter=',', unpack=True, usecols=[0,1]) 
    #f_gas_Santini = fgas_dat 
    if z > 0.2 and z < 0.6: 
        m_dat = np.array([10., 10.5, 11.0, 11.5]) 
        alpha = np.array([-1.53, -1.53, -1.34, -1.58]) 
        beta = np.array([-0.52, -0.52, -0.53, -0.85]) 
    elif z <= 0.2: 
        m_dat = np.array([10.25, 10.5, 11.0, 11.5]) 
        alpha = np.array([-2.17, -2.17, -2.17, -1.53]) 
        beta = np.array([-1.04, -1.04, -1.04, -0.52]) 
    f_gas_Santini = 10.**(alpha + beta * (m_dat - 11.)) 
    M_gas_Santini = np.log10(1./(1./f_gas_Santini-1.) * 10**m_dat)
    t_gas_Santini = 10**M_gas_Santini / 10**(sfr_evol.AverageLogSFR_sfms(m_dat, z, sfms_prop={'name': 'linear', 'zslope':1.14}) + 9.) 

    # Boselli+2014
    f_gas_Boselli = 1. - 1./(1. + 10**(-0.69 * m_arr + 6.63))
    M_gas_Boselli = np.log10(1./(1./f_gas_Boselli-1.) * 10**m_arr)
    t_gas_Boselli = 10**M_gas_Boselli / 10**(sfr_evol.AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14}) + 9.) 

    # Popping+(2015) SHAMed gas results 
    # import Popping data 
    popping_file = ''.join(['dat/observations/popping2015/', 'mstar.vs.SFE.dat'])
    z_pop, logm_pop, logsfe_pop, sig_logsfr_pop = np.loadtxt(popping_file, unpack=True, skiprows=4, usecols=[0,1,2,3])
    tdep_pop = (1./10**logsfe_pop)/1e9    # convert to Gyrs

    popping_m_lowz = 10**logm_pop[np.where(z_pop == 0.0)]
    popping_delt_lowz = tdep_pop[np.where(z_pop == 0.0)]
    popping_m_highz = 10**logm_pop[np.where(z_pop == 0.5)]
    popping_delt_highz = tdep_pop[np.where(z_pop == 0.5)]

    t_lowz, t_highz = 13.8099, 8.628
    popping_delt_tinterp_lowz = ((popping_delt_lowz - popping_delt_highz)/(t_lowz - t_highz)) *\
            (11.8271 - t_highz) + popping_delt_highz
    popping_delt_tinterp_highz = ((popping_delt_lowz - popping_delt_highz)/(t_lowz - t_highz))*\
            (10.5893 - t_highz) + popping_delt_highz
    popping_delt_tinterp = ((popping_delt_lowz - popping_delt_highz)/(t_lowz - t_highz)) *\
            (11.1980 - t_highz) + popping_delt_highz

    # make the figure  
    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(6, 4.5))
    sub = fig.add_subplot(111)
    #stewart_plot, = sub.plot(10**m_arr, t_gas_Stewart, lw=3, c='k', 
    #        label=r'$\mathtt{\hat{t}_{gas}}$ Stewart+2009')
    #boselli_plot, = sub.plot(10**m_arr, t_gas_Boselli, lw=3, c='k', ls='-.',
    #        label=r'$\mathtt{\hat{t}_{gas}}$ Boselli+2014')
    #santini_plot, = sub.plot(10**m_dat, t_gas_Santini, lw=3, c='k', ls='dotted', 
    #        label=r'$\mathtt{\hat{t}_{gas}}$ Santini+2014')
    popping_plot = sub.fill_between(popping_m_highz, 
            popping_delt_tinterp_highz, popping_delt_tinterp_lowz, 
            color=pretty_colors[1], alpha=0.5, edgecolor='none', 
            label=r'Strangulation') 
    sub.plot(popping_m_lowz, popping_delt_tinterp, lw=2, alpha=0.75,
            c=pretty_colors[1], ls='-.')

    #peng_m = np.array([10.121919584954604, 10.348897535667964, 10.71206225680934, 10.890402075226978, 11.110894941634243])
    #peng_m = np.array([9.92088197146563, 10.488326848249027, 10.757457846952011, 11.033073929961091, 11.123865110246435])
    #peng_m = np.array([9.429542645241039, 10.11248454882571, 10.328800988875155, 10.65945611866502, 10.869592088998765, 11.088998763906057]) 
    #peng_delt = np.array([5., 5., 4., 3., 2., 1.])

    #peng_plot, = sub.plot(10**np.arange(9.5, 11.1, 0.1), 
    #        np.repeat(4., len(np.arange(9.5, 11.1, 0.1))), lw=3, c='k', ls='--',
    #        label=r'$\mathtt{\hat{t}_{Q}}$ Peng+2015') 
    #peng_plot, = sub.plot(10**peng_m, 
    #        peng_delt, lw=3, c='k', ls='--',
    #        label=r'$\mathtt{t_{mig}}$ Peng+2015') 
    #q_plot, = sub.plot(10**m_arr, t_quench, lw=3, c=pretty_colors[3], 
    #        label=r'$\mathtt{\hat{t}_{Q}^{cen}}$ Hahn+2016') 
    #sub.plot(10**m_arr, t_quench_low, lw=3, c=pretty_colors[3])
    #sub.plot(10**m_arr, t_quench_high, lw=3, c=pretty_colors[3])
    q_plot = sub.fill_between(10**m_arr, t_quench_low, t_quench_high, color=pretty_colors[3], 
            alpha=0.5, edgecolor='none', label=r'Hahn+2016') 
    

    # Haywood+2016
    del_logSFR = sfr_evol.AverageLogSFR_sfms(np.array([np.log10(6.*10**10)]), z, sfms_prop={'name': 'linear', 'zslope': 1.14}) + 2. 
    f_sfr = 10**(-1. * del_logSFR)                      # fractional decrease of SFR 
    t_q_haywood = -1. * np.log(f_sfr) * (1.5/np.log(10.))
    # Martig+2009
    mq_plot = sub.scatter([6.*10**10, 2.*10**11], [t_q_haywood, 2.5], c='k', s=200, marker='*',
            label=r'Morphological')
    #martig_plot = sub.scatter([2.*10**11], [2.5], c='k', s=75, marker='*',
    #        label=r'Martig+2009(Morph)') 
    #haywood_plot = sub.scatter([6.*10**10], t_q_haywood, c='k', s=75, marker='^', zorder=5,
    #        label=r'Haywood+2016(Morph)') 

    sub.set_xlim([10**9.7, 10**11.5]) 
    sub.set_xscale('log') 
    sub.set_xlabel(r'$\mathtt{M}_*$ [$\mathtt{M_\odot}$]', fontsize=35) 
    sub.set_ylim([2., 6.75]) 
    sub.set_ylabel(r'$\mathtt{t_{quenching}}$ [Gyr]', fontsize=35) 
    first_legend = plt.legend(handles=[mq_plot], 
            loc='lower left', scatterpoints=1, 
            handletextpad=0.1, markerscale=1, scatteryoffsets=[0.5,0.5]) 
    
    ax = plt.gca().add_artist(first_legend)
    plt.legend(handles=[q_plot, popping_plot], #santini_plot, 
            handletextpad=0.1, loc='upper right')

    #first_legend = plt.legend(handles=[stewart_plot, boselli_plot, peng_plot], #santini_plot, 
    #        handletextpad=-0.1, loc='upper right')
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/paper/',
        't_quenching_comparison', 
        '.z', str(z), 
        '.proposal', '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()
    Util.png2pdf(fig_file)
    return None 



def splashback(tf, abcrun=None, prior_name=None): 
    ppp = PlotABC(tf, abcrun=abcrun, prior_name=prior_name)
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = ppp.med_theta

    sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
    sim_kwargs = sfinherit_kwargs.copy()
    sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
    sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
    sim_kwargs['evol_prop']['tau'] = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
    inh = Inherit([1], 
            nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
            subhalo_prop=sim_kwargs['subhalo_prop'], 
            sfr_prop=sim_kwargs['sfr_prop'], 
            evol_prop=sim_kwargs['evol_prop'])
    des_dict = inh() 
    des = des_dict['1'] 

    # ancestor
    anc = inh.ancestor 
    anc_snap_index = anc.snap_index
    succession, will = Util.intersection_index(
            getattr(des, 'ancestor'+str(sim_kwargs['nsnap_ancestor'])), 
            anc_snap_index
            )
    M_sham_hist = anc.Msham_evol[0][will]
    splash = np.zeros(len(succession))
    
    n_tot = 0 
    n_splash = 0 
    n_pure_splash = 0
    for ii in range(len(M_sham_hist)): 
        #if des.mass[succession[ii]] > mass_limit: 
        #    n_tot += 1 
        if np.min((M_sham_hist[ii])[:des.nsnap_genesis[succession[ii]]]) <= 0.:  
            splash[ii] = 1
            n_splash += 1
            if np.sum((M_sham_hist[ii])[:des.nsnap_genesis[succession[ii]]] <= 0.) > 1: 
                n_pure_splash += 1
    print np.float(np.sum(splash))/np.float(len(splash))
    
    model_bin_mid, model_ssfr_dist = des.Ssfr()
        
    # no splash backs
    des.data_columns = ['mass', 'ssfr'] 
    des.sample_trim(np.where(splash == 0))
    model_bin_mid, model_ssfr_dist_nosplash = des.Ssfr()

    # group catalog
    groupcat = GroupCat(Mrcut=18, position='central')
    sdss_bin_mid, sdss_ssfr_dist = groupcat.Ssfr()

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(20, 5))
    bkgd = fig.add_subplot(111, frameon=False)

    panel_mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    for i_m, mass_bin in enumerate(panel_mass_bins): 
        sub = fig.add_subplot(1, 4, i_m+1)

        if i_m == 3:  
            model_label = 'Hahn+(2016)'
            sdss_label = 'SDSS Centrals'
        else: 
            sdss_label = None 
            model_label = None

        sub.plot(model_bin_mid[i_m], model_ssfr_dist[i_m], 
                lw=3, ls='-', c=pretty_colors[3], label=model_label)

        sub.plot(model_bin_mid[i_m], model_ssfr_dist_nosplash[i_m], 
                lw=3, ls='-', c=pretty_colors[7], label='No splashback')

        sub.plot(sdss_bin_mid[i_m], sdss_ssfr_dist[i_m], 
                lw=2, ls='--', c='k', label=sdss_label)

        massbin_str = ''.join([ 
            r'$\mathtt{log \; M_{*} = [', 
            str(mass_bin[0]), ',\;', 
            str(mass_bin[1]), ']}$'
            ])
        sub.text(-12., 1.4, massbin_str, fontsize=20)
    
        # x-axis
        if i_m == 3:
            sub.set_xticks([-13, -12, -11, -10, -9])
        else: 
            sub.set_xticks([-13, -12, -11, -10])
        sub.set_xlim([-13., -9.])
        # y-axis 
        sub.set_ylim([0.0, 1.7])
        sub.set_yticks([0.0, 0.5, 1.0, 1.5])
        if i_m == 0: 
            sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
        else: 
            sub.set_yticklabels([])
        
        ax = plt.gca()
        leg = sub.legend(bbox_to_anchor=(-8.5, 1.55), loc='upper right', prop={'size': 20}, borderpad=2, 
                bbox_transform=ax.transData, handletextpad=0.5)
    
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
        
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/paper/',
        'SSFR.splashback_test.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close() 
    return None 


def non_decreasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))


def keep_non_descreasing(L): 
    out = [L[0]]
    indices = [0]
    for i_ll in range(1, len(L)):
        ll = L[i_ll]
        if ll >= out[-1]: 
            out.append(ll)
            indices.append(i_ll)
    return out, indices


if __name__=='__main__': 
    #fig_SSFR_Poster()
    #fig_SMFevol()
    #fig_fQcen_evol()
    #fig_SFMSevol()
    #for z in [0.1, 0.3, 0.5, 0.7, 0.9]: 
    #    fig_gas_depletion(z=z)
    #fig_gas_depletion()
    #fig_gas_depletion_Santini(z=0.5)
    #figSFH_SchematicDemo(7, 'RHOssfrfq_TinkerFq_Std', prior_name='updated')
    #fig_SSFR_SDSS()
    #fig_quenching_comparison(z = 0.2)
    #fig_quenching_comparison(z = 0.5)
    #fig_quenching_comp_proposal(z = 0.5)
    #fig_tau_SMFevol(
    #        standard_run='RHOssfrfq_TinkerFq_Std', standard_tf=7, 
    #        noSMF_run='RHOssfrfq_TinkerFq_NOSMFevol', noSMF_tf=8, 
    #        extraSMF_run='RHOssfrfq_TinkerFq_XtraSMF', extraSMF_tf=9)
    #splashback(7, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated')
    fig_SSFRevol(7, 'multirho_inh', prior_name='try0', orientation='portrait')
    #fig_SSFRevol(6, 'RHOssfrfq_TinkerFq_Std', prior_name='updated')
    #fig_SFRassign(7, 'RHOssfrfq_TinkerFq_Std', prior_name='updated')
    #figSFH_demo(7, 'multirho_inh', prior_name='try0')
    #fig_ABC_posterior(7, abcrun='multifq_wideprior', prior_name='updated')
    #fig_SSFR_ABC_post(7, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated')
    #fig_tau_ABC_post(7, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated')
    #fig_SSFR_tau_satellite(10, abcrun='SatABC_TinkerFq', prior_name='satellite')
