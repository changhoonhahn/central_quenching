'''

Calculations and figures pertaining to the 
central quenching timescale paper

'''
import os 
import pickle
import numpy as np 

import util.util as Util 
import sfr_evol
from inherit import Inherit
from abcee import PlotABC
from abcee import ReadABCrun

from observations import GroupCat 

import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors

def fig_SFRassign(t, abcrun, prior_name='try0'): 
    ''' Demonstrate the SFR assignment scheme based through the 
    SSFR distribution function P(SSFR) with different GV prescriptions.
    '''
    abc_plot = PlotABC(t, abcrun=abcrun, prior_name=prior_name) 
    # median theta values 
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = abc_plot.med_theta

    no_gv = [0., 0.]
    med_gv = [gv_slope, gv_offset]
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
            label = r'$\mathtt{f_{GV} = '+str(round(med_gv[0], 1))+' (log\;M_* - 10.5) + '+str(round(med_gv[1], 1))+'}$'
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


def fig_SSFRevol(t, abcrun, prior_name='try0'): 
    ''' Demonstrate the SFR assignment scheme based through the 
    SSFR distribution function P(SSFR) with different GV prescriptions.
    '''
    abc_plot = PlotABC(t, abcrun=abcrun, prior_name=prior_name) 
    # median theta values 
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = abc_plot.med_theta

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
    inh = Inherit(nsnaps, 
            nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
            subhalo_prop=sim_kwargs['subhalo_prop'], 
            sfr_prop=sim_kwargs['sfr_prop'], 
            evol_prop=sim_kwargs['evol_prop'])
    anc = inh.ancestor
    des_dict = inh()

    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(7, 7))
    sub = fig.add_subplot(111)

    bin_mid, ssfr_dist = anc.Ssfr()
    sub.plot(bin_mid[1], ssfr_dist[1], c='k', lw=2, ls='--', label=r"$\mathtt{z=1.0}$")
    
    for i_snap in nsnaps[::-1]: 
        bin_mid, ssfr_dist = des_dict[str(i_snap)].Ssfr()
        label = r"$\mathtt{\;\;"+str(round(des_dict[str(i_snap)].zsnap,1))+"}$"

        sub.plot(bin_mid[1], ssfr_dist[1], c=pretty_colors[i_snap], lw=2, ls='-', label=label)

        sub.text(-10.75, 1.4, r"$\mathtt{log\;M_* = [10.1, 10.5]}$", fontsize=20)

    sub.set_xlim([-13.0, -8.0])
    sub.set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=25) 
    sub.set_ylim([0.0, 1.6])
    sub.set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=25) 
    sub.legend(loc='upper left', prop={'size': 20})

    fig_file = ''.join(['figure/paper/', 'SSFRevol.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    plt.close() 
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
    fig = plt.figure(figsize=(20, 6))

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
    fig = plt.figure(figsize=(20, 6))

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
        'SSFR.ABC_posterior',
        '.', abcrun, 
        '.', prior_name, '_prior', 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)

    #ppp.Ssfr(filename=fig_file)
    plt.close() 
    Util.png2pdf(fig_file) 
    return None


def fig_tau_ABC_post(tf, abcrun=None, prior_name='try0'): 
    ''' The SSFR distribution from the median value of the ABC posterior
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
    satplot, = sub.plot(10**m_arr, sfr_evol.getTauQ(m_arr, tau_prop={'name': 'satellite'}), 
            c='k', ls='--', lw=3, label='Satellites (Wetzel+2014)')
    sub.set_xlabel(r"$\mathtt{log(M_*\;[M_\odot])}$", fontsize=25)
    sub.set_xscale('log') 
    sub.set_xlim([10**9.5, 10**11.5])

    sub.set_ylim([0.0, 1.7])
    sub.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
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
    #fig_SSFRevol(7, 'multirho_inh', prior_name='try0')
    #fig_SFRassign(7, 'multirho_inh', prior_name='try0')
    #figSFH_demo(7, 'multirho_inh', prior_name='try0')
    #fig_ABC_posterior(7, abcrun='multifq_wideprior', prior_name='updated')
    #fig_SSFR_ABC_post(7, abcrun='multifq_wideprior', prior_name='updated')
    #fig_tau_ABC_post(7, abcrun='multifq_wideprior', prior_name='updated')
    fig_SSFR_tau_satellite(12, abcrun='rhofq_tausat', prior_name='satellite')

