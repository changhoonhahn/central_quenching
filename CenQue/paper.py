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

import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors



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

    figdata_file = ''.join(['/data1/hahn/', 'figSFH_demo.data_file.p']) 
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
        q_Msham_i = (anc.Msham_evol[0])[i_ran]
        if q_Msham_i.min() > 9.:  
            q_SFRs_i = anc.ssfr[i_ran] + q_Msham_i
            onlyone+= 1
    # star-forming 
    onlyone = 0 
    while onlyone == 0: 
        i_ran = np.random.choice(sf_a_index1)
        sf_Msham_i = (anc.Msham_evol[0])[i_ran][::-1]
        dutycycle_prop = sfr_evol.dutycycle_param(1, 
                dutycycle_prop={'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3})
        dutycycle_prop['delta_sfr'] = 0.3

        if sf_Msham_i.min() > 9.5:  
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
                sf_SFRs_i[i_z] = logsfr_sfms + logsfr_sfduty 
            onlyone += 1
    # quenching  
    onlyone = 0 
    tqq = [7.]
    while onlyone == 0.: 
        i_ran = np.random.choice(sf_a_index1)
        qing_Msham_i = (anc.Msham_evol[0])[i_ran][::-1]
        dutycycle_prop = sfr_evol.dutycycle_param(1, 
                dutycycle_prop={'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3})
        dutycycle_prop['delta_sfr'] = -0.1

        if qing_Msham_i.min() > 9.0 and qing_Msham_i[0] < 10.:  
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

                qing_SFRs_i[i_z] = logsfr_sfms + logsfr_sfduty + logsfr_quench
            onlyone += 1

    # overquenching 
    avg_q_ssfr = sfr_evol.AverageLogSSFR_q_peak(qing_Msham_i[-1])
    sigma_q_ssfr = sfr_evol.ScatterLogSSFR_q_peak(qing_Msham_i[-1])
    min_q_ssfr = sigma_q_ssfr * np.random.randn(1) + avg_q_ssfr 
    qing_ssfr_i = qing_SFRs_i - qing_Msham_i
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
    q_SFRs_i = q_SFRs_i[::-1]
    sub1.quiver(t_cosmics[:-1], q_SFRs_i[:-1], t_cosmics[1:]-t_cosmics[:-1], q_SFRs_i[1:]-q_SFRs_i[:-1], 
            color=pretty_colors[4], scale_units='xy', angles='xy', scale=1)
    sub1.quiver(t_cosmics[:-1], sf_SFRs_i[:-1], t_cosmics[1:]-t_cosmics[:-1], sf_SFRs_i[1:]-sf_SFRs_i[:-1], 
            color=pretty_colors[1], scale_units='xy', angles='xy', scale=1)
    sub1.quiver(t_cosmics[:-1], qing_SFRs_i[:-1], t_cosmics[1:]-t_cosmics[:-1], qing_SFRs_i[1:]-qing_SFRs_i[:-1], 
            color=pretty_colors[6], scale_units='xy', angles='xy', scale=1)
    #sub1.plot(t_cosmics, q_SFRs[i,:], c=pretty_colors[4])
    
    sub1.set_xlim([5., 14.])
    sub1.set_xlabel(r'$\mathtt{t_{cosmic}}\;[\mathtt{Gyr}]$', fontsize=25) 

    sub1.set_yticklabels([])
    sub1.set_ylim([-5., 2.0]) 

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/test/', 'paper_SFH_demo.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    return None 


if __name__=='__main__': 
    figSFH_demo(7, 'multirho_inh', prior_name='try0')
