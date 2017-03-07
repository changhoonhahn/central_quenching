
import numpy as np
import time
import warnings
warnings.filterwarnings('ignore')

from gal_prop import Fq 
from gal_prop import SMF
from abcee import PlotABC
from inherit import Inherit
from lineage import Lineage
from util.util import code_dir
from util.util import get_t_nsnap
from observations import GroupCat
from central_subhalo import CentralSubhalos

from plotting.plots import PlotSMF

# --- plotting --- 
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors



def DescendantQAplot(nsnap_descendants, nsnap_ancestor=15, gv=None, tau=None, fudge=None): 
    ''' The ultimate QAplot to rule them all. 4 panels showing all the properties.
    '''
    inh_time = time.time() 
    inh = Inherit(nsnap_descendants,  nsnap_ancestor=nsnap_ancestor, gv=gv, tau=tau, fudge=fudge, quiet=True)
    des_dict = inh() 
    print time.time() - inh_time 
    
    for n_d in des_dict.keys():
        descendant = des_dict[str(n_d)]
        started_here = np.where(descendant.nsnap_genesis == nsnap_ancestor)
        start_mass = descendant.mass_genesis[started_here]

        # Mass bins
        mass_bins = np.arange(7., 12.5, 0.5)
        mass_bin_low = mass_bins[:-1]
        mass_bin_high = mass_bins[1:]

        plt.close()
        prettyplot()
        pretty_colors = prettycolors()
        fig = plt.figure(1, figsize=[25,6])
        for i_sub in range(1,5): 
            sub_i = fig.add_subplot(1,4,i_sub) 
            
            if i_sub == 1:  # SMF
                mf = SMF()
                mass, phi = mf.Obj(descendant)

                sub_i.plot(mass, phi, lw=4, c=pretty_colors[descendant.nsnap], label=r'Simulated')

                censub = CentralSubhalos()
                censub.Read(descendant.nsnap, 
                        scatter=inh.subhalo_prop['scatter'], 
                        source=inh.subhalo_prop['source'], 
                        nsnap_ancestor=nsnap_ancestor)

                mass, phi = mf._smf(censub.mass)
                sub_i.plot(mass, phi, c='k', lw=4, ls='--', label='Central Subhalos') 

                sub_i.set_ylim([10**-5, 10**-1])
                sub_i.set_xlim([7.5, 12.0])
                plt.xticks([8., 9., 10., 11., 12.])
                sub_i.set_yscale('log')

                # x,y labels
                sub_i.set_xlabel(r'Mass $\mathtt{M_*}$', fontsize=25) 
                sub_i.set_ylabel(r'Stellar Mass Function $\mathtt{\Phi}$', fontsize=25) 
                sub_i.legend(loc='upper right', frameon=False)

            elif i_sub == 2: # SFMS
                # SMF composition based on their initial mass at nsnap_ancestor 
                for i_m in range(len(mass_bin_low)): 
                    mbin = np.where((start_mass > mass_bin_low[i_m]) & (start_mass <= mass_bin_high[i_m]))

                    sub_i.scatter(
                            descendant.mass[started_here[0][mbin]], 
                            descendant.sfr[started_here[0][mbin]], 
                            color=pretty_colors[i_m], 
                            label=r"$\mathtt{M_{*,i}=}$"+str(round(mass_bin_low[i_m],2))+"-"+str(round(mass_bin_high[i_m],2))
                            )
                
                qfrac = Fq()
                m_arr = np.arange(9.0, 12.5, 0.5)
                sub_i.plot(
                        m_arr, 
                        qfrac.SFRcut(m_arr, descendant.zsnap, sfms_prop=inh.sfr_prop['sfms']), 
                        c='k', ls='--', lw=4)

                sub_i.set_xlim([9.0, 12.0])
                sub_i.set_ylim([-5.0, 2.0])
                sub_i.set_xlabel(r'$\mathtt{log\;M_*}$', fontsize=25)
                sub_i.set_ylabel(r'$\mathtt{log\;SFR}$', fontsize=25)

            elif i_sub == 3: #SMHM
                for i_m in range(len(mass_bin_low)): 
                    mbin = np.where((start_mass > mass_bin_low[i_m]) & (start_mass <= mass_bin_high[i_m]))

                    sub_i.scatter(
                            descendant.halo_mass[started_here[0][mbin]], 
                            descendant.mass[started_here[0][mbin]], 
                            color=pretty_colors[i_m], 
                            label=r"$\mathtt{M_{*,i}=}$"+str(round(mass_bin_low[i_m],2))+"-"+str(round(mass_bin_high[i_m],2))
                            )

                stellarmass = descendant.mass[started_here]
                halomass = descendant.halo_mass[started_here]
                mbin = np.arange(halomass.min(), halomass.max(), 0.25) 
                mlow = mbin[:-1]
                mhigh = mbin[1:]
                
                muMstar = np.zeros(len(mlow))
                sigMstar = np.zeros(len(mlow))

                for im in range(len(mlow)): 
                    mbin = np.where((halomass > mlow[im]) & (halomass <= mhigh[im]))
                    muMstar[im] = np.mean(stellarmass[mbin])
                    sigMstar[im] = np.std(stellarmass[mbin])
                sub_i.errorbar(0.5*(mlow+mhigh), muMstar, yerr=sigMstar, color='k', lw=3, fmt='o', capthick=2)

                sub_i.set_ylim([9.0, 12.0])
                sub_i.set_xlim([10.0, 15.0])
                sub_i.set_ylabel(r'Stellar Mass $\mathtt{M_*}$', fontsize=25) 
                sub_i.set_xlabel(r'Halo Mass $\mathtt{M_{Halo}}$', fontsize=25) 

                #sub_i.legend(loc='upper left', frameon=False, scatterpoints=1)
            elif i_sub == 4: # Fq
                #mass, fq = descendant.Fq()
                sfq = qfrac.Classify(descendant.mass, descendant.sfr, descendant.zsnap, 
                        sfms_prop=descendant.sfr_prop['sfms'])
                gc = GroupCat(Mrcut=18, position='central')
                gc.Read()
                gc_sfq = qfrac.Classify(gc.mass, gc.sfr, np.mean(gc.z), 
                        sfms_prop=descendant.sfr_prop['sfms'])
                #sub_i.plot(mass, fq, color=pretty_colors[descendant.nsnap], lw=3, ls='--', 
                #        label=r'$\mathtt{z =} '+str(descendant.zsnap)+'$')
                M_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
                M_low = M_bin[:-1]
                M_high = M_bin[1:]
                M_mid = 0.5 * (M_low + M_high)

                fq = np.zeros(len(M_low))
                gc_fq = np.zeros(len(M_low))
                for i_m in xrange(len(M_low)):
                    mlim = np.where((descendant.mass > M_low[i_m]) & (descendant.mass <= M_high[i_m]))
                    gc_mlim = np.where((gc.mass > M_low[i_m]) & (gc.mass <= M_high[i_m]))
                    ngal = np.float(len(mlim[0]))
                    gc_ngal = np.float(len(gc_mlim[0]))

                    if ngal != 0:  # no galaxy in mass bin 
                        ngal_q = np.float(len(np.where(sfq[mlim] == 'quiescent')[0]))
                        fq[i_m] = ngal_q/ngal
                    if gc_ngal != 0:
                        gc_ngal_q = np.float(len(np.where(gc_sfq[gc_mlim] == 'quiescent')[0]))
                        gc_fq[i_m] = gc_ngal_q/gc_ngal

                sub_i.plot(M_mid, fq, color=pretty_colors[descendant.nsnap], lw=3, label=r'$\mathtt{z =} '+str(descendant.zsnap)+'$')
                
                fq_model = qfrac.model(M_bin, descendant.zsnap, lit=inh.sfr_prop['fq']['name'])
                sub_i.plot(M_bin, fq_model, color='k', lw=4, ls='--', label=inh.sfr_prop['fq']['name'])
                sub_i.scatter(M_mid, gc_fq, color='k', s=100, lw=0, label='Group Catalog')

                sub_i.set_xlim([9.0, 12.0])
                sub_i.set_ylim([0.0, 1.0])

                sub_i.set_xlabel(r'Mass $\mathtt{M_*}$') 
                sub_i.set_ylabel(r'Quiescent Fraction $\mathtt{f_Q}$', fontsize=20) 

                sub_i.legend(loc='upper left', frameon=False, scatterpoints=1, markerscale=0.75)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0) 
        fig_name = ''.join([
            'figure/test/', 
            'QAplot',
            '.nsnap', str(n_d), 
            '.', ''.join([str(nd) for nd in nsnap_descendants]),
            '.png'])
        fig.savefig(fig_name, bbox_inches='tight') 
        plt.close()
    return None 


def DescendantSSFR(nsnap_descendants, nsnap_ancestor=15, gv=None, tau=None, fudge=None): 
    '''     
    '''
    inh_time = time.time() 
    inh = Inherit(nsnap_descendants,  nsnap_ancestor=nsnap_ancestor, gv=gv, tau=tau, fudge=fudge, quiet=True)
    des_dict = inh() 
    print time.time() - inh_time 

    for n_d in des_dict.keys():
        descendant = des_dict[str(n_d)]

        fig_name = ''.join([
            'figure/test/', 
            'SSFR',
            '.nsnap', str(n_d), 
            '.', ''.join([str(nd) for nd in nsnap_descendants]),
            '.png'])

        descendant.plotSsfr(line_color='red', line_width=4, 
                sfms_prop=inh.sfr_prop['sfms'], z=descendant.zsnap, 
                groupcat=True, savefig=fig_name)
        plt.close()
    return None


def dutycycle_profile(): 
    ''' Profile how long simulations take with versus without SF duty cycle 
    '''
    
    duty_list = [
            {'name': 'notperiodic'}, 
            {'freq_range': [6.283185307179586, 62.83185307179586], 'name': 'newamp_squarewave', 'sigma': 0.3}]

    abc_plot = PlotABC(7, abcrun='multifq_wideprior', prior_name='updated')    
    med_theta = abc_plot.med_theta
    
    for duty_dict in duty_list: 
        start_time = time.time() 

        gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = med_theta 
        inh = Inherit([1,3,6], 
                nsnap_ancestor=15,
                subhalo_prop={'scatter': 0.2, 'source': 'li-march'}, 
                sfr_prop={
                    'fq': {'name': 'wetzel'}, 
                    'sfms': {'name': 'linear', 'zslope': 1.14}, 
                    'gv': {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
                    }, 
                evol_prop={
                    'mass': {'name': 'sham'}, 
                    'sfr': {'dutycycle': duty_dict},
                    'type': 'simult', 
                    'fudge': {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset},
                    'tau': {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
                    }
                )
        des_dict = inh() 
        print duty_dict['name'], ' takes ',  time.time() - start_time, ' seconds'



if __name__=='__main__': 


    dutycycle_profile()
    #DescendantQAplot
    #        gv=[0.3, 0.3], tau=[-0.72827483846612151, 0.57825125514919362], fudge=[-1.25, 1.75])
    splashback([1], nsnap_ancestor=15, 
            gv=[0.3, 0.3], tau=[-0.72827483846612151, 0.57825125514919362], fudge=[-1.25, 1.75])
    #DescendantSSFR([1], nsnap_ancestor=15, 
    #        gv=[0.3, 0.3], tau=[-0.72827483846612151, 0.57825125514919362], fudge=[-1.25, 1.75])
    #DescendantQAplot([5], nsnap_ancestor=15, 
    #        gv=[0.3, 0.3], tau=[-0.72827483846612151, 0.57825125514919362], fudge=[-1.25, 1.75])
    #DescendantSSFR([5], nsnap_ancestor=15, 
    #        gv=[0.3, 0.3], tau=[-0.72827483846612151, 0.57825125514919362], fudge=[-1.25, 1.75])
    #
    #DescendantQAplot([1, 5], nsnap_ancestor=15, 
    #        gv=[0.3, 0.3], tau=[-0.72827483846612151, 0.57825125514919362], fudge=[-1.25, 1.75])
    #DescendantSSFR([1, 5], nsnap_ancestor=15, 
    #        gv=[0.3, 0.3], tau=[-0.72827483846612151, 0.57825125514919362], fudge=[-1.25, 1.75])
