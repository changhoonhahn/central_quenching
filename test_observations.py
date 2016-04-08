'''

Module to test groupcat.py module 


'''
import numpy as np

from observations import GroupCat
from observations import PrimusSDSS
from observations import ObservedSFMS
from observations import FitObservedSFMS
from observations import ObservedSSFR
from observations import ObservedSSFR_Peaks
from observations import FitObservedSSFR_Peaks
from observations import Lee2015_SFMS_zslope

from sfr_evol import AverageLogSFR_sfms
from sfr_evol import ScatterLogSFR_sfms

from gal_prop import Fq

from util import util

# plotting
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def BuildGroupCat(Mrcut=18, position='central'): 
    ''' Build Group Catalog 
    '''
    grpcat = GroupCat(Mrcut=Mrcut, position=position)
    grpcat.CompileGroupCat()
    grpcat.Write()

    return None

def GroupCat_iSEDfitMatch(Mrcut=18, position='central'): 
    ''' Test _iSEDfitMatch function of Group Catalog class object
    '''
    grpcat = GroupCat(Mrcut=Mrcut, position=position)
    grpcat.Read()
    grpcat._iSEDfitMatch()

    return None

def PlotObservedSFMS(Mfid=10.5, isedfit=False, sfms_prop=None):
    ''' Plot the observed SFMS relation and the linear best fit to the relation. 
    '''
    prettyplot()
    fig = plt.figure(figsize=(24,5))
    pltsub = fig.add_subplot(111)
    
    zbins = [0.03, 0.1, 0.3, 0.5, 0.7, 0.9]
    for i_z, z in enumerate(zbins): 
        # preset kwargs for group and SDSS+PRIMUS catalogs
        if z == 0.03: 
            observable = 'groupcat'
            obs_str = 'Group Catalog'
            if not isedfit: 
                kwargs = {'Mrcut': 18, 'position': 'central'}
                file_flag = ''
            else: 
                kwargs = {'Mrcut': 18, 'position': 'central', 'isedfit': True}
                obs_str += ' iSEDfit'
                file_flag = 'isedfit'
        else: 
            observable = 'sdssprimus'
            obs_str = 'iSEDfit z='+str(round(z, 2))
            kwargs = {'redshift': z, 'environment': 'no'} 

        mass, muSFR, sigSFR, ngal = ObservedSFMS(observable, **kwargs) 
        bestfit = FitObservedSFMS(observable, Mfid=Mfid, **kwargs) 

        sub = fig.add_subplot(1, len(zbins), i_z+1)

        sub.fill_between(mass, 
                muSFR-sigSFR, 
                muSFR+sigSFR, color=prettycolors()[i_z], label='Observed')
        mrange = np.arange(9.0, 12.0, 0.1)
        sub.plot(mrange, util.line(mrange-Mfid, bestfit), c='k', ls='--', lw=4, 
                label='Best-fit') 
        if sfms_prop is not None: 
            plt.plot(mrange, AverageLogSFR_sfms(mrange, z, sfms_prop=sfms_prop), 
                    c='k', ls='-.', lw=2, label='Model')
        # plot description 
        sub.text(9.15, -4.5, obs_str, fontsize=15)

        sub.set_ylim([-5.0, 2.5]) 
        sub.set_xlim([9.0, 12.0])

        plt.xticks([9., 10., 11.])
        if i_z != 0: 
            sub.set_yticklabels([])
        else: 
            sub.legend(loc='upper left') 
   
    pltsub.spines['top'].set_color('none')
    pltsub.spines['bottom'].set_color('none')
    pltsub.spines['left'].set_color('none')
    pltsub.spines['right'].set_color('none')
    pltsub.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    pltsub.set_ylabel("log SFR", fontsize=25)
    pltsub.set_xlabel('log M*', fontsize=25)
    
    fig.subplots_adjust(wspace=0., hspace=0.)
    fig_file = ''.join([
        'figure/test/'
        'observedSFMS', 
        '.groupcat', file_flag,
        '.sdssprimus',
        '.bestfits.png'
        ])

    fig.savefig(fig_file, bbox_inches='tight', dpi=150)

def PlotObservedSSFR(observable, isedfit=False, Peak=False): 
    '''
    '''
    prettyplot()
    if Peak and observable != 'groupcat': 
        raise ValeuError

    mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    if observable == 'groupcat': 
        zbins = [0.03]    
    elif observable == 'sdssprimus': 
        zbins = [0.1, 0.3, 0.5, 0.7, 0.9]

    qf = Fq() 

    for i_z, z in enumerate(zbins): 
        fig = plt.figure(figsize=(16,16))
        fig.subplots_adjust(hspace=0., wspace=0.)
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        # preset kwargs for group and SDSS+PRIMUS catalogs
        if observable == 'groupcat':
            obs_str = 'Group Catalog'
            file_flag = observable
            if not isedfit: 
                kwargs = {'Mrcut': 18, 'position': 'central'}
                file_flag += ''
            else: 
                kwargs = {'Mrcut': 18, 'position': 'central', 'isedfit': True}
                obs_str += ' iSEDfit'
                file_flag += 'isedfit'
        else: 
            obs_str = 'iSEDfit z='+str(round(z, 2))
            kwargs = {'redshift': z, 'environment': 'no'} 
        # SSFR SF peak fit
        if Peak: 
            sfpeakfit = FitObservedSSFR_Peaks(observable=observable, sfq='star-forming', **kwargs)
            qpeakfit = FitObservedSSFR_Peaks(observable=observable, sfq='quiescent', **kwargs)

            print 'SF Slope ', sfpeakfit[0], ', Offset ', sfpeakfit[1]
            print 'Q Slope ', qpeakfit[0], ', Offset ', qpeakfit[1]

        ssfr_bin_mid, ssfr_dist = ObservedSSFR(observable, **kwargs) 

        for i_mass, mbin in enumerate(mass_bins): 
            subs[i_mass].plot(ssfr_bin_mid[i_mass], ssfr_dist[i_mass], color='k', lw=4, ls='-', label=None)
            if Peak: 
                ssfr_sfpeak = sfpeakfit[0]*(np.mean(mbin)-10.5) + sfpeakfit[1] - np.mean(mbin)
                subs[i_mass].vlines(ssfr_sfpeak, 0.0, 100., lw=3, linestyle='--', color='blue')
                
                ssfr_qpeak = qpeakfit[0]*(np.mean(mbin)-10.5) + qpeakfit[1] - np.mean(mbin)
                subs[i_mass].vlines(ssfr_qpeak, 0.0, 100., lw=3, linestyle='--', color='red')

            subs[i_mass].vlines(
                    qf.SFRcut(np.array([np.mean(mbin)]), z, 
                        sfms_prop={'name': 'kinked', 'mslope_lowmass': 0.7, 'zslope': 1.5}
                        )-np.mean(mbin), 
                    0., 100., 
                    lw=3, linestyle='--', color='k')

            subs[i_mass].set_xlim([-13.0, -7.0])
            subs[i_mass].set_ylim([0.0, 1.6])
            
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M_{*} = [', 
                str(mass_bins[i_mass][0]), ',\;', 
                str(mass_bins[i_mass][1]), ']}$'
                ])
            subs[i_mass].text(-10.5, 1.4, massbin_str,
                    fontsize=24
                    )

            if i_mass == 0: 
                subs[i_mass].set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=20) 
                subs[i_mass].set_xticklabels([])
            elif i_mass == 1: 
                subs[i_mass].set_xticklabels([])
                subs[i_mass].set_yticklabels([])
            elif i_mass == 2:
                #sub.set_yticklabels([])
                subs[i_mass].set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=20) 
                subs[i_mass].set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=20) 
            else: 
                subs[i_mass].set_yticklabels([])
                subs[i_mass].set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=20) 

                #subs[i_mass].legend(loc='lower right', frameon=False)

        fig_file = ''.join([
            'figure/test/'
            'observedSSFR', 
            observable, 
            '.z', str(round(z,2)),
            '.png'
            ])

        fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    return None

def PlotLee2015_SFMS_zdep(): 
    ''' Plot the S0 term (redshift dependence term of the SFMS parmaterization) 
    of the Lee et al. (2015) SFMS fits. 
    '''
    z_mid = np.array([0.36, 0.55, 0.70, 0.85, 0.99, 1.19])
    S0 = np.array([0.80, 0.99, 1.23, 1.35, 1.53, 1.72])
    S0_err = np.array([0.019, 0.015, 0.016, 0.014, 0.017, 0.024])

    zslope, const = Lee2015_SFMS_zslope()

    prettyplot()
    pretty_colors = prettycolors() 

    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.errorbar(z_mid, S0, yerr=S0_err, c='k', elinewidth=3, capsize=10, label='Lee et al. (2015)') 
    sub.plot(z_mid, zslope * (z_mid - 0.05) + const, c=pretty_colors[1], lw=2, ls='--', label='Best fit')
    
    ztext = '\n'.join(['Redshift slope', r'$\mathtt{A_z='+str(round(zslope,2))+'}$'])

    sub.text(0.1, 1.5, ztext, fontsize=20)

    sub.legend(loc='lower right') 
    sub.set_ylabel('Lee et al. (2015) $S_0$ term', fontsize=25) 
    sub.set_xlim([0.0, 1.5]) 
    sub.set_xlabel('Redshift $(\mathtt{z})$', fontsize=25)

    fig_file = ''.join(['figure/', 'Lee2015_SFMS_zdep.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    return None

    
if __name__=="__main__": 
    PlotLee2015_SFMS_zdep()

    #PlotObservedSSFR('groupcat', isedfit=False, Peak=True)

    #GroupCat_iSEDfitMatch(Mrcut=18, position='central')
    #[BuildGroupCat(Mrcut=Mr, position='central') for Mr in [18, 19, 20]]
    #PlotObservedSSFR(isedfit=False)
    #print PlotObservedSFMS(isedfit=True,
    #        sfms_prop={'name': 'linear', 'mslope': 0.55, 'zslope': 1.1}
    #        )
