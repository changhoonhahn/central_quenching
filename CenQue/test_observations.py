'''

Module to test groupcat.py module 


'''
import numpy as np
from scipy import interpolate

from observations import GroupCat
from observations import PrimusSDSS
from observations import ObservedSFMS
from observations import FitObservedSFMS
from observations import ObservedSSFR
from observations import ObservedSSFR_Peaks
from observations import FitObservedSSFR_Peaks
from observations import Lee2015_SFMS_zslope
from observations import FqCen_bestfit

from sfr_evol import AverageLogSFR_sfms
from sfr_evol import ScatterLogSFR_sfms
from sfr_evol import AverageLogSSFR_q_peak
from sfr_evol import ScatterLogSSFR_q_peak

from gal_prop import Fq
from gal_prop import Ssfr

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

#def Build_CombinedCatalog(Mrcut=18):
#    ''' Compile galaxy properties derived from different methods for the VAGC catalog from Jeremy 
#    and, for now, output to pickle file
#    '''
#    grpcat = GroupCat(Mrcut=Mrcut, position='all')  # all group catalog 
#    grpcat.Read() 
#    grpcat._iSEDfitMatch()
#    grpcat._Match_OtherSFR(lit=lit)
#
#    return None


def GroupCat_iSEDfitMatch(Mrcut=18, position='central'): 
    ''' Test _iSEDfitMatch function of Group Catalog class object
    '''
    grpcat = GroupCat(Mrcut=Mrcut, position=position)
    grpcat.Read()
    grpcat._iSEDfitMatch()

    return None

def GroupCat_SalimMatch(galprop, lit='uv', Mrcut=18, position='central'): 
    ''' 
    '''
    grpcat = GroupCat(Mrcut=Mrcut, position=position)
    grpcat.Read()
    grpcat._Match_OtherSFR(lit=lit)
    hasmatch = np.where(getattr(grpcat, lit+'_match') >= 0)[0]

    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(10,10))
    bkgd = fig.add_subplot(111, frameon=False)
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    fig.subplots_adjust(hspace=0., wspace=0.)

    if galprop == 'mass': 
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

        for i_mass, mbin in enumerate(mass_bins): 
            mbin = np.where((grpcat.mass[hasmatch] > mass_bins[i_mass][0]) & (grpcat.mass[hasmatch] < mass_bins[i_mass][-1])) 
            subs[i_mass].scatter(grpcat.mass[hasmatch[mbin]], getattr(grpcat, lit+'_mass')[hasmatch[mbin]], color=pretty_colors[3], s=4) 
            subs[i_mass].plot([9., 12.], [9., 12.], color='k', lw=4, ls='--') 
            subs[i_mass].set_xlim([9., 12.]) 
            subs[i_mass].set_ylim([9., 12.]) 
            subs[i_mass].set_xticks([9., 10., 11., 12.])
            subs[i_mass].set_yticks([9., 10., 11., 12.])
    
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M^{VAGC}_{*} = [', 
                str(mass_bins[i_mass][0]), ',\;', 
                str(mass_bins[i_mass][1]), ']}$'
                ])
            subs[i_mass].text(9.2, 11.6, massbin_str, fontsize=15)

            if i_mass == 0: 
                subs[i_mass].set_xticklabels([])
            elif i_mass == 1: 
                subs[i_mass].set_xticklabels([])
                subs[i_mass].set_yticklabels([])
            elif i_mass == 3: 
                subs[i_mass].set_yticklabels([])
    
        bkgd.set_ylabel(lit+' $\mathcal{M}_*$', fontsize=30) 
        bkgd.set_xlabel('VAGC $\mathcal{M}_*$', fontsize=30) 

    elif galprop == 'ssfr':
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

        for i_mass, mbin in enumerate(mass_bins): 
            mbin = np.where((grpcat.mass[hasmatch] > mass_bins[i_mass][0]) & (grpcat.mass[hasmatch] < mass_bins[i_mass][-1])) 
            subs[i_mass].scatter(grpcat.ssfr[hasmatch[mbin]], getattr(grpcat, lit+'_ssfr')[hasmatch[mbin]], color=pretty_colors[3], s=4) 
            subs[i_mass].plot([-13., -9.], [-13., -9.], color='k', lw=4, ls='--') 
            subs[i_mass].set_xlim([-13., -9.]) 
            subs[i_mass].set_ylim([-13., -9.]) 
    
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M^{VAGC}_{*} = [', 
                str(mass_bins[i_mass][0]), ',\;', 
                str(mass_bins[i_mass][1]), ']}$'
                ])
            subs[i_mass].text(-12.8, -9.5, massbin_str, fontsize=15)
            subs[i_mass].set_xticks([-13, -12, -11, -10, -9])
            subs[i_mass].set_yticks([-13, -12, -11, -10, -9])

            if i_mass == 0: 
                subs[i_mass].set_xticklabels([])
            elif i_mass == 1: 
                subs[i_mass].set_xticklabels([])
                subs[i_mass].set_yticklabels([])
            elif i_mass == 3: 
                subs[i_mass].set_yticklabels([])
    
        bkgd.set_ylabel(lit+' SSFR', fontsize=30) 
        bkgd.set_xlabel('VAGC SSFR ', fontsize=30) 
    
    elif galprop == 'sfr':
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

        for i_mass, mbin in enumerate(mass_bins): 
            mbin = np.where((grpcat.mass[hasmatch] > mass_bins[i_mass][0]) & (grpcat.mass[hasmatch] < mass_bins[i_mass][-1])) 
            subs[i_mass].scatter(grpcat.sfr[hasmatch[mbin]], getattr(grpcat, lit+'_sfr')[hasmatch[mbin]], color=pretty_colors[3], s=4) 
            subs[i_mass].plot([-3., 2.], [-3., 2.], color='k', lw=4, ls='--') 
            subs[i_mass].set_xlim([-3., 3.]) 
            subs[i_mass].set_ylim([-3., 3.]) 
    
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M^{VAGC}_{*} = [', 
                str(mass_bins[i_mass][0]), ',\;', 
                str(mass_bins[i_mass][1]), ']}$'
                ])
            subs[i_mass].text(-12.8, -9.5, massbin_str, fontsize=15)
            subs[i_mass].set_xticks([-3, -2, -1, 0, 1, 2, 3])
            subs[i_mass].set_yticks([-3, -2, -1, 0, 1, 2, 3])

            if i_mass == 0: 
                subs[i_mass].set_xticklabels([])
            elif i_mass == 1: 
                subs[i_mass].set_xticklabels([])
                subs[i_mass].set_yticklabels([])
            elif i_mass == 3: 
                subs[i_mass].set_yticklabels([])
    
        bkgd.set_ylabel(lit+' SFR', fontsize=30) 
        bkgd.set_xlabel('VAGC SFR ', fontsize=30) 

    fig_file = ''.join([
        'figure/test/', 
        'Salim_Comparison',
        '.Mr', str(Mrcut), 
        '.', position, 
        '.', galprop, 
        '.', lit, 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    return None

def GroupCat_SalimMatch_SFMS(galprop, lit='uv', Mrcut=18, position='central'): 
    ''' 
    '''
    grpcat = GroupCat(Mrcut=Mrcut, position=position)
    grpcat.Read()
    grpcat._Match_OtherSFR(lit=lit)
    #grpcat._match_NSA_UVssfr()
    hasmatch = np.where(getattr(grpcat, lit+'_match') >= 0)[0]

    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(10,10))
    bkgd = fig.add_subplot(111, frameon=False)
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    fig.subplots_adjust(hspace=0., wspace=0.)

    if galprop == 'mass': 
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

        for i_mass, mbin in enumerate(mass_bins): 
            mbin = np.where((grpcat.mass[hasmatch] > mass_bins[i_mass][0]) & (grpcat.mass[hasmatch] < mass_bins[i_mass][-1])) 
            subs[i_mass].scatter(grpcat.mass[hasmatch[mbin]], getattr(grpcat, lit+'_mass')[hasmatch[mbin]], color=pretty_colors[3], s=4) 
            subs[i_mass].plot([9., 12.], [9., 12.], color='k', lw=4, ls='--') 
            subs[i_mass].set_xlim([9., 12.]) 
            subs[i_mass].set_ylim([9., 12.]) 
            subs[i_mass].set_xticks([9., 10., 11., 12.])
            subs[i_mass].set_yticks([9., 10., 11., 12.])
    
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M^{VAGC}_{*} = [', 
                str(mass_bins[i_mass][0]), ',\;', 
                str(mass_bins[i_mass][1]), ']}$'
                ])
            subs[i_mass].text(9.2, 11.6, massbin_str, fontsize=15)

            if i_mass == 0: 
                subs[i_mass].set_xticklabels([])
            elif i_mass == 1: 
                subs[i_mass].set_xticklabels([])
                subs[i_mass].set_yticklabels([])
            elif i_mass == 3: 
                subs[i_mass].set_yticklabels([])
    
        bkgd.set_ylabel(lit+' $\mathcal{M}_*$', fontsize=30) 
        bkgd.set_xlabel('VAGC $\mathcal{M}_*$', fontsize=30) 

    elif galprop == 'ssfr':
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

        for i_mass, mbin in enumerate(mass_bins): 
            mbin = np.where((grpcat.mass[hasmatch] > mass_bins[i_mass][0]) & (grpcat.mass[hasmatch] < mass_bins[i_mass][-1])) 
            subs[i_mass].scatter(grpcat.ssfr[hasmatch[mbin]], getattr(grpcat, lit+'_ssfr')[hasmatch[mbin]], color=pretty_colors[3], s=4) 
            subs[i_mass].plot([-13., -9.], [-13., -9.], color='k', lw=4, ls='--') 
            subs[i_mass].set_xlim([-13., -9.]) 
            subs[i_mass].set_ylim([-13., -9.]) 
    
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M^{VAGC}_{*} = [', 
                str(mass_bins[i_mass][0]), ',\;', 
                str(mass_bins[i_mass][1]), ']}$'
                ])
            subs[i_mass].text(-12.8, -9.5, massbin_str, fontsize=15)
            subs[i_mass].set_xticks([-13, -12, -11, -10, -9])
            subs[i_mass].set_yticks([-13, -12, -11, -10, -9])

            if i_mass == 0: 
                subs[i_mass].set_xticklabels([])
            elif i_mass == 1: 
                subs[i_mass].set_xticklabels([])
                subs[i_mass].set_yticklabels([])
            elif i_mass == 3: 
                subs[i_mass].set_yticklabels([])
    
        bkgd.set_ylabel(lit+' SSFR', fontsize=30) 
        bkgd.set_xlabel('VAGC SSFR ', fontsize=30) 
    
    elif galprop == 'sfr':
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]

        for i_mass, mbin in enumerate(mass_bins): 
            mbin = np.where((grpcat.mass[hasmatch] > mass_bins[i_mass][0]) & (grpcat.mass[hasmatch] < mass_bins[i_mass][-1])) 
            subs[i_mass].scatter(grpcat.sfr[hasmatch[mbin]], getattr(grpcat, lit+'_sfr')[hasmatch[mbin]], color=pretty_colors[3], s=4) 
            subs[i_mass].plot([-3., 2.], [-3., 2.], color='k', lw=4, ls='--') 
            subs[i_mass].set_xlim([-3., 3.]) 
            subs[i_mass].set_ylim([-3., 3.]) 
    
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M^{VAGC}_{*} = [', 
                str(mass_bins[i_mass][0]), ',\;', 
                str(mass_bins[i_mass][1]), ']}$'
                ])
            subs[i_mass].text(-12.8, -9.5, massbin_str, fontsize=15)
            subs[i_mass].set_xticks([-3, -2, -1, 0, 1, 2, 3])
            subs[i_mass].set_yticks([-3, -2, -1, 0, 1, 2, 3])

            if i_mass == 0: 
                subs[i_mass].set_xticklabels([])
            elif i_mass == 1: 
                subs[i_mass].set_xticklabels([])
                subs[i_mass].set_yticklabels([])
            elif i_mass == 3: 
                subs[i_mass].set_yticklabels([])
    
        bkgd.set_ylabel(lit+' SFR', fontsize=30) 
        bkgd.set_xlabel('VAGC SFR ', fontsize=30) 

    fig_file = ''.join([
        'figure/test/', 
        'Salim_Comparison',
        '.Mr', str(Mrcut), 
        '.', position, 
        '.', galprop, 
        '.', lit, 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    return None


def GroupCat_SalimMatch_SSFR(Mrcut=18, position='central'):
    ''' Compare the SSFr distribution
    '''
    grpcat = GroupCat(Mrcut=Mrcut, position=position)
    grpcat.Read()
    grpcat._Match_OtherSFR(lit='salim2016')
    grpcat._Match_OtherSFR(lit='uv')
    
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(16,16))
    fig.subplots_adjust(hspace=0., wspace=0.)
    subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

    mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    for i_mass, mbin in enumerate(mass_bins): 
        in_massbin = np.where(
                (grpcat.mass > mass_bins[i_mass][0]) & 
                (grpcat.mass <= mass_bins[i_mass][1]) & 
                (grpcat.salim2016_match != -999) & 
                (grpcat.uv_match != -999))
            
        vagc_dist, bin_edges = np.histogram(grpcat.ssfr[in_massbin], range=[-13., -7.], bins=40, normed=True)
        salim_dist, bin_edges = np.histogram(grpcat.salim2016_ssfr[in_massbin], range=[-13., -7.], bins=40, normed=True)
        uv_dist, bin_edges = np.histogram(grpcat.uv_ssfr[in_massbin], range=[-13., -7.], bins=40, normed=True)

        subs[i_mass].plot(0.5*(bin_edges[:-1] + bin_edges[1:]), vagc_dist, color='k', lw=4, ls='-') 
        subs[i_mass].plot(0.5*(bin_edges[:-1] + bin_edges[1:]), salim_dist, color=pretty_colors[3], lw=4, ls='--') 
        subs[i_mass].plot(0.5*(bin_edges[:-1] + bin_edges[1:]), uv_dist, color=pretty_colors[5], lw=4, ls=':') 
        #subs[i_mass].plot(salim_bin[i_mass], salim_ssfr[i_mass], color='r', lw=4, ls='--') 
        subs[i_mass].set_xlim([-13.0, -9.])
        subs[i_mass].set_ylim([0.0, 1.6])
        
        massbin_str = ''.join([ 
            r'$\mathtt{log \; M_{*} = [', 
            str(mass_bins[i_mass][0]), ',\;', 
            str(mass_bins[i_mass][1]), ']}$'
            ])
        subs[i_mass].text(-11.5, 1.4, massbin_str,
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
        
    fig_file = ''.join([
        'figure/test/', 
        'SSFR_comparison',
        '.Mr', str(Mrcut), 
        '.', position, 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    #plt.show() 
    return None

def GroupCat_SalimMatch_SFR(Mrcut=18, position='central'):
    ''' Compare the SSFr distribution
    '''
    grpcat = GroupCat(Mrcut=Mrcut, position=position)
    grpcat.Read()
    grpcat._Match_OtherSFR(lit='salim2016')
    grpcat._Match_OtherSFR(lit='uv')
    print len(np.where(grpcat.uv_match != -999)[0])

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(16,16))
    fig.subplots_adjust(hspace=0., wspace=0.)
    subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

    mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    for i_mass, mbin in enumerate(mass_bins): 
        in_massbin = np.where(
                (grpcat.mass > mass_bins[i_mass][0]) & 
                (grpcat.mass <= mass_bins[i_mass][1]) & 
                (grpcat.salim2016_match != -999) & 
                (grpcat.uv_match != -999))
            
        vagc_dist, bin_edges = np.histogram(grpcat.sfr[in_massbin], range=[-3., 3.], bins=40, normed=True)
        salim_dist, bin_edges = np.histogram(grpcat.salim2016_sfr[in_massbin], range=[-3., 3.], bins=40, normed=True)
        uv_dist, bin_edges = np.histogram(grpcat.uv_sfr[in_massbin], range=[-3., 3.], bins=40, normed=True)

        subs[i_mass].plot(0.5*(bin_edges[:-1] + bin_edges[1:]), vagc_dist, color='k', lw=4, ls='-') 
        subs[i_mass].plot(0.5*(bin_edges[:-1] + bin_edges[1:]), salim_dist, color=pretty_colors[3], lw=4, ls='--') 
        subs[i_mass].plot(0.5*(bin_edges[:-1] + bin_edges[1:]), uv_dist, color=pretty_colors[5], lw=4, ls=':') 
        #subs[i_mass].plot(salim_bin[i_mass], salim_ssfr[i_mass], color='r', lw=4, ls='--') 
        subs[i_mass].set_xlim([-3., 3.])
        subs[i_mass].set_ylim([0.0, 1.6])
        
        massbin_str = ''.join([ 
            r'$\mathtt{log \; M_{*} = [', 
            str(mass_bins[i_mass][0]), ',\;', 
            str(mass_bins[i_mass][1]), ']}$'
            ])
        subs[i_mass].text(-1.5, 1.4, massbin_str,
                fontsize=24
                )

        if i_mass == 0: 
            subs[i_mass].set_ylabel(r'$\mathtt{P(log \; SFR)}$', fontsize=20) 
            subs[i_mass].set_xticklabels([])
        elif i_mass == 1: 
            subs[i_mass].set_xticklabels([])
            subs[i_mass].set_yticklabels([])
        elif i_mass == 2:
            #sub.set_yticklabels([])
            subs[i_mass].set_ylabel(r'$\mathtt{P(log \; SFR)}$', fontsize=20) 
            subs[i_mass].set_xlabel(r'$\mathtt{log \; SFR \;[yr^{-1}]}$', fontsize=20) 
        else: 
            subs[i_mass].set_yticklabels([])
            subs[i_mass].set_xlabel(r'$\mathtt{log \; SFR \;[yr^{-1}]}$', fontsize=20) 
        
    fig_file = ''.join([
        'figure/test/', 
        'SFR_comparison',
        '.Mr', str(Mrcut), 
        '.', position, 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    #plt.show() 
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

        mass, muSFR, sigSFR, ngal = ObservedSFMS(observable, sfms_prop=sfms_prop, **kwargs) 
        #bestfit = FitObservedSFMS(observable, Mfid=Mfid, **kwargs) 

        sub = fig.add_subplot(1, len(zbins), i_z+1)

        sub.fill_between(mass, 
                muSFR-sigSFR, 
                muSFR+sigSFR, color=prettycolors()[i_z], label='Observed')
        mrange = np.arange(9.0, 12.0, 0.1)
        #sub.plot(mrange, util.line(mrange-Mfid, bestfit), c='k', ls='--', lw=4, 
        #        label='Best-fit') 
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
    if Peak and 'groupcat' not in observable: 
        raise ValeuError

    mass_bins = [[9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]]
    if 'groupcat' in observable: 
        zbins = [0.03]    
    elif observable == 'sdssprimus': 
        zbins = [0.1, 0.3, 0.5, 0.7, 0.9]

    qf = Fq() 

    for i_z, z in enumerate(zbins): 
        fig = plt.figure(figsize=(16,16))
        fig.subplots_adjust(hspace=0., wspace=0.)
        subs = [fig.add_subplot(2, 2, i_mass+1) for i_mass in xrange(4)]  

        # preset kwargs for group and SDSS+PRIMUS catalogs
        if 'groupcat' in observable:
            if 'sat' in observable: 
                obs_str = 'Group Catalog Satellites'
                censat = 'satellite'
            elif 'cen' in observable: 
                obs_str = 'Group Catalog Centrals' 
                censat = 'central'
            file_flag = observable

            if not isedfit: 
                kwargs = {'Mrcut': 18, 'position': censat}
                file_flag += ''
            else: 
                kwargs = {'Mrcut': 18, 'position': censat, 'isedfit': True}
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

def Plot_fQcentrals():
    ''' Plot the quiescent fractions from Tinker et al. (2013)
    '''
    # mass binnning we impose 
    m_low = np.array([9.5, 10., 10.5, 11., 11.5]) 
    m_high = np.array([10., 10.5, 11., 11.5, 12.0])
    m_mid = 0.5 * (m_low + m_high) 

    # SDSS 
    fq_file = ''.join(['dat/observations/cosmos_fq/', 'fcen_red_sdss_scatter.dat']) 
    m_sdss, fqcen_sdss, N_sdss = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2])
    
    fqcen_sdss_rebin = [] 
    for im, m_mid_i in enumerate(m_mid): 
        sdss_mbin = np.where(
                (m_sdss >= m_low[im]) & 
                (m_sdss < m_high[im])) 
    
        fqcen_sdss_rebin.append(
                np.sum(fqcen_sdss[sdss_mbin] * N_sdss[sdss_mbin].astype('float'))/np.sum(N_sdss[sdss_mbin].astype('float'))
                )
    fqcen_sdss_rebin = np.array(fqcen_sdss_rebin)

    prettyplot()
    pretty_colors = prettycolors()
    
    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.plot([9.]+list(m_mid), [0.]+list(fqcen_sdss_rebin), c='k') 
    
    for iz, z in enumerate([0.36, 0.66, 0.88]): 
        fq_file = ''.join(['dat/observations/cosmos_fq/', 
            'stats_z', str(iz+1), '.fq_cen']) 
        
        m_cosmos, fqcen_cosmos, fqcen_cosmos_low, fqcen_cosmos_high = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2,3])
        m_cosmos = np.log10(m_cosmos)

        fqcen_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos) 
        fqcen_low_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos_low) 
        fqcen_high_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos_high) 
    
        fqcen_cosmos_rebin = fqcen_interp(m_mid)
        fqcen_low_cosmos_rebin = fqcen_low_interp(m_mid)
        fqcen_high_cosmos_rebin = fqcen_high_interp(m_mid)

        sub.fill_between([9.]+list(m_mid), 
                [0.]+list(fqcen_low_cosmos_rebin), [0.]+list(fqcen_high_cosmos_rebin), 
                color=pretty_colors[2*iz+1]) 
    
    sub.set_xlim([8.75, 12.0])
    sub.set_xlabel(r'$\mathtt{M_*}$', fontsize=30) 
    sub.set_ylim([0., 1.0])
    sub.set_ylabel(r'$\mathtt{f_Q^{cen}}$', fontsize=30) 
    plt.show() 
        
def Plot_fQcen_parameterized():
    ''' Plot the bestfit parameterization of the quiescent fractions 
    from Tinker et al. (2013)
    '''
    # best fit alpha(M*) values 
    m_mid, alpha_m = FqCen_bestfit(clobber=True) 

    # mass binnning we impose 
    m_low = np.array([9.5, 10., 10.5, 11., 11.5]) 
    m_high = np.array([10., 10.5, 11., 11.5, 12.0])
    m_mid = 0.5 * (m_low + m_high) 

    # SDSS 
    fq_file = ''.join(['dat/observations/cosmos_fq/', 'fcen_red_sdss_scatter.dat']) 
    m_sdss, fqcen_sdss, N_sdss = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2])
    
    fqcen_sdss_rebin = [] 
    for im, m_mid_i in enumerate(m_mid): 
        sdss_mbin = np.where(
                (m_sdss >= m_low[im]) & 
                (m_sdss < m_high[im])) 
    
        fqcen_sdss_rebin.append(
                np.sum(fqcen_sdss[sdss_mbin] * N_sdss[sdss_mbin].astype('float'))/np.sum(N_sdss[sdss_mbin].astype('float'))
                )
    fqcen_sdss_rebin = np.array(fqcen_sdss_rebin)

    fqcen_cosmos_rebin = [] 
    fqcen_low_cosmos_rebin = [] 
    fqcen_high_cosmos_rebin = [] 
    for iz, z in enumerate([0.36, 0.66, 0.88]): 
        fq_file = ''.join(['dat/observations/cosmos_fq/', 
            'stats_z', str(iz+1), '.fq_cen']) 
        
        m_cosmos, fqcen_cosmos, fqcen_cosmos_low, fqcen_cosmos_high = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2,3])
        m_cosmos = np.log10(m_cosmos)

        fqcen_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos) 
        fqcen_low_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos_low) 
        fqcen_high_interp = interpolate.interp1d(m_cosmos, fqcen_cosmos_high) 
    
        fqcen_cosmos_rebin.append( fqcen_interp(m_mid) ) 
        fqcen_low_cosmos_rebin.append( fqcen_low_interp(m_mid) ) 
        fqcen_high_cosmos_rebin.append( fqcen_high_interp(m_mid) ) 

    fqcen_cosmos_rebin = np.array(fqcen_cosmos_rebin) 
    fqcen_low_cosmos_rebin = np.array(fqcen_low_cosmos_rebin)
    fqcen_high_cosmos_rebin = np.array(fqcen_high_cosmos_rebin) 

    z_arr = np.array([0.0, 0.36, 0.66, 0.88])
    
    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure()#figsize=(20,5)) 
    sub = fig.add_subplot(111)#, len(m_mid), im+1) 
    for im in range(len(m_mid)): 
        sub.fill_between(z_arr, 
                np.array([fqcen_sdss_rebin[im]]+list(fqcen_low_cosmos_rebin[:,im])),
                np.array([fqcen_sdss_rebin[im]]+list(fqcen_high_cosmos_rebin[:,im])),
                color=pretty_colors[2*im], 
                alpha=0.5, 
                label=r'$\mathtt{log\;M_*} = '+str(m_mid[im])+'$')
        sub.plot(
                z_arr,
                fqcen_sdss_rebin[im] * (1. + z_arr)**alpha_m[im], 
                c='k', ls='--', lw=2)

    sub.legend(loc='best') 
    
    sub.set_xlim([0.0, 1.0])
    #sub.set_ylim([0.0, 1.0])
    sub.set_ylim([0.01, 1.0]) 
    sub.set_yscale('log') 
    sub.set_xlabel(r'$\mathtt{M_*}$', fontsize=30) 
    sub.set_ylabel(r'$\mathtt{f_Q^{cen}}$', fontsize=30) 
    fig_file = ''.join(['figure/test/', 
        'test_fQcen_parameterized.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close() 

def Plot_fQcen_SDSS(): 
    ''' Compare the quiescent fraction of SDSS from the *corrected* SDSS fQ^cen 
    from Tinker et al. (2013) versus the Wetzel et al. (2013) parameterization. 
    '''
    # mass binnning we impose 
    m_low = np.array([9.5, 10., 10.5, 11., 11.5]) 
    m_high = np.array([10., 10.5, 11., 11.5, 12.0])
    m_mid = 0.5 * (m_low + m_high) 

    # SDSS 
    fq_file = ''.join(['dat/observations/cosmos_fq/', 'fcen_red_sdss_scatter.dat']) 
    m_sdss, fqcen_sdss, N_sdss = np.loadtxt(fq_file, unpack=True, usecols=[0,1,2])
    
    fqcen_sdss_rebin = [] 
    for im, m_mid_i in enumerate(m_mid): 
        sdss_mbin = np.where(
                (m_sdss >= m_low[im]) & 
                (m_sdss < m_high[im])) 
    
        fqcen_sdss_rebin.append(
                np.sum(fqcen_sdss[sdss_mbin] * N_sdss[sdss_mbin].astype('float'))/np.sum(N_sdss[sdss_mbin].astype('float'))
                )
    fqcen_sdss_rebin = np.array(fqcen_sdss_rebin)

    prettyplot() 
    pretty_colors = prettycolors()  
    fig = plt.figure() 
    sub = fig.add_subplot(111)
    qf = Fq()
    fqcen_model = qf.model(m_mid, 0.05, lit='cosmos_tinker') 
    sub.plot(m_mid, fqcen_model, 
            c=pretty_colors[3], lw=3, label=r'Wetzel et al. (2013) fit') 
    sub.scatter(m_mid, fqcen_sdss_rebin, 
            color=pretty_colors[0], lw=0, s=40, label=r'Tinker et al. (2013)') 

    sub.set_xlim([9.0, 12.0]) 
    sub.set_xlabel(r'$\mathtt{log\;M_*}$', fontsize=25) 
    sub.set_ylim([0.0, 1.0]) 
    sub.set_ylabel(r'$\mathtt{f_Q^{cen}}$', fontsize=25) 
    sub.legend(loc='upper left', scatterpoints=1, markerscale=3) 
    fig_file = ''.join(['figure/test/', 
        'Fq_central_SDSS.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close() 

def fgas_comparison(): 
    '''
    Compare f_gas (gas fraction) from different literature 
    '''
    m_arr = np.arange(9.0, 12.1, 0.1)   # mass array 

    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(15,7)) 
    fig1 = plt.figure(figsize=(15,7)) 
    fig2 = plt.figure(figsize=(15,7)) 
    for iz, z in enumerate([0.1, 0.4, 0.8]): 
        sub = fig.add_subplot(1, 3, iz+1)
        sub1 = fig1.add_subplot(1, 3, iz+1)
        sub2 = fig2.add_subplot(1, 3, iz+1)
        # Stewart+2009
        f_stargas = 0.04 * ((10.**m_arr)/(4.5*10.**11.))**(-0.59 * (1. + z)**0.45)
        f_gas_Stewart = 1. - 1./(1. + f_stargas)

        # Santini+2014
        if z < 0.2: 
            dat_file = util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.1.dat'
        elif (z >= 0.2) & (z < 0.6): 
            dat_file = util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.4.dat'
        elif (z >= 0.6) & (z < 1.): 
            dat_file = util.code_dir().split('CenQue')[0]+'dat/santini_fgas_z0.8.dat'
        else: 
            raise ValueError
        m_dat, fgas_dat = np.loadtxt(dat_file, delimiter=',', unpack=True, usecols=[0,1]) 
        f_gas_Santini = fgas_dat 
        
        # Boselli+2014
        f_gas_Boselli = 1. - 1./(1. + 10**(-0.69 * m_arr + 6.63))
        
        M_gas_Stewart = np.log10(1./(1./f_gas_Stewart-1.) * 10**m_arr)
        M_gas_Santini = np.log10(1./(1./f_gas_Santini-1.) * 10**m_dat)
        M_gas_Boselli = np.log10(1./(1./f_gas_Boselli-1.) * 10**m_arr)
    
        t_gas_stewart = 10**M_gas_Stewart / 10**AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14})/10**9
        t_gas_santini = 10**M_gas_Santini / 10**AverageLogSFR_sfms(m_dat, z, sfms_prop={'name': 'linear', 'zslope':1.14})/10**9
        t_gas_boselli = 10**M_gas_Boselli / 10**AverageLogSFR_sfms(m_arr, z, sfms_prop={'name': 'linear', 'zslope':1.14})/10**9


        sub.plot(m_arr, f_gas_Stewart, lw=3, c=pretty_colors[1], label='Stewart+2009') 
        sub.plot(m_dat, f_gas_Santini, lw=3, c=pretty_colors[3], label='Santini+2014') 
        if z < 0.2: 
            sub.plot(m_arr, f_gas_Boselli, lw=3, c=pretty_colors[5], label='Boselli+2014') 
        sub.text(10.0, 0.05, r"$\mathtt{z = "+str(z)+"}$", fontsize=25) 
        sub.set_xlim([9.7, 11.5]) 
        if iz == 1: 
            sub.set_xlabel(r'Stellar Mass', fontsize=25) 
        sub.set_ylim([0., 0.5]) 
        if iz == 0: 
            sub.set_ylabel(r'$\mathtt{f_{gas}} = \mathtt{M_{gas}/(M_{gas} + M_{star})}$', fontsize=25) 
            sub.legend(loc='upper right') 
        else: 
            sub.set_yticklabels([]) 


        sub1.plot(m_arr, M_gas_Stewart, lw=3, c=pretty_colors[1], label='Stewart+2009') 
        sub1.plot(m_dat, M_gas_Santini, lw=3, c=pretty_colors[3], label='Santini+2014') 
        if z < 0.2: 
            sub1.plot(m_arr, M_gas_Boselli, lw=3, c=pretty_colors[5], label='Boselli+2014') 
    
        sub1.text(10.0, 9.25, r"$\mathtt{z = "+str(z)+"}$", fontsize=25) 
        sub1.set_xlim([9.7, 11.5]) 
        if iz == 1: 
            sub1.set_xlabel(r'Stellar Mass', fontsize=25) 
        sub1.set_ylim([9., 10.5]) 
        if iz == 0: 
            sub1.set_ylabel(r'$\mathtt{M_{gas}}$', fontsize=25) 
            sub1.legend(loc='upper right') 
        else: 
            sub1.set_yticklabels([]) 


        sub2.plot(m_arr, t_gas_stewart, lw=3, c=pretty_colors[1], label='Stewart+2009') 
        sub2.plot(m_dat, t_gas_santini, lw=3, c=pretty_colors[3], label='Santini+2014') 
        if z < 0.2: 
            sub2.plot(m_arr, t_gas_boselli, lw=3, c=pretty_colors[5], label='Boselli+2014, z=0') 
        
        sub2.text(10.0, 0.05, r"$\mathtt{z = "+str(z)+"}$", fontsize=25) 
        sub2.set_xlim([9.7, 11.5]) 
        if iz == 1: 
            sub2.set_xlabel(r'Stellar Mass', fontsize=25) 
        sub2.set_ylim([0., 10.]) 
        if iz == 0: 
            sub2.set_ylabel(r'$\mathtt{M_{gas}(z)/SFR^{SFMS}(z)}$ [Gyr]', fontsize=25) 
            sub2.legend(loc='upper right') 
        else: 
            sub2.set_yticklabels([]) 

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/',
        'f_gas_comparison',
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()

    fig1.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/',
        'M_gas_comparison',
        '.png'])
    fig1.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()

    fig2.subplots_adjust(wspace=0.0, hspace=0.0)
    fig_file = ''.join(['figure/',
        't_gas_comparison',
        '.png'])
    fig2.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()
    return None

def CenSat_GroupCat(): 
    '''
    '''
    # import centrals
    centrals = GroupCat(Mrcut=18, position='central')
    centrals.Read()

    # import satellites 
    satellites = GroupCat(Mrcut=18, position='satellite')
    satellites.Read()

    N_cen = np.sum(centrals.mass >= 9.7) 
    N_sat = np.sum(satellites.mass >= 9.7) 
    print np.float(N_cen)/np.float(N_cen + N_sat)


if __name__=="__main__": 
    #GroupCat_SalimMatch_SSFR(Mrcut=18, position='central')
    #GroupCat_SalimMatch_SFR(Mrcut=18, position='central')
    #GroupCat_SalimMatch('sfr', lit='uv', Mrcut=18, position='central')
    #GroupCat_SalimMatch('ssfr', lit='uv', Mrcut=18, position='central')
    #GroupCat_SalimMatch('mass', lit='salim2016', Mrcut=18, position='central')
    #GroupCat_SalimMatch('ssfr', lit='salim2016', Mrcut=18, position='central')
    #GroupCat_SalimMatch(Mrcut=18, position='satellite')
    #CenSat_GroupCat()
    #Plot_fQcen_SDSS()
    #Plot_fQcen_parameterized()

    #Plot_fQcentrals()

    #grpcat = GroupCat(Mrcut=18, position='central')
    #grpcat.Read()
    #print np.min(grpcat.z), np.max(grpcat.z)
    #PlotLee2015_SFMS_zdep()
    #fgas_comparison()

    [BuildGroupCat(Mrcut=Mr, position='all') for Mr in [18, 19, 20]]
    #[BuildGroupCat(Mrcut=Mr, position='central') for Mr in [18, 19, 20]]
    #[BuildGroupCat(Mrcut=Mr, position='satellite') for Mr in [18, 19, 20]]
    #PlotObservedSSFR('groupcat_cen', isedfit=False, Peak=True)
    #PlotObservedSSFR('groupcat_sat', isedfit=False, Peak=True)

    #GroupCat_iSEDfitMatch(Mrcut=18, position='central')
    #PlotObservedSSFR(isedfit=False)
    #print PlotObservedSFMS(isedfit=False,
    #        sfms_prop={'name': 'linear', 'mslope': 0.53, 'zslope': 1.14}
    #        )
