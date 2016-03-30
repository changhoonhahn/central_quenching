'''

Module to test groupcat.py module 


'''
import numpy as np

from observations import GroupCat
from observations import PrimusSDSS
from observations import ObservedSFMS
from observations import FitObservedSFMS

from sfr_evol import AverageLogSFR_sfms
from sfr_evol import ScatterLogSFR_sfms

from util import util

# plotting
import matplotlib.pyplot as plt
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors


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



if __name__=="__main__": 
    #GroupCat_iSEDfitMatch(Mrcut=18, position='central')
    #[BuildGroupCat(Mrcut=Mr, position='central') for Mr in [18, 19, 20]]
    print PlotObservedSFMS(isedfit=True,
            sfms_prop={'name': 'linear', 'mslope': 0.55, 'zslope': 1.1}
            )
