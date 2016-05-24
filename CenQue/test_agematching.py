'''
'''

import numpy as np

# Local ---
from gal_prop import SMF
from satellite import SGPop 
from agematching import AllSubhalos

import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors 



def AllSubhalosSMF(scatter=0.2, source='li-march', age_indicator=None): 
    ''' Compare the SMF of the AllSubhalo object to analytic SMF 
    '''
    # import AllSubhalo object
    owl = AllSubhalos(scatter=scatter, source=source)
    owl.Build(age_indicator=age_indicator)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    smf = SMF()
    
    # AllSubhalo SMF
    owl.mass = getattr(owl, 'm.star') 
    subh_mass, subh_phi = smf.Obj(owl, dlogm=0.1, LF=False)
    sub.plot(subh_mass, subh_phi, lw=4, c=pretty_colors[3]) 
    subh_mass, subh_phi = smf._smf(owl.mass[np.where(owl.sfr != -999.)], dlogm=0.1)
    sub.plot(subh_mass, subh_phi, lw=1, c='k') 
    sub.vlines(9.7, subh_phi.min(), subh_phi.max(), linestyle='--', color='k')
    
    # Analytic SMF
    analytic_mass, analytic_phi = smf.analytic(0.05, source=source) 
    sub.plot(analytic_mass, analytic_phi, lw=4, ls='--', c='k', label='Analytic')
    
    # y-axis
    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_ylabel(r'$\Phi$', fontsize=25)
    # x-axis
    sub.set_xlim([6.0, 12.0])
    sub.set_xlabel(r'$\mathtt{M_*}$', fontsize=25)

    sub.legend(loc='lower left') 
    
    if age_indicator is None: 
        age_str = 'NOagematch'
    else: 
        age_str = 'agematch_'+age_indicator

    fig_file = ''.join(['figure/test/', 
        'AllSubhalosSMF', 
        '.sham', 
        '.', str(scatter),
        '.', source, 
        '.', age_str, 
        '.png']) 

    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close() 
    return None 


def AllSubhaloSFMS(scatter=0.2, source='li-march', age_indicator=None): 
    '''
    '''
    # import AllSubhalo object
    owl = AllSubhalos(scatter=scatter, source=source)
    owl.Build(age_indicator=age_indicator)

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(16,8))

    if age_indicator is None: 
        sub1 = fig.add_subplot(121) # centrals
        sub1.scatter(getattr(owl, 'm.star')[owl.cen_index], owl.sfr[owl.cen_index], c=pretty_colors[3], lw=0)
        sub2 = fig.add_subplot(122) # satellites
        sub2.scatter(getattr(owl, 'm.star')[owl.sat_index], owl.sfr[owl.sat_index], c=pretty_colors[5], lw=0)
    else: 
        sub1 = fig.add_subplot(121) # centrals
        cen = sub1.scatter(getattr(owl, 'm.star')[owl.cen_index], owl.sfr[owl.cen_index], c=owl.central_age, cmap='viridis', lw=0) 
        plt.colorbar(cen) 
        sub2 = fig.add_subplot(122) # satellites 
        sub2.scatter(getattr(owl, 'm.star')[owl.sat_index], owl.sfr[owl.sat_index], c=pretty_colors[5], lw=0) 
    
    sub1.text(9.5, -4., 'Central', fontsize=25)
    sub2.text(9.5, -4., 'Satellite', fontsize=25)

    sub1.set_xlim([9.0, 12.0])
    sub1.set_xlabel(r'$\mathtt{M}_*$', fontsize=25)
    sub1.set_ylim([-5.0, 2.0])
    sub1.set_xlabel('SFR', fontsize=25)

    sub2.set_xlim([9.0, 12.0])
    sub2.set_xlabel(r'$\mathtt{M}_*$', fontsize=25)
    sub2.set_ylim([-5.0, 2.0])
    sub2.set_yticklabels([]) 

    if age_indicator is None: 
        age_str = 'NOagematch'
    else: 
        age_str = 'agematch_'+age_indicator

    fig_file = ''.join(['figure/test/', 
        'AllSubhalosSFMS', 
        '.sham', 
        '.', str(scatter),
        '.', source, 
        '.', age_str, 
        '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None



if __name__=='__main__': 
    #AllSubhaloSFMS(age_indicator=None)
    AllSubhaloSFMS(age_indicator='acc80')
