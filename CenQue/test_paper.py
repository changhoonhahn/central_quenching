'''

Test pertaining to the central quenching timescale paper

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

from gal_prop import SMF
from observations import GroupCat 
from observations import FqCen_bestfit

import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def SMFevol(): 
    ''' Compare the Marchesini interpolation to others in the literature 
    '''
    #zrange = [0.9, 0.5, 0.05]
    zrange = [0.9]
    source_list = ['li-march', 'li-march-extreme']

    prettyplot() 
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(7, 7))
    sub = fig.add_subplot(111)

    smf = SMF()
    for i_s, source in enumerate(source_list): 
        for iz, z in enumerate(zrange): 
            mass, phi = smf.analytic(z, source=source)  # SMF phi(M*) 
            phi = np.log10(phi)   
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

    for z in [0.9]:#0.75, 0.9, 1.25]: 
        mass, phi = smf.analytic(z, source='muzzin')  # SMF phi(M*) 
        phi = np.log10(phi)   
        if z == 0.9: 
            phi_low = (1.2/16.25+1.)*phi
            phi_high = (13.91-1.23)/(13.91)*phi
            sub.fill_between(mass, phi_low, phi_high)
        else: 
            sub.plot(mass, phi, lw=6, ls=':', c='k')
            
    #sub.set_ylim([10**-3.75, 10**-1.75])
    sub.set_ylim([-3.75, -1.75])
    sub.set_xlim([8.8, 11.5])
    sub.set_xticks([9.0, 10.0, 11.0])
    sub.minorticks_on()
    #sub.set_yscale('log')
    sub.set_xlabel(r'log($\mathtt{M_*} /\mathtt{M_\odot}$)', fontsize=25) 
    sub.set_ylabel(r'log($\mathtt{\Phi / Mpc^{-3}\;dex^{-1}}$)', fontsize=25) 
    sub.legend(loc='lower left', scatteryoffsets=[0.6])
    
    plt.show() 
    
    #fig_file = ''.join(['figure/', 'test_SMF_evol.png'])
    #fig.savefig(fig_file, bbox_inches='tight', dpi=150) 
    return None


if __name__=="__main__": 
    SMFevol() 
