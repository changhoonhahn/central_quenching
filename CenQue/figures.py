'''

Figures for presentations or proposals


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
from gal_prop import Ssfr
from observations import GroupCat 
from observations import FqCen_bestfit

import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def SSFRevol(t, abcrun, prior_name='try0'): 
    ''' Demonstrate the SFR assignment scheme through the 
    SSFR distribution function P(SSFR) evolution over multiple
    snapshots. Also include extreme versions of the quenching
    timescale for comparison. 
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

    nsnaps = range(1,13)#[1,4,7,10,12]
    color_scheme = ['#6c1500', '#b42400', '#c34f32', '#d27b66', '#e1a799'][::-1]
    #['#cc2800', '#991e00', '#661400', '#330a00', '#000000']

    #nsnaps = [1,4,10]
    print 'starting inherit'
    inh = Inherit(nsnaps, 
            nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
            subhalo_prop=sim_kwargs['subhalo_prop'], 
            sfr_prop=sim_kwargs['sfr_prop'], 
            evol_prop=sim_kwargs['evol_prop'])
    print 'finished inherit1'
    anc = inh.ancestor
    print 'finished inherit2'
    des_dict = inh()
    print 'finished inherit3'

    #prettyplot() 
    #pretty_colors = prettycolors()
        
    print 'ssfr'
    bin_mid, ssfr_dist = anc.Ssfr()
    print 'ssfr calculated'
    output = {} 
    output['anc_ssfr'] =  [bin_mid, ssfr_dist]

    for ii, i_snap in enumerate(nsnaps[::-1]): 
        bin_mid, ssfr_dist = des_dict[str(i_snap)].Ssfr()
        output['snap'+str(i_snap)+'_ssfr'] = [bin_mid, ssfr_dist]
        output['snap'+str(i_snap)+'_z'] = str(round(des_dict[str(i_snap)].zsnap,1))

    pickle.dump(output, open('SSFRevol_demo.p', 'wb'))
    return None 


if __name__=="__main__": 
    SSFRevol(7, 'multirho_inh', prior_name='try0')
    #SSFRevol(6, 'RHOssfrfq_TinkerFq_Std', prior_name='updated')
