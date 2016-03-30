'''

Test the SFR evolution functions 


'''
import numpy as np

import util.util as Util
from observations import GroupCat
from observations import PrimusSDSS

from sfr_evol import getTauQ
from sfr_evol import AverageLogSFR_sfms
from sfr_evol import ScatterLogSFR_sfms
from sfr_evol import DeltaLogSFR_SF_Q_peak 

from plotting.plots import PlotSFMS

import matplotlib.pyplot as plt
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

def ModelSFMS(observable='sdssprimus'): 
    '''
    '''
    if observable == 'sdssprimus': 
        z_bins = [0.1, 0.3, 0.5, 0.7, 0.9]
    elif observable == 'groupcat': 
        z_bins = [0.1]
    
    for i_z, z_mid in enumerate(z_bins): 

        if observable == 'sdssprimus': 
            obsdata = PrimusSDSS(z_mid)
            obsdata.Read()
            z_bin = z_mid 

        elif observable == 'groupcat': 
            obsdata = GroupCat(Mrcut=18, position='central')
            obsdata.Read()
            z_bin = np.mean(obsdata.z)

        # SDSS/PRIMUS SFMS of 0 environment galaxies
        sfms_plot = PlotSFMS()
        sfms_plot.plot(
                mass=obsdata.mass, 
                sfr=obsdata.sfr, 
                allgal=True, 
                color='blue')#sfms_plot.pretty_colors[i_z])
        
        # Model SFMS 
        mstar = np.arange(7.0, 12.01, 0.01)
        plt.plot(
                mstar, 
                AverageLogSFR_sfms(mstar, z_bin), 
                c='k', ls='-', 
                lw=4, label='Model SFMS')
        plt.plot(mstar, 
                AverageLogSFR_sfms(mstar, z_bin) - ScatterLogSFR_sfms(mstar, z_bin), 
                ls='--', 
                lw='3', 
                color='k')
        plt.plot(mstar, 
                AverageLogSFR_sfms(mstar, z_bin) + ScatterLogSFR_sfms(mstar, z_bin), 
                ls='--', 
                lw='3', 
                color='k')
        plt.text(9.25, -4.0, r"$\mathtt{z = "+str(round(z_bin,2))+"}$", fontsize=25)
        plt.legend(loc='lower right')

        fig_file = ''.join([
            'figure/test/', 
            'SFRevol.ModelSFMS', 
            '.', observable, '_z', str(round(z_bin, 2)), 
            '.png'])
        sfms_plot.save_fig(fig_file)
        plt.close()

def Plot_Fq_timedelay(lit='wetzelsmooth'): 
    '''
    '''
    prettyplot()
    pretty_colors = prettycolors()
    M_arr = np.arange(9.0, 12.5, 0.5)
    t_arr = np.arange(8.0, 13.8, 0.2)
    
    fig = plt.figure()
    sub = fig.add_subplot(111)
    for Mstar in M_arr:  
        time_delay = [] 
        print Mstar, getTauQ(np.array([Mstar]), tau_prop={'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6})
        print DeltaLogSFR_SF_Q_peak(np.array([Mstar]), Util.get_zsnap(np.mean(t_arr)), sfms_prop={'name': 'linear', 'mslope': 0.55, 'zslope': 1.1})
        for tt in t_arr: 
            z_in = Util.get_zsnap(tt)            
            #time_delay.append(
            #        DeltaLogSFR_SF_Q_peak(np.array([Mstar]), z_in, sfms_prop={'name': 'linear', 'mslope': 0.55, 'zslope': 1.1})
            #        )
            time_delay.append(
                    (DeltaLogSFR_SF_Q_peak(np.array([Mstar]), z_in, sfms_prop={'name': 'linear', 'mslope': 0.55, 'zslope': 1.1})*\
                            getTauQ(np.array([Mstar]), tau_prop={'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}))/0.434
                    )
        sub.plot(t_arr, time_delay, label=r"$\mathtt{M}_* = "+str(round(Mstar,2))+"$") 
    sub.set_xlim([t_arr.min(), t_arr.max()])
    sub.legend(loc='upper left')
    #sub.set_ylim([0., 1.])
    fig_file = ''.join(['figure/test/','fq_timedelay.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    return None


if __name__=='__main__': 
    Plot_Fq_timedelay(lit='wetzelsmooth')
    #ModelSFMS(observable='sdssprimus')
    #ModelSFMS(observable='groupcat')
