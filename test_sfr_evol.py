'''

Test the SFR evolution functions 


'''
import numpy as np

from util import mpfit
import util.util as Util
from observations import GroupCat
from observations import PrimusSDSS

from sfr_evol import getTauQ
from sfr_evol import AverageLogSFR_sfms
from sfr_evol import ScatterLogSFR_sfms
from sfr_evol import DeltaLogSFR_SF_Q_peak 

from plotting.plots import PlotSFMS

import matplotlib.pyplot as plt
from ChangTools.utility.plotting import prettyplot
from ChangTools.utility.plotting import prettycolors 

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


if __name__=='__main__': 
    Plot_Fq_timedelay(lit='wetzelsmooth')
    #ModelSFMS(observable='sdssprimus')
    #ModelSFMS(observable='groupcat')
