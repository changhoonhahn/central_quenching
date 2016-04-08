'''

Plotting stellar mass to halo mass ratio

'''
import numpy as np
import matplotlib.pyplot as plt

#----- Local -----
from plots import Plots

import bovy_plot as bovy
from util import cenque_utility as util

class PlotMstarMhalo(Plots): 

    def __init__(self, **kwargs): 
        ''' Class that describes the stellar mass to halo mass ratio plots 
        for CenQue objects
        '''

        super(PlotMstarMhalo, self).__init__(**kwargs)

    def cenque(self, cenque, **mkwargs):
        ''' 
        Plot stellar mass to halo mass ratio of input CenQue object
        '''

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            label = kwargs['label']
        else: 
            label = 'z = '+str(cenque.zsnap) 
        
        if 'color' in kwargs: 
            color = kwargs['color']
        else: 
            try: 
                color = self.pretty_colors[cenque.nsnap]
            except TypeError: 
                color = 'black'

        bovy.scatterplot(
                cenque.halo_mass, 
                cenque.mass, 
                scatter=True, 
                color=color, 
                s=3, 
                xrange=[10.0, 15.0],
                yrange=[9.0, 12.0], 
                ylabel='\mathtt{M_*}', 
                xlabel='\mathtt{M_{Halo}}'
                )

        return None   

    def bloodline(self, bloodline, snapshot, **mkwargs):
        '''
        '''
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(snapshot))

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            label = kwargs['label']
        else: 
            label = ''.join([
                'z =', 
                str(round(util.get_z_nsnap(descendant.nsnap), 2))
                ])
            print label
        
        if 'color' in kwargs: 
            color = kwargs['color']
        else: 
            try: 
                color = self.pretty_colors[descendant.nsnap]
            except TypeError: 
                color = 'black'

        bovy.scatterplot(
                descendant.halo_mass, 
                descendant.mass, 
                scatter=True, 
                color=color, 
                s=3, 
                xrange=[10.0, 15.0],
                yrange=[9.0, 12.0], 
                ylabel='\mathtt{M_*}', 
                xlabel='\mathtt{M_{Halo}}'
                )

        return None   

    def central_subhalos(self, centsub, **mkwargs): 
        '''
        Plot the stellar mass to halo mass relation for CentralSubhalos object
        '''

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            label = kwargs['label']
        else: 
            label = ''.join([
                'z =', 
                str(round(util.get_z_nsnap(centsub.snapshot), 2))
                ])
            print label
        
        if 'color' in kwargs: 
            color = kwargs['color']
        else: 
            try: 
                color = self.pretty_colors[centsub.snapshot]
            except TypeError: 
                color = 'black'

        halo_mass = getattr(centsub, 'halo.m.max')
        
        relevant_mstar_range = np.where(centsub.mass > 9.0)
        #print len(halo_mass[relevant_mstar_range])

        bovy.scatterplot(
                halo_mass[relevant_mstar_range], 
                centsub.mass[relevant_mstar_range], 
                scatter=True, 
                color=color, 
                s=3, 
                xrange=[10.0, 15.0],
                yrange=[9.0, 12.0], 
                ylabel='\mathtt{M_*}', 
                xlabel='\mathtt{M_{Halo}}',
                levels=[0.68, 0.95]
                )
        #self.subs.scatter(
        #        halo_mass, 
        #        centsub.mass, 
        #        color = color, 
        #        label = label
        #        ) 

        return None   

    def set_axes(self): 
        ''' 
        Set up axes
        '''
        self.subs.set_ylim([9.0, 12.0])
        self.subs.set_xlim([10.0, 15.0])
        
        self.subs.set_ylabel(r'Stellar Mass $\mathtt{M_*}$', fontsize=25) 
        self.subs.set_xlabel(r'Halo Mass $\mathtt{M_{Halo}}$', fontsize=25) 

        self.subs.legend(loc='upper left', frameon=False, scatterpoints=1)

        return None

    def save_fig(self, file_name): 
        '''
        Save figure to file_name
        '''

        plt.savefig(file_name, bbox_inches="tight")
        
        return None

def plot_mstar_mhalo_scatter_central_subhalo(): 
    for scat in [0.0, 0.2]: 
        for i_snap in [1,3,5,7,9]: 
            censub = CentralSubhalos(scatter = scat)
            censub.read(i_snap)
            mstar_mhalo_plot = PlotMstarMhalo()
            mstar_mhalo_plot.central_subhalos(censub)

            #mstar_mhalo_plot.set_axes()
            fig_file = ''.join([
                '/home/users/hahn/research/pro/tinker/central_quenching/figure/', 
                'plot_central_subhalo_mstar_mhalo_', 
                str(round(scat,1)), 
                'scatter_snapshot', 
                str(i_snap),
                '.png'
                ])
            print fig_file

            plt.savefig(fig_file, bbox_inches="tight")
            plt.close()

if __name__=="__main__":
    plot_mstar_mhalo_scatter_central_subhalo()
