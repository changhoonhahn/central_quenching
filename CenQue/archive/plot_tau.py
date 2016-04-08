'''

Plotting class object for ploting quenching timescale 

'''

import numpy as np

from plots import Plots
from util.tau_quenching import get_quenching_efold

class PlotTau(Plots): 

    def __init__(self, tau_params, **kwargs):

        super(PlotTau, self).__init__(**kwargs)

        mass_bin = np.arange(9.0, 12.0, 0.1)
    
        tau_str = ''
        for i_tau, tau_param in enumerate(tau_params): 

            tau_mass = get_quenching_efold(mass_bin, tau_param = tau_param) 
            
            if tau_param['name'] == 'line': 
                tau_name = "Hahn+ (in prep)"

                #"Best-fit Central"
                #tau_name = ''.join([
                #    tau_param['name'], 
                #    '_slope', str(tau_param['slope']), 
                #    '_yint', str(tau_param['yint'])
                #    ])
                #tau_str += tau_name
                tau_str += tau_name
            else:
                tau_str += tau_param['name']

            self.subs.plot(
                    10**mass_bin, 
                    tau_mass, 
                    color=pretty_colors[i_tau+1], 
                    lw=4, 
                    label=tau_str
                    ) 
        
        # satellite quenching fraction for comparison 
        tau_mass = get_quenching_efold(mass_bin, tau_param = {'name': 'satellite'}) 

        self.subs.plot(10**mass_bin, tau_mass, color='black', lw=4, ls='--',
                label=r'Satellite (Wetzel+ 2013)') 

        self.subs.set_xlim([10**9.5, 10**11.75])
        self.subs.set_ylim([0.0, 2.0])
        self.subs.set_xscale('log')

        self.subs.set_xlabel(r'Stellar Mass $[\mathtt{M_\odot}]}$',fontsize=28) 
        self.subs.legend(loc='upper right') 

        self.subs.set_ylabel(r'Quenching Timescales $(\tau_\mathtt{Q})$ [Gyr]',fontsize=28) 

        self.fig_name = ''.join([
            '/home/users/hahn/research/pro/tinker/central_quenching/', 
            'figure/tau_quenching_efold', tau_str, '.png'
            ])

if __name__=="__main__":
    tau_params = [
            {'name': 'line', 'fid_mass': 10.75, 'slope': -0.57, 'yint': 0.5},
            {'name': 'line', 'fid_mass': 10.75, 'slope': -1.0, 'yint': 0.5}
            ]
    plot_tau(tau_params)
