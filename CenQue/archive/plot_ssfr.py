'''

Plot CenQue class object


'''
import numpy as np
import matplotlib.pyplot as plt

# --- Local --- 
from plots import Plots
from ssfr import Ssfr
from util.cenque_utility import get_q_ssfr_mean
from sfms.fitting import get_param_sfr_mstar_z
from group_catalog.group_catalog import central_catalog

class PlotSSFR(Plots): 

    def __init__(self, cenque=None, **kwargs): 
        ''' 
        Child class of Plots class that describes the sSFR distribution plots for 
        CenQue and GroupCat objects 
        '''
        super(PlotSSFR, self).__init__(**kwargs)

        self.fig = plt.figure(figsize=(16,16))
        self.fig.subplots_adjust(hspace=0., wspace=0.)
    
        self.panel_mass_bins = self.mass_panels()
        self.n_mass_bins = len(self.panel_mass_bins)

        # panel subplot
        self.subs = [
                self.fig.add_subplot(2, 2, i_mass+1) 
                for i_mass in xrange(self.n_mass_bins)
                ]  

    def cenque(self, cenque, quenching=False, **mkwargs): 
        ''' 
        Plot sSFR distribution for CenQue class object
        
        Parameters
        ----------
        cenque : 
            Central Quenching Class Object
        quenching : 
            True : Plot stacked sSFR distribution for CenQue data that 
                differentiates the quiescent, quenching, and star-forming populations
            False : Plot SSFR distribution for the entire CenQue data

        '''
        
        if not quenching: 

            ssfr_bin_mid, ssfr_hist = cenque.Ssfr()  # CenQue SSFR distribution 

            if self.kwargs == {}: 
                kwargs = mkwargs.copy() 
            else: 
                kwargs = (self.kwargs).copy()
                kwargs.update(mkwargs)
        
            for i_mass, panel_mass in enumerate(self.panel_mass_bins):       # loop through each panel 

                if 'label' in kwargs: 
                    ssfr_hist_label = kwargs['label']
                else: 
                    ssfr_hist_label = 'z ='+str(cenque.zsnap) 
            
                if 'line_color' in kwargs: 
                    line_color = kwargs['line_color']
                else: 
                    try: 
                        line_color = self.pretty_colors[cenque.nsnap]
                    except TypeError: 
                        line_color = 'black'

                if 'line_style' in kwargs:
                    line_style = kwargs['line_style'] 
                else: 
                    line_style = '-'

                if 'lw' in kwargs: 
                    line_width = kwargs['lw'] 
                else:
                    line_width = 4

                self.subs[i_mass].plot(
                        ssfr_bin_mid[i_mass], 
                        ssfr_hist[i_mass], 
                        color = line_color, 
                        lw = line_width, 
                        ls = line_style, 
                        label = ssfr_hist_label
                        ) 
        else: 
            for i_mass, panel_mass in enumerate(self.panel_mass_bins):       # loop through each panel 
                
                sf_massbin = np.where(
                        (cenque.mass >= panel_mass[0]) & 
                        (cenque.mass < panel_mass[1]) & 
                        (cenque.gal_type == 'star-forming') 
                        )

                quenching_massbin = np.where(
                        (cenque.mass >= panel_mass[0]) & 
                        (cenque.mass < panel_mass[1]) & 
                        (cenque.tau > 0.0) 
                        )

                q_notquenching_massbin = np.where(
                        (cenque.mass >= panel_mass[0]) & 
                        (cenque.mass < panel_mass[1]) & 
                        (cenque.gal_type == 'quiescent') &
                        (cenque.tau < 0.0) 
                        )

                data_list, color_list, label_list = [], [], []

                if len(quenching_massbin[0]) > 0: 
                    data_list.append(
                            cenque.ssfr[quenching_massbin]
                            )
                    color_list.append('orange')
                    label_list.append('Quenching')

                if len(sf_massbin[0]) > 0: 
                    data_list.append(
                            cenque.ssfr[sf_massbin]
                            )
                    color_list.append('blue')
                    label_list.append('Star-Forming')

                if len(q_notquenching_massbin[0]) > 0:
                    data_list.append(
                            cenque.ssfr[q_notquenching_massbin]
                            )
                    color_list.append('red')
                    label_list.append('Quiescent')
                #print 'quenching ', len(quenching_massbin[0])
                #print 'starforming ', len(sf_massbin[0])
                #print 'quiescent ', len(q_notquenching_massbin[0])

                self.subs[i_mass].hist(
                        data_list, 
                        50, 
                        normed = 1, 
                        stacked=True, 
                        color = color_list,
                        label = label_list 
                        )

                if panel_mass == self.panel_mass_bins[-1]: 
                    self.subs[i_mass].text(-9.25, 1., 'z ='+str(cenque.zsnap),
                        fontsize=24
                        )

        return None   

    def groupcat(self, Mrcut=18, **kwargs): 
        ''' 
        Plot sSFR distribution for Group Catalog data
        '''
    
        ssfr_dist = Ssfr()
        ssfr_bin_mid, ssfr_hist = ssfr_dist.groupcat(Mrcut=Mrcut)
    
        for i_mass, panel_mass in enumerate(self.panel_mass_bins):       # loop through each panel 

            if Mrcut == 18: 
                z_med = 0.03
            elif Mrcut == 19: 
                z_med = 0.05
            elif Mrcut == 20:
                z_med = 0.08

            ssfr_label= ''.join([
                r'SDSS $\mathtt{M_r =', str(Mrcut), '}$'])
            #, ', z = ', str(z_med)
            #    ]) 

            if 'lw' in kwargs.keys(): 
                lwid = kwargs['lw']
            else: 
                lwid = 4

            if 'ls' in kwargs.keys(): 
                lsty = kwargs['ls']
            else: 
                lsty = '--'

            if 'color' in kwargs.keys(): 
                col = kwargs['color']
            else: 
                col = 'k' 

            self.subs[i_mass].plot(
                    ssfr_bin_mid[i_mass],
                    ssfr_hist[i_mass], 
                    color = col, 
                    lw = lwid, 
                    ls = lsty, 
                    label = ssfr_label) 

        return None   
    
    def mass_panels(self): 
        """ Mass bin panels of the figure
        """
        panel_mass_bins = [ 
                [9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]
                ]
        return panel_mass_bins

    def set_axes(self): 
        """ Set up axes
        """
        for i_sub in xrange(len(self.subs)): 
            self.subs[i_sub].set_xlim([-13.0, -7.0])
            self.subs[i_sub].set_ylim([0.0, 1.6])
            
            massbin_str = ''.join([ 
                r'$\mathtt{log \; M_{*} = [', 
                str(self.panel_mass_bins[i_sub][0]), ',\;', 
                str(self.panel_mass_bins[i_sub][1]), ']}$'
                ])
            self.subs[i_sub].text(-10.5, 1.4, massbin_str,
                    fontsize=24
                    )

            if i_sub == 0: 
                self.subs[i_sub].set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=20) 
                self.subs[i_sub].set_xticklabels([])
            elif i_sub == 1: 
                self.subs[i_sub].set_xticklabels([])
                self.subs[i_sub].set_yticklabels([])
            elif i_sub == 2:
                #sub.set_yticklabels([])
                self.subs[i_sub].set_ylabel(r'$\mathtt{P(log \; SSFR)}$', fontsize=20) 
                self.subs[i_sub].set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=20) 
            else: 
                self.subs[i_sub].set_yticklabels([])
                self.subs[i_sub].set_xlabel(r'$\mathtt{log \; SSFR \;[yr^{-1}]}$', fontsize=20) 

                self.subs[i_sub].legend(loc='lower right', frameon=False)

        return None
    
"""
def plot_cenque_ssfr_dist(cenque, fig=None, **kwargs): 
    ''' wrapper to plot sSFR distribution for CenQue data
    '''
    if fig == None: 
        pltcq = PlotCenque()
    else: 
        pltcq = fig 

    pltcq.cenque_ssfr_dist(cenque)
    
    for i_mass, panel_mass in enumerate(pltcq.panel_mass_bins):       # loop through each panel 

        ssfr_cut = -11.35 + 0.76*(cenque.zsnap-0.05) - 0.35*((0.5 * np.sum(panel_mass))-10.5)

        pltcq.subs[i_mass].vlines(ssfr_cut, 0.0, 10.0, lw=4, ls='--')

    pltcq.set_axes()

    return pltcq   

def plot_ssfr_cenque(n_snaps = [1], original_nsnap = None, cenque_type = 'evol_from13', 
        fq_prop = {'name': 'wetzelsmooth'}, 
        sf_prop = {'name': 'average'}, 
        tau_prop = {'name': 'satellite'}, 
        Mrcut = 18): 
    ''' Plot evolved CenQue galaxy catalog with the group catalog for comparison 
    '''

    prettyplot()
    
    if original_nsnap != None: 
        # plot original nsnap or not 
        snap = CenQue(n_snap = original_nsnap, cenque_type = 'sf_assigned') 
        snap.readin()

        ssfr_fig = PlotCenque(cenque = snap, lw = 2, linestyle = '--') 
    
    for i_nsnap in n_snaps: 
        next_snap = CenQue(n_snap = i_nsnap, cenque_type = cenque_type) 
        next_snap.tau_prop = {'name': 'satellite'} 
        next_snap.readin()
        
        try: 
            ssfr_fig.cenque_ssfr_dist(next_snap, label=r'Satellite $\tau$')
        except UnboundLocalError: 
            ssfr_fig = PlotCenque(cenque=next_snap, label=r'Satellite $\tau$')

    # plot specified snaphsots
    for i_nsnap in n_snaps: 
        next_snap = CenQue(n_snap = i_nsnap, cenque_type = cenque_type) 
        next_snap.tau_prop = tau_prop
        next_snap.readin()
        
        try: 
            ssfr_fig.cenque_ssfr_dist(next_snap, line_color='red', label=r'Central $\tau$')
        except UnboundLocalError: 
            ssfr_fig = PlotCenque(cenque=next_snap, line_color='red', label=r'Central $\tau$')
    
    # plot SDSS group catalog
    if Mrcut != None: 
        ssfr_fig.groupcat_ssfr_dist(Mrcut=Mrcut)

        if Mrcut == 18: 
            z_med = 0.03
        elif Mrcut == 19: 
            z_med = 0.05
        else: 
            raise NotImplementedError()
    
    for i_mass, panel_mass in enumerate(ssfr_fig.panel_mass_bins):      
        # sSFR cutoff for SF/Q galaxies 
        ssfr_cut = -11.35 + 0.76*(next_snap.zsnap-0.05) - 0.35*((0.5 * np.sum(panel_mass))-10.5)

        #ssfr_fig.subs[i_mass].vlines(ssfr_cut, 0.0, 10.0, lw=2, linestyles='dotted')
        
        # Green valley region (between the peaks of the two distributions) 
        q_ssfr_massbin = np.mean(get_q_ssfr_mean(panel_mass)) 

        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

        #sf_ssfr_massbin_min = sfr_mstar_z(panel_mass[0], z_med) - panel_mass[0]
        sf_ssfr_massbin_max = sfr_mstar_z(panel_mass[1], z_med) - panel_mass[1]
        
        #ssfr_fig.subs[i_mass].vlines(q_ssfr_massbin, 0.0, 10.0, lw=2, color='green')
        #ssfr_fig.subs[i_mass].vlines(sf_ssfr_massbin_min, 0.0, 10.0, lw=2, color='green')
        #ssfr_fig.subs[i_mass].vlines(sf_ssfr_massbin_max, 0.0, 10.0, lw=2, color='green')

    ssfr_fig.set_axes()
    
    if tau_prop['name'] == 'line': 
        tau_str = ''.join([
            tau_prop['name'], 
            '_slope', str(tau_prop['slope']), 
            '_yint', str(tau_prop['yint'])
            ])
    else: 
        tau_str = tau_prop['name']

    fig_file = ''.join([
        '/home/users/hahn/research/pro/tinker/central_quenching/figure/', 
        'plot_ssfr_nsnap',
        ''.join([str(i_nsnap) for i_nsnap in n_snaps]),
        '_', cenque_type, 
        '_', tau_str, 'tau', 
        '_Mr', str(Mrcut), 
        '.png'
        ]) 
        
    ssfr_fig.fig.savefig(fig_file) 
    ssfr_fig.fig.clear()
    plt.close()
"""
