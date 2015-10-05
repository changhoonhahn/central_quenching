'''

Plot CenQue class object


'''
import numpy as np
import matplotlib.pyplot as plt

# --- Local --- 
from cenque import CenQue
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from group_catalog.group_catalog import central_catalog

class PlotCenque: 

    def __init__(self, cenque=None, **kwargs): 
        """ Class that describes the sSFR distribution plots for 
        CenQue objects
        """
        self.kwargs = kwargs 

        self.fig = plt.figure(1, figsize=(16,16))
        self.fig.subplots_adjust(hspace=0., wspace=0.)

        self.panel_mass_bins = self.mass_panels()
        self.subs = [
                self.fig.add_subplot(2, 2, i_mass+1) 
                for i_mass in xrange(len(self.panel_mass_bins))
                ] # panel subplot 
        self.pretty_colors = prettycolors()

        if cenque != None: 
            self.cenque_ssfr_dist(cenque)
    
    def cenque_ssfr_dist(self, cenque): 
        ''' Plot sSFR distribution for CenQue data
        '''
    
        for i_mass, panel_mass in enumerate(self.panel_mass_bins):       # loop through each panel 

            mass_limit = np.where(
                    (cenque.mass >= panel_mass[0]) & 
                    (cenque.mass < panel_mass[1])
                    )
            ngal_bin = len(mass_limit[0])

            # SSFR histogram 
            ssfr_hist, ssfr_bin_edges = np.histogram(
                    cenque.ssfr[mass_limit], 
                    range = [-13.0, -7], 
                    bins = 40, 
                    normed = True
                    )
            ssfr_bin_low = ssfr_bin_edges[:-1]
            ssfr_bin_high = ssfr_bin_edges[1:]
            ssfr_bin_mid = 0.5 * (ssfr_bin_low + ssfr_bin_high) 
        
            if 'label' in self.kwargs: 
                ssfr_hist_label = self.kwargs['label']
            else: 
                try: 
                    ssfr_hist_label = 'z ='+str(cenque.zsnap) 
                except AttributeError: 
                    ssfr_hist_label = 'Centrals'
        
            if 'line_color' in self.kwargs: 
                line_color = self.kwargs['line_color']
            else: 
                try: 
                    line_color = self.pretty_colors[cenque.nsnap]
                except TypeError: 
                    line_color = 'black'

            if 'line_style' in self.kwargs:
                line_style = self.kwargs['line_style'] 
            else: 
                line_style = '-'

            if 'lw' in self.kwargs: 
                line_width = self.kwargs['lw'] 
            else:
                line_width = 4

            self.subs[i_mass].plot(
                    ssfr_bin_mid, 
                    ssfr_hist, 
                    color = line_color, 
                    lw = line_width, 
                    ls = line_style, 
                    label = ssfr_hist_label) 

        return None   

    def groupcat_ssfr_dist(self, Mrcut=18): 
        ''' Plot sSFR distribution for Group Catalog data
        '''
    
        groupcat = central_catalog(Mrcut=Mrcut)
    
        for i_mass, panel_mass in enumerate(self.panel_mass_bins):       # loop through each panel 

            mass_limit = np.where(
                    (groupcat.mass >= panel_mass[0]) & 
                    (groupcat.mass < panel_mass[1])
                    )
            ngal_bin = len(mass_limit[0])

            # SSFR histogram 
            ssfr_hist, ssfr_bin_edges = np.histogram(
                    groupcat.ssfr[mass_limit], 
                    range = [-13.0, -7], 
                    bins = 40, 
                    normed = True
                    )
            ssfr_bin_low = ssfr_bin_edges[:-1]
            ssfr_bin_high = ssfr_bin_edges[1:]
            ssfr_bin_mid = 0.5 * (ssfr_bin_low + ssfr_bin_high) 

            ssfr_label= 'SDSS Group Catalog Mr ='+str(Mrcut)
        
            self.subs[i_mass].plot(
                    ssfr_bin_mid, 
                    ssfr_hist, 
                    color = 'k', 
                    lw = 4, 
                    ls = '--', 
                    label = ssfr_label) 

        return None   
    
    def cenque_quenching_ssfr_dist(self, cenque): 
        """ Plot stacked sSFR distribution for CenQue data that 
        differentiates the quiescent, quenching, and star-forming populations
        """
    
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

            data_list = []  
            color_list = [] 
            label_list = [] 
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
