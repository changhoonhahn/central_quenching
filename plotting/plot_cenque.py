'''

Plot CenQue class object


'''
import numpy as np
import matplotlib.pyplot as plt

# --- Local --- 
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

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
            cenque_ssfr_dist(self, cenque)
    
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
                    ssfr_hist_label = r'$\mathtt{z='+str(cenque.zsnap)+'}$' 
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

        return None

def plot_cenque_ssfr_dist(cenque, fig=None, **kwargs): 
    ''' Plot sSFR distribution for CenQue data
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

def plot_cenque_ssfr_dist_evolution(n_snaps = [12,11,10,9,8,7,6,5,4,3,2,1], Mrcut=18, **kwargs): 
    ''' Plot evolution of the CenQue SSFR distribution 

    Parameters
    ----------
    Mrcut : Absolute magnitude cut that specifies the group catalog 
    nsnaps : List of snapshot #s to plot  


    '''
    snap = cq.CenQue(n_snap=13, cenque_type = 'sf_assgined') 
    snap.readin()  

    ssfr_fig = PlotCenque(cenque=snap, lw=2, linestyle='--')
    
    # Overplot CenQue of specified Snapshots 
    for i_nsnap in n_snaps:  
        next_snap = cq.CenQue(n_snap = i_nsnap, cenque_type = 'evol_from13') 
        next_snap.readin()
        
        ssfr_fig.cenque_ssfr_dist(next_snap)
    
    # overplot SDSS group catalog sSFR dist
    central_ssfr = cq_group.central_catalog(Mrcut=Mrcut) 
    #ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, label= r'$M_\mathtt{r,cut} = '+str(Mrcut)+'$', **kwargs) 
    ssfr_fig = plot_cenque_ssfr_dist(central_ssfr, fig=ssfr_fig, label= r'SDSS Group Catalog', **kwargs) 
    """ 
    # file name ----------------------------------------------------------------------------
    if 'sfms_slope' in kwargs.keys():       # sfms specifier
        slope_str = str("%.2f" % kwargs['sfms_slope']) 
        yint_str = str("%.2f" % kwargs['sfms_yint']) 
        sfms_str = '_sfms_slope'+slope_str+'_yint'+yint_str
    else: 
        sfms_str = ''

    # Quenching Fraction specifier 
    if 'fqing_slope' in kwargs.keys(): 
        fqing_slope_str = str(kwargs['fqing_slope'])
    else: 
        fqing_slope_str = str(0.63)

    if 'fqing_yint' in kwargs.keys(): 
        fqing_yint_str = str(kwargs['fqing_yint'])
    else: 
        fqing_yint_str = str(-6.04) 

    fqing_str = ''.join([fqing_slope_str, '_', fqing_yint_str, 'fqing']) 

    # tau specifier
    if kwargs['tau'] == 'discrete': 
        tau_str = '_'.join( [str("%.1f" % t) for t in kwargs['tau_param']] )+'tau'
    elif kwargs['tau'] == 'linefit':
        tau_str = '_'.join( [str("%.2f" % t) for t in kwargs['tau_param']] )+'tau'
    else: 
        tau_str = kwargs['tau']+'tau'
    
    # Stellar mass specifier 
    if kwargs['stellmass'].lower() == 'integrated': 
        mass_str = '_integ'
    elif kwargs['stellmass'].lower() == 'sham': 
        mass_str = '_sham'
    else: 
        raise NotImplementedError('asdfalkjlkjasdf') 

    # SFR specifier
    if kwargs['sfr'] == 'sfr_avg': 
        file_type_str = mass_str+'_sfravg'
    elif kwargs['sfr'] == 'sfr_func': 
        file_type_str = mass_str+'_sfrfunc'
    else: 
        raise NotImplementedError('asdfasdflkjasdf;lkjasdf') 

    fig_file = ''.join(['figure/', 
        'cenque_ssfr_evol_', tau_str, '_', fqing_str, '_Mrcut', str(Mrcut), 
        file_type_str, '.png'])
    ssfr_fig.savefig(fig_file, bbox_inches='tight') 
    ssfr_fig.clear() 
    plt.close(ssfr_fig)
    """
