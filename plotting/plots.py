'''

Plotting class objects

'''
import numpy as np 
import matplotlib.pyplot as plt

# SMF
from smf import SMF
# M*-Mhalo plot
import bovy_plot as bovy
from util import cenque_utility as util
# SF-MS plot
from util.cenque_utility import get_z_nsnap
from sfms.fitting import get_param_sfr_mstar_z
# Tau 
from util.tau_quenching import get_quenching_efold

from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

class Plots(object): 
    def __init__(self, **kwargs): 
        '''
        Class that generally encompasses all of the plotting I will do for the CenQue project
        '''
        self.kwargs = kwargs

        prettyplot()
        self.pretty_colors = prettycolors()

        self.fig = plt.figure(figsize=[10, 10])
        self.sub = self.fig.add_subplot(1,1,1)

    def save_fig(self, file_name): 
        '''
        save figure  to file_name
        '''
        self.fig.savefig(file_name, bbox_inches = 'tight')
        return None

    def show(self): 
        plt.show()
        plt.close()
        return None


# Quiescent Fraction 
class PlotFq(Plots):            
    def __init__(self, **kwargs): 
        ''' 
        Child class of Plots class that plots the quiescent fraction 
        for CenQue/GroupCat objects
        '''
       
        super(PlotFq, self).__init__(**kwargs) 

        self.z = None
        
        mass_bin = np.arange(9.0, 12.0, 0.2) 
        self.masses = 0.5 * (mass_bin[:-1] + mass_bin[1:])   # default masses
        
        self.fig = plt.figure(figsize=[10,10])
        self.subs = self.fig.add_subplot(1,1,1)
    def cenque(self, cenque, **mkwargs):
        ''' 
        Plot 'observed' fQ for CenQue data. Observed fQ is 
        computed using cq_fq function, which calculates fQ based
        on an evolving sSFR(M*,z) cut. 
        '''

        masses, fq = cenque.Fq()

        self.masses = masses
        self.z = cenque.zsnap

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            fq_label = kwargs['label']
        else: 
            fq_label = 'z ='+str(cenque.zsnap) 
        
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

        self.subs.plot(
                masses, 
                fq, 
                color = line_color, 
                lw = line_width, 
                ls = line_style, 
                label = fq_label
                ) 

        return None   
    def param_fq(self, fq_prop = {'name': 'wetzelsmooth'}, z = None, **mkwargs):
        ''' 
        Parameterized queiscent fraction  
        '''
        if z is None: 
            if self.z is None: 
                raise ValeuError
            else: 
                redshift = self.z
        else: 
            redshift = z
            self.z = z 

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            fq_label = kwargs['label']
        else: 
            fq_label = fq_prop['name']+'; z = '+str(redshift) 
        
        if 'line_color' in kwargs: 
            line_color = kwargs['line_color']
        else: 
            line_color = 'black'

        if 'line_style' in kwargs:
            line_style = kwargs['line_style'] 
        else: 
            line_style = '-'

        if 'lw' in kwargs: 
            line_width = kwargs['lw'] 
        else:
            line_width = 4

        # parameterized fq 
        fq = get_fq(
                self.masses, 
                redshift, 
                lit = fq_prop['name']
                )

        self.subs.plot(
                self.masses, fq, 
                color = line_color,  
                lw = 4, 
                ls = '--', 
                label = fq_label 
                ) 

        return None
    def set_axes(self): 
        ''' Set up axes
        '''
        self.subs.set_xlim([9.0, 12.0])
        self.subs.set_ylim([0.0, 1.0])
        
        self.subs.set_xlabel(r'Mass $\mathtt{M_*}$') 
        self.subs.set_ylabel(r'Quiescent Fraction $\mathtt{f_Q}$', fontsize=20) 

        self.subs.legend(loc='upper left', frameon=False)

        return None


# M* - Mahlo plots
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


# Star Forming Main Sequence
class PlotSFMS(Plots): 
    def __init__(self, **kwargs): 
        ''' 
        Child class of Plots class that plots the StarForming Main Sequence for 
        different class objects (CenQue/GroupCat)
        '''
        super(PlotSFMS, self).__init__(**kwargs)
        self.cenque_nsnap = None
        self.scatter = None

    def cenque(self, cq_obj, scatter=False, bovyplot=True, **mkwargs): 
        '''
        Plot SF-MS (M* vs SFR) for CenQue object
        '''
        if self.kwargs == {}: 
            self.kwargs = mkwargs.copy() 
        else: 
            self.kwargs = (self.kwargs).copy()
            self.kwargs.update(mkwargs)
    
        if 'label' in self.kwargs: 
            label = self.kwargs['label']
        else: 
            label = 'z = '+str(cq_obj.zsnap) 
        
        if 'color' in self.kwargs: 
            color = self.kwargs['color']
        else: 
            try: 
                color = self.pretty_colors[cq_obj.nsnap]
                self.cenque_nsnap = cq_obj.nsnap
            except TypeError: 
                color = 'black'

        if 'justsf' in self.kwargs: 
            if self.kwargs['justsf']: 
                justsf = np.where(cq_obj.gal_type == 'star-forming')

                cq_mass = cq_obj.mass[justsf]
                cq_sfr = cq_obj.sfr[justsf]
            else: 
                cq_mass = cq_obj.mass
                cq_sfr = cq_obj.sfr
        else: 
            cq_mass = cq_obj.mass
            cq_sfr = cq_obj.sfr
        
        if scatter: 
            self.sub.scatter(
                    cq_mass,
                    cq_sfr,
                    color=color, 
                    s=3
                    )
            self.sub.set_xlim([9.0, 12.0])
            self.sub.set_ylim([-5.0, 2.0])
            self.sub.set_xlabel(r'$\mathtt{M_*}$')
            self.sub.set_ylabel(r'$\mathtt{log\;SFR}$')
            self.scatter = True

        else: 
            bovy.scatterplot(
                    cq_mass, 
                    cq_sfr,
                    scatter=True, 
                    color=color, 
                    s=3, 
                    xrange=[9.0, 12.0], 
                    yrange=[-5.0, 2.0], 
                    xlabel='\mathtt{M_*}', 
                    ylabel='\mathtt{log\;SFR}'
                    )
            self.scatter = False

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
                descendant.mass, 
                descendant.sfr, 
                scatter=True, 
                color=color, 
                s=3, 
                xrange=[9.0, 12.0], 
                yrange=[-3.0, 3.0], 
                xlabel='\mathtt{M_*}', 
                ylabel='\mathtt{log\;SFR}'
                )

        return None   

    def param_sfms(self, nsnap = None, **pltkwargs): 
        '''
        Plot parameterized SFMS
        '''
        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

        if nsnap is None: 
            if self.cenque_nsnap is not None: 
                nsnap = self.cenque_nsnap
            else: 
                raise ValueError('nsnap needs to be specified')
        
        if 'justsf' in self.kwargs: 
            self.sub.plot(
                np.arange(8.5, 12.0, 0.1), 
                sfr_mstar_z(np.arange(8.5, 12.0, 0.1), get_z_nsnap(nsnap)), 
                c = 'k', 
                ls = '--', 
                lw = 3 
                )
        else:
            plt.plot(
                    np.arange(8.5, 12.0, 0.1), 
                    sfr_mstar_z(np.arange(8.5, 12.0, 0.1), get_z_nsnap(nsnap)), 
                    c = 'k', 
                    ls = '--', 
                    lw = 3 
                    )

        return None
        
    def save_fig(self, file_name): 
        '''
        Save figure to file 
        '''
        if 'justsf' in self.kwargs: 
            if self.kwargs['justsf']: 
                if self.scatter: 
                    return super(PlotSFMS, self).save_fig(file_name) 
                else: 
                    plt.savefig(file_name, bbox_inches='tight') 
        else: 
            return super(PlotSFMS, self).save_fig(file_name) 
        return None


# Quenching e-folding time
class PlotTau(Plots): 
    def __init__(self, tau_params, **kwargs):
        '''
        '''
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


# Specific Star Formation Rates
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
        '''
        Mass bin panels of the figure
        '''
        panel_mass_bins = [ 
                [9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]
                ]
        return panel_mass_bins

    def set_axes(self): 
        ''' Set up axes
        '''
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


# Stellar Mass Function 
class PlotSMF(Plots): 
    def __init__(self, **kwargs): 
        ''' 
        Child class of Plots class that describes SMF plots for 
        different class objects including CenQue and GroupCat
        '''
        super(PlotSMF, self).__init__(**kwargs)

    def cenque(self, cenque, type='total', **kwargs): 
        ''' 
        Plot SMF distribution for CenQue class object
        
        Parameters
        ----------
        cenque : 
            Central Quenching Class Object
        type : 
            type of SMF plot. If type = 'total', then plot the total SMF. 
            if tyep='sfq' then plot the SF and Q components separately in the
            SMF. 
        '''
        if 'label' in kwargs: 
            smf_label = kwargs['label']
            kwargs.pop('label', None)
        else: 
            smf_label = 'z ='+str(cenque.zsnap) 
    
        if 'line_color' in kwargs: 
            line_color = kwargs['line_color']
            kwargs.pop('line_color', None)
        else: 
            try: 
                line_color = self.pretty_colors[cenque.nsnap]
            except TypeError: 
                line_color = 'black'

        if 'line_style' in kwargs:
            line_style = kwargs['line_style'] 
            kwargs.pop('line_style', None)
        else: 
            line_style = '-'

        if 'lw' in kwargs: 
            line_width = kwargs['lw'] 
            kwargs.pop('lw', None)
        else:
            line_width = 4

        #mf = SMF()
        mass, phi = cenque.SMF(**kwargs) #mf.cenque(cenque, **mkwargs)

        self.sub.plot(mass, phi, 
                c=line_color, ls=line_style, lw=line_width,
                label=smf_label)
        return mass, phi 
    
    def analytic(self, redshift, source='li-drory-march', **kwargs): 
        '''
        Analytic SMF at redshift 
        '''
        if 'label' in kwargs: 
            smf_label = kwargs['label']
            kwargs.pop('label', None)
        else: 
            smf_label = 'Analytic'
    
        if 'line_color' in kwargs: 
            line_color = kwargs['line_color']
            kwargs.pop('line_color', None)
        else: 
            line_color = 'black'

        if 'line_style' in kwargs:
            line_style = kwargs['line_style'] 
            kwargs.pop('line_style', None)
        else: 
            line_style = '--'

        if 'lw' in kwargs: 
            line_width = kwargs['lw'] 
            kwargs.pop('lw', None)
        else:
            line_width = 3
        smf = SMF()
        mass, phi = smf.analytic(redshift, source=source) 
        
        self.sub.plot(mass, phi, 
                c=line_color, ls=line_style, lw=line_width,
                label=smf_label
                )

    def groupcat(self, Mrcut=18, **kwargs): 
        ''' 
        Plot sSFR distribution for Group Catalog data
        '''
        raise NotImplementedError

    def set_axes(self): 
        ''' 
        Set up axes of figures
        '''
        self.sub.set_ylim([10**-5, 10**-1])
        self.sub.set_xlim([7.5, 12.0])
        self.sub.set_yscale('log')

        # x and y labels
        self.sub.set_xlabel(r'Mass $\mathtt{M_*}$') 
        self.sub.set_ylabel(r'Stellar Mass Function $\mathtt{\Phi}$', fontsize=20) 
        
        self.sub.legend(loc='upper left', frameon=False)
        return None
