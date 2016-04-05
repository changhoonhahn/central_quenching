'''

Plotting class objects

'''
import numpy as np 
import matplotlib.pyplot as plt

from observations import GroupCat
from central_subhalo import CentralSubhalos

from gal_prop import Fq 
from gal_prop import SMF
from gal_prop import Ssfr
from sfr_evol import getTauQ
from sfr_evol import AverageLogSFR_sfms

import bovy_plot as bovy

# SMF
#from smf import SMF
#from quiescent_fraction import get_fq
#from quiescent_fraction import sfr_cut 
## M*-Mhalo plot
#from util import cenque_utility as util
## SF-MS plot
#from util.cenque_utility import get_z_nsnap
#from sfms.fitting import get_param_sfr_mstar_z
## Tau 
#from util.tau_quenching import get_quenching_efold

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
    
        Ssfr_obj = Ssfr()
        self.panel_mass_bins = Ssfr_obj.mass_bins
        self.n_mass_bins = len(self.panel_mass_bins)

        # panel subplot
        self.subs = [
                self.fig.add_subplot(2, 2, i_mass+1) 
                for i_mass in xrange(self.n_mass_bins)
                ]  

    def plot(self, mass=None, ssfr=None, SSFRdist=None, sfms_prop=None, z=None, **pltkwargs): 
        ''' Plot the sSFR distriubiton 
        '''
        if self.kwargs == {}: 
            kwargs = pltkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(pltkwargs)

        if SSFRdist is not None: 
            ssfr_bin_mid, ssfr_dist = SSFRdist
        else: 
            if mass is None or ssfr is None: 
                raise ValueError
            else: 
                Ssfr_obj = Ssfr()
                ssfr_bin_mid, ssfr_dist = Ssfr_obj.Calculate(mass, ssfr)
        
        # loop through mass bin 
        for i_mass, panel_mass in enumerate(self.panel_mass_bins):       
        
            # plot configuration
            if 'label' in kwargs:   # label 
                hist_label = kwargs['label']
            else: 
                hist_label = None 
            # line color 
            if 'line_color' in kwargs: 
                if isinstance(kwargs['line_color'], str): 
                    line_color = kwargs['line_color']
                elif isinstance(kwargs['line_color'], int): 
                    line_color = self.pretty_colors[kwargs['line_color']]
            else: 
                line_color = 'black'
            # line style
            if 'line_style' in kwargs:
                line_style = kwargs['line_style'] 
            else: 
                line_style = '-'
            # line width
            if 'lw' in kwargs: 
                line_width = kwargs['lw'] 
            else:
                line_width = 4

            self.subs[i_mass].plot(
                    ssfr_bin_mid[i_mass], 
                    ssfr_dist[i_mass], 
                    color=line_color, 
                    lw=line_width, 
                    ls=line_style, 
                    label=hist_label
                    ) 
            
            if sfms_prop is not None: 
                if z is None: 
                    raise ValueError
                qf = Fq()
                self.subs[i_mass].vlines(
                        qf.SFRcut(np.array([np.mean(panel_mass)]), z, 
                            sfms_prop=sfms_prop 
                            )-np.mean(panel_mass), 
                        0., 100., 
                        lw=3, linestyle='--', color='k')

            if i_mass == 2: 
                plt.sca(self.subs[i_mass])
                plt.xticks([-13, -12, -11, -10, -9, -8])
                plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

        return None   

    def GroupCat(self, Mrcut=18, **kwargs): 
        ''' Plot sSFR distribution for Group Catalog data
        '''
        groupcat = GroupCat(Mrcut=Mrcut, position='central')
        ssfr_bin_mid, ssfr_hist =groupcat.Ssfr()

        # loop through each panel 
        for i_mass, panel_mass in enumerate(self.panel_mass_bins):       
            if Mrcut == 18: 
                z_med = 0.03
            elif Mrcut == 19: 
                z_med = 0.05
            elif Mrcut == 20:
                z_med = 0.08

            ssfr_label= ''.join([
                r'SDSS $\mathtt{M_r =', str(Mrcut), '}$'])

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
    
    def set_axes(self): 
        ''' Set up axes for the SSFR distribution figure
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


# Quiescent Fraction 
class PlotFq(Plots):            
    def __init__(self, **kwargs): 
        ''' Child class of Plots class that plots the quiescent fraction 
        as a function of stellar mass.
        '''
        super(PlotFq, self).__init__(**kwargs) 

        self.z = None
        
        mass_bin = np.arange(9.0, 12.0, 0.2) 
        self.masses = 0.5 * (mass_bin[:-1] + mass_bin[1:])   # default masses
        
        self.fig = plt.figure(figsize=[10,10])
        self.subs = self.fig.add_subplot(1,1,1)

    def plot(self, mass=None, sfr=None, z=None, FQdist=None, **mkwargs):
        ''' Plot the quiescent fraction as a function of stellar mass. 
        The quiescent fraction is calculated based on an evolving sSFR(M*,z) cut. 

        Parameters
        ----------
        FQdist : list
            List that specifies the mass bins and the quiescent fraction 
            (eq. [massbin, quiescentfraction])
        '''
        if FQdist is None: 
            if sfr is None or z is None: 
                raise ValueError
            fq_obj = Fq()
            masses, fq = fq_obj.Calculate(mass=mass, sfr=sfr, z=z) 
        else:
            masses, fq = FQdist

        self.masses = masses
        self.z = z  

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            fq_label = kwargs['label']
        else: 
            fq_label = None
        
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

        self.subs.plot(masses, fq, 
                color = line_color, 
                lw = line_width, 
                ls = line_style, 
                label = fq_label) 

        return None   

    def model(self, fq_prop={'name': 'wetzelsmooth'}, z=None, **mkwargs):
        ''' Plot the model Parameterized queiscent fraction as a function 
        of stellar mass
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
    
        if 'label' in kwargs and kwargs['label'] is not None: 
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
        fq_obj = Fq()
        fq = fq_obj.model(self.masses, self.z, lit=fq_prop['name'])
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


# Quenching e-folding time
class PlotTau(Plots): 
    def __init__(self, **kwargs):
        ''' Plots subclass that plots the quenching e-fold time as a function 
        of stellar mass 
        '''
        super(PlotTau, self).__init__(**kwargs)
    
    def plot(self, tau_dicts, **kwargs):
        self.mass_bin = np.arange(9.0, 12.0, 0.1)

        if not isinstance(tau_dicts, list): 
            tau_dicts = [tau_dicts]
    
        tau_str = ''
        for i_tau, tau_param in enumerate(tau_dicts): 
            tau_m = getTauQ(self.mass_bin, tau_prop=tau_param) 
            
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

            self.sub.plot(
                    10**self.mass_bin, 
                    tau_m, 
                    color=self.pretty_colors[i_tau+1], 
                    lw=4, 
                    label=tau_str
                    ) 
        
        # satellite quenching fraction for comparison 
        tau_sat = getTauQ(self.mass_bin, tau_prop={'name': 'satellite'}) 

        self.sub.plot(10**self.mass_bin, tau_sat, color='black', lw=4, ls='--',
                label=r'Satellite (Wetzel+ 2013)') 

        self.sub.set_xlim([10**9.5, 10**11.75])
        self.sub.set_ylim([0.0, 2.0])
        self.sub.set_xscale('log')

        self.sub.set_xlabel(r'Stellar Mass $[\mathtt{M_\odot}]}$',fontsize=28) 
        self.sub.legend(loc='upper right') 

        self.sub.set_ylabel(r'Quenching Timescales $(\tau_\mathtt{Q})$ [Gyr]',fontsize=28) 

        self.fig_name = ''.join([
            '/home/users/hahn/research/pro/tinker/central_quenching/', 
            'figure/tau_quenching_efold', tau_str, '.png'
            ])
        return None 

        
# M* - Mahlo plots
class PlotSMHM(Plots): 
    def __init__(self, **kwargs): 
        ''' Class that describes the Stellar Mass to Halo Mass relation figures.
        '''
        super(PlotSMHM, self).__init__(**kwargs)

    def plot(self, stellarmass=None, halomass=None, bovyplot=True, **mkwargs):
        ''' Plot SMHM ratio of input stellar mass and halo mass
        '''
        if stellarmass is None or halomass is None: 
            raise ValueError

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            label = kwargs['label']
        else: 
            label = None 

        if 'color' in kwargs: 
            if isinstance(kwargs['color'], str):
                color = kwargs['color']
            elif isinstance(kwargs['color'], int): 
                color = self.pretty_colors[kwargs['color']]
        else: 
            color = 'black'
        
        if bovyplot: 
            bovy.scatterplot(
                    halomass, 
                    stellarmass,
                    scatter=True, 
                    color=color, 
                    s=3, 
                    xrange=[10.0, 15.0],
                    yrange=[9.0, 12.0], 
                    ylabel='\mathtt{M_*}', 
                    xlabel='\mathtt{M_{Halo}}'
                    )
        else: 
            self.sub.scatter(
                    halomass, 
                    stellarmass,
                    color=color,
                    s=3, 
                    label=label)
        return None   

    def plotSummary(self, stellarmass=None, halomass=None, **mkwargs):
        ''' Plot the summary statistis of the SMHM relation. The average M* 
        for halo mass bins along with the scatter. 
        '''
        if halomass is None or stellarmass is None: 
            raise ValueError
        mbin = np.arange(halomass.min(), halomass.max(), 0.25) 
        mlow = mbin[:-1]
        mhigh = mbin[1:]

        muMstar, sigMstar = [], [] 
        for im in range(len(mlow)): 
            mbin = np.where((halomass > mlow[im]) & (halomass <= mhigh[im]))
            muMstar.append(np.mean(stellarmass[mbin]))
            sigMstar.append(np.std(stellarmass[mbin]))
        self.sub.plot(0.5*(mlow+mhigh), muMstar, color='k', ls='--', lw=4)
        self.sub.errorbar(0.5*(mlow+mhigh), muMstar, yerr=sigMstar, color='k', lw=3, fmt='o', capthick=2)

    def set_axes(self): 
        ''' 
        Set up axes
        '''
        self.sub.set_ylim([9.0, 12.0])
        self.sub.set_xlim([10.0, 15.0])
        
        self.sub.set_ylabel(r'Stellar Mass $\mathtt{M_*}$', fontsize=25) 
        self.sub.set_xlabel(r'Halo Mass $\mathtt{M_{Halo}}$', fontsize=25) 

        self.sub.legend(loc='upper left', frameon=False, scatterpoints=1)

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
        ''' Child class of Plots class that plots the StarForming Main Sequence 
        '''
        super(PlotSFMS, self).__init__(**kwargs)

    def plot(self, mass=None, sfr=None, sfr_class=None, 
            gal_class='star-forming', bovyplot=True, sigSFR=True, **mkwargs): 
        ''' Plot SFR as a function of stellar mass for star-forming 
        galaxies 
        '''
        if mass is None or sfr is None:
            raise ValueError

        if sfr_class is None and gal_class != 'all': 
            raise ValueError
        
        if self.kwargs == {}: 
            self.kwargs = mkwargs.copy() 
        else: 
            self.kwargs = (self.kwargs).copy()
            self.kwargs.update(mkwargs)
        # label 
        if 'label' in self.kwargs: 
            label = self.kwargs['label']
        else: 
            label = None
        # color 
        if 'color' in self.kwargs: 
            if isinstance(self.kwargs['color'], str): 
                color = self.kwargs['color']
            elif isinstance(self.kwargs['color'], int): 
                color = self.pretty_colors[self.kwargs['color']]
        else: 
            color = 'black'
        
        if gal_class != 'all': 
            just_gal_class = np.where(sfr_class == gal_class)
            mass = mass[just_gal_class]
            sfr = sfr[just_gal_class] 
        
        if not bovyplot: 
            self.sub.scatter(mass, sfr, color=color, s=3, label=label)
            
            if sigSFR: 
                mbin = np.arange(7.0, 12.25, 0.25)
                mlow = mbin[:-1]
                mhigh = mbin[1:]
                muSFR, sigSFR = [], [] 
                for im in range(len(mlow)): 
                    mbin = np.where((mass > mlow[im]) & (mass <= mhigh[im]))
                    muSFR.append(np.mean(sfr[mbin]))
                    sigSFR.append(np.std(sfr[mbin]))
                self.sub.plot(0.5*(mlow+mhigh), muSFR, color='k', ls='--', lw=4)
                self.sub.errorbar(0.5*(mlow+mhigh), muSFR, yerr=sigSFR, color='k')

            self.sub.set_xlim([9.0, 12.0])
            self.sub.set_ylim([-5.0, 2.0])
            self.sub.set_xlabel(r'$\mathtt{log\;M_*}$', fontsize=25)
            self.sub.set_ylabel(r'$\mathtt{log\;SFR}$', fontsize=25)
        else: 
            bovy.scatterplot(mass, sfr, scatter=True, color=color, 
                    s=3, 
                    levels=[0.68, 0.95, 0.997], #1,2,3 sigmas
                    xrange=[9.0, 12.0], 
                    yrange=[-5.0, 2.0], 
                    xlabel='\mathtt{M_*}', 
                    ylabel='\mathtt{log\;SFR}'
                    )
            #plt.plot(np.arange(9.0, 12.0, 0.1), 
            #        sfr_cut(np.arange(9.0, 12.0, 0.1), cq_obj.zsnap), 
            #        ls='--', lw=3, c='r')
        return None   

    def model(self, z=None, bovyplot=True, **pltkwargs): 
        ''' Plot the model SFMS
        '''
        if z is None: 
            raise ValueError("Specify Redshift") 
        
        mstar = np.arange(8.5, 12.0, 0.1)

        if not bovyplot: 
            self.sub.plot(mstar, AverageLogSFR_sfms(mstar, z),
                    c = 'k', 
                    ls = '--', 
                    lw = 3)
        else:
            plt.plot(
                    mstar, 
                    AverageLogSFR_sfms(mstar, z),
                    c = 'k', 
                    ls = '--', 
                    lw = 3)

        return None
        
    def save_fig(self, file_name): 
        '''
        Save figure to file 
        '''
        plt.savefig(file_name, bbox_inches='tight') 
        return None



# Stellar Mass Function 
class PlotSMF(Plots): 
    def __init__(self, **kwargs): 
        ''' Child class of Plots class that describes SMF plots
        '''
        super(PlotSMF, self).__init__(**kwargs)

    def plot(self, mass=None, smf=None, **kwargs): 
        ''' Plot Stellar Mass Function  
        '''
        if 'label' in kwargs: 
            smf_label = kwargs['label']
            kwargs.pop('label', None)
        else: 
            smf_label = ''
    
        if 'line_color' in kwargs: 
            if isinstance(kwargs['line_color'], str):
                line_color = kwargs['line_color']
            elif isinstance(kwargs['line_color'], int):
                line_color = self.pretty_colors[kwargs['line_color']]
            kwargs.pop('line_color', None)
        else: 
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
        
        if smf is None: 
            mf = SMF()
            masses, phi = mf._smf(mass, **kwargs)
        else: 
            masses, phi = smf

        self.sub.plot(masses, phi, 
                c=line_color, ls=line_style, lw=line_width,
                label=smf_label)
        return mass, phi 

    def model(self, redshift, source='li-drory-march', **kwargs): 
        ''' Plot the model (Analytic) SMF at given redshift 
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

    def CentralSubhalo(self, nsnap, subhalo_prop=None):
        ''' Plot SMF of Central Subhalos
        '''
        if subhalo_prop is None: 
            raise ValueError

        censub = CentralSubhalos()
        censub.Read(nsnap, 
                scatter=subhalo_prop['scatter'], 
                source=subhalo_prop['source'], 
                nsnap_ancestor=subhalo_prop['nsnap_ancestor'])
        
        self.plot(mass=censub.mass, 
                line_color='k', line_style='-', label='Central Subhalos')

        return None

    def set_axes(self): 
        ''' Set up default figure axes 
        '''
        self.sub.set_ylim([10**-5, 10**-1])
        self.sub.set_xlim([7.5, 12.0])
        self.sub.set_yscale('log')

        # x and y labels
        self.sub.set_xlabel(r'Mass $\mathtt{M_*}$', fontsize=25) 
        self.sub.set_ylabel(r'Stellar Mass Function $\mathtt{\Phi}$', fontsize=25) 
        
        self.sub.legend(loc='upper right', frameon=False)
        return None



'''
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
'''

