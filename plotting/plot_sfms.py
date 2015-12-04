'''

Plotting for SF-MS module 


'''
import numpy as np
import scipy as sp
import os.path
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# --- Local ---
from plots import Plots
import bovy_plot as bovy
from group_catalog.group_catalog import sf_centrals

class PlotSFMS(Plots): 

    def __init__(self, **kwargs): 
        ''' 
        Child class of Plots class that plots the StarForming Main Sequence for 
        different class objects (CenQue/GroupCat)
        '''

        super(PlotSFMS, self).__init__(**kwargs)

    def cenque(self, cq_obj, **mkwargs): 
        '''
        Plot SF-MS (M* vs SFR) for CenQue object
        '''

        if self.kwargs == {}: 
            kwargs = mkwargs.copy() 
        else: 
            kwargs = (self.kwargs).copy()
            kwargs.update(mkwargs)
    
        if 'label' in kwargs: 
            label = kwargs['label']
        else: 
            label = 'z = '+str(cq_obj.zsnap) 
        
        if 'color' in kwargs: 
            color = kwargs['color']
        else: 
            try: 
                color = self.pretty_colors[cq_obj.nsnap]
            except TypeError: 
                color = 'black'

        if 'justsf' in kwargs: 

            if not kwargs['justsf']: 
                raise ValueError(
                        'Why would you specify it; to make it false, just leave it blank'
                        )
            
            justsf = np.where(cq_obj.gal_type == 'star-forming')
            
            self.sub.scatter(
                    cq_obj.mass[justsf],
                    cq_obj.sfr[justsf],
                    color=color, 
                    s=3
                    )
            self.sub.set_xlim([9.0, 12.0])
            self.sub.set_ylim([-5.0, 2.0])
            self.sub.set_xlabel(r'$\mathtt{M_*}$')
            self.sub.set_ylabel(r'$\mathtt{log\;SFR}$')

        else: 

            bovy.scatterplot(
                    cq_obj.mass, 
                    cq_obj.sfr,
                    scatter=True, 
                    color=color, 
                    s=3, 
                    xrange=[9.0, 12.0], 
                    yrange=[-5.0, 2.0], 
                    xlabel='\mathtt{M_*}', 
                    ylabel='\mathtt{log\;SFR}'
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

    def save_fig(self, file_name): 
        '''
        Save figure to file 
        '''
        
        plt.savefig(file_name, bbox_inches='tight') 

        return None

def plot_sfms_data(Mrcut=18): 
    ''' Plot StarForming Main Sequence from iSEDfit data and flexible SFMS fit function 

    Parameters
    ----------
    Mrcut : absolute magnitude to specify the group catalog 

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 
    zbins = [ (0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8)]  # zbin
    
    # figure with zbin panels 
    fig, subs = plt.subplots(1, len(zbins), figsize=[20, 5]) 
    subs = subs.ravel() 

    fits = [] 
    for i_z, zbin in enumerate(zbins): 
        
        if i_z == 0:    # 2D histogram of SF centrals  

            centrals = sfms.sf_centrals(Mrcut=Mrcut)    # SF centrals
            subs[i_z].hist2d(centrals.mass, centrals.sfr, bins=[30, 50])
        
        # avg SFR of envcount SF galaxies 
        mass, avg_sfrs, var_sfrs = [], [], [] 
        for i_mass in range(len(mass_bin.mass_low)): 

            avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z_envcount(mass_bin.mass_mid[i_mass], 
                    0.5*(zbin[0] + zbin[1]))

            if ngal < 10: 
                continue    # skip mass bin if there aren't many galaxies
            if mass_bin.mass_mid[i_mass] < 9.5: 
                continue    # skip low mass

            mass.append(mass_bin.mass_mid[i_mass]) 
            avg_sfrs.append(avg_sfr)
            var_sfrs.append(var_sfr)

        subs[i_z].errorbar(mass, avg_sfrs, yerr=var_sfrs, 
                lw=4, c=pretty_colors[1], label='Avg SFR (EnvCount)')
    
        # bestfit line for envcount SF-MS 
        if i_z == 0: 
            param = sfms.get_bestfit_envcount_sfms()
            subs[i_z].plot(mass, param[0]*(np.array(mass)-10.5) + param[1], 
                    c=pretty_colors[2])

        # avg SFR for SF Group Catalog 
        if i_z == 0:  
            mass, avg_sfrs, var_sfrs = [], [], [] 
            for i_mass in range(len(mass_bin.mass_low)): 

                avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z_groupcat(mass_bin.mass_mid[i_mass], Mrcut=Mrcut)

                if ngal < 10: 
                    continue    # skip mass bin if there aren't many galaxies
                if mass_bin.mass_mid[i_mass] < 9.5: 
                    continue    # skip low mass

                mass.append(mass_bin.mass_mid[i_mass]) 
                avg_sfrs.append(avg_sfr)
                var_sfrs.append(var_sfr)

            subs[i_z].errorbar(mass, avg_sfrs, yerr=var_sfrs, 
                    lw=4, c=pretty_colors[3], label='Avg SFR (GroupCat)')

        # bestfit line for SF Group Catalog 
        bestfit_sfr = [] 
        for m in mass: 
            sfcen_bestfit = util.get_sfr_mstar_z_bestfit(m, 0.5*(zbin[0]+zbin[1]), clobber=True)
            bestfit_sfr.append(sfcen_bestfit[0]) 

        subs[i_z].plot(mass, bestfit_sfr, c=pretty_colors[3])
        
        subs[i_z].text(10.0, 1.25, '$\mathtt{z \sim '+str(0.5*(zbin[0]+zbin[1]))+'}$') 

        subs[i_z].set_xlim([9.5, 12.0]) 
        subs[i_z].set_ylim([-1.5, 1.5]) 
        if i_z in (1, 2):
            subs[i_z].set_xlabel('log(M)') 
        if i_z == 0: 
            subs[i_z].set_ylabel('log(SFR)') 
        #if i_z == 3: 
        #    subs[i_z].legend()
    
    fig_file = ''.join(['figure/', 'sf_ms_data.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear() 

def plot_sfms_groupcat(Mrcut=18, writeout=True):
    ''' StarForming Main Sequence (mass vs SFR) of Star-Forming SDSS group catalog  

    Parameters
    ----------
    Mrcut : Absolute magnitude cut off of SDSS group catalog 

    '''
    # star-forming central galaxies from the SDSS group catalog
    sf_cen = sf_centrals(Mrcut=Mrcut) 

    #prettyplot()        
    pretty_colors = prettycolors() 

    #mass_bin = util.simple_mass_bin()       # mass bin 
    bovy.scatterplot(
            sf_cen.mass,
            sf_cen.sfr,
            scatter = True,
            levels = [0.68, 0.95, 0.997],
            color = pretty_colors[1],
            s = 3,
            xrange = [9.5, 12.0],
            yrange = [-1.5, 1.5],
            xlabel = 'log \; M_{*}',
            ylabel = 'log \; SFR'
            )

    #fig = plt.figure(figsize = (8,8))
    #sub = fig.add_subplot(111)

    #sub.hist2d(
    #        sf_cen.mass, 
    #        sf_cen.sfr, 
    #        bins = [30, 50]
    #        )

    #sub.set_xlim([9.5, 12.0])
    #sub.set_ylim([-1.5, 1.5])
    #sub.set_xlabel(r'$\mathtt{log \; M_{*}}$', fontsize=20)
    #sub.set_ylabel(r'$\mathtt{log \; SFR}$', fontsize=20)
    
    if writeout: 
        fig_name = ''.join([
            'figure/', 
            'sfms_SDSS_groupcat_', str(Mrcut), '.png'
            ])
        plt.savefig(fig_name, bbox_inches='tight')

    return None

def plot_ssfms_groupcat(Mrcut=18):
    ''' Plot Specific SF Main Sequence (mass vs SFR) from Star-Forming SDSS group catalog  

    Parameters
    ----------
    None

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 

    fig = plt.figure(1)#, figsize=[7,5])
    subs = fig.add_subplot(111) 

    mass, avg_sfrs, var_sfrs, avg_ssfrs = [], [], [], [] 
    for i_mass in range(len(mass_bin.mass_low)): 
        avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z_groupcat(mass_bin.mass_mid[i_mass], 
                Mrcut=Mrcut) 

        if ngal < 10: 
            continue    # skip mass bin if there aren't many galaxies

        if mass_bin.mass_mid[i_mass] < 9.5: 
            continue    # skip low mass

        mass.append(mass_bin.mass_mid[i_mass]) 
        avg_sfrs.append(avg_sfr)
        var_sfrs.append(var_sfr)
        avg_ssfrs.append(avg_sfr-mass_bin.mass_mid[i_mass])
    
    # 2D histogram of centrals
    sf_cen = sfms.sf_centrals(Mrcut=Mrcut) 
    subs.hist2d(sf_cen.mass, sf_cen.ssfr, bins=[30, 50])

    subs.errorbar(mass, avg_ssfrs, yerr=var_sfrs, 
            lw=4, c=pretty_colors[1], label='Average SFR')

    bestfit_params = sfms.get_bestfit_groupcat_sfms(Mrcut=Mrcut, clobber=True) 
    subs.plot(mass, util.line(np.array(mass)-10.5, bestfit_params) - np.array(mass), 
            lw=4, ls='--', c=pretty_colors[2], label='MPfit line') 
    subs.text(11.25, -10.0, 'Slope = '+str("%.2f" % bestfit_params[0]))
    subs.text(11.25, -10.5, 'Y-int = '+str("%.2f" % bestfit_params[1]))

    subs.set_xlim([9.5, 12.0]) 
    subs.set_ylim([-9.5, -11.75]) 
    subs.set_xlabel('log(M)') 
    subs.set_ylabel('log(sSFR)') 
    subs.legend(loc='upper right') 
    
    fig_name = ''.join(['figure/', 'ssf_ms_groupcat_', str(Mrcut), '.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear()

def plot_q_groupcat(Mrcut=18): 
    ''' Plot mass vs sSFR for the Quiescet SDSS group catalog  

    Parameters
    ----------
    Mrcut : Absolute magnitude cut-off that specifies the group catalog

    Notes
    -----
    * Mainly used for determing the final sSFR of the quiescent population 
        * Because the SFR for the quiescent populations are just upper limits, this is just to make sure that the offset with the quiescent mode in the SSFR mode does not hinder the fitting procedure. 

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    mass_bin = util.simple_mass_bin()       # mass bin 

    fig = plt.figure(1)
    subs = fig.add_subplot(111) 
        
    q_cen = sfms.q_centrals(Mrcut=Mrcut) # Quiescent centrals  
    subs.hist2d(q_cen.mass, q_cen.ssfr, bins=[30, 1000])  # 2D histogram 
    
    med_ssfrs = [] 
    masses = np.arange(9.75, 11.5, 0.25) 
    for mass in masses: 
        med_ssfr, var, ngal = sfms.get_ssfr_mstar_qgroupcat(mass) 
        med_ssfrs.append(med_ssfr) 

    subs.plot(masses, med_ssfrs, 
            lw=4, c=pretty_colors[1]) 

    fit_line_param = sfms.get_bestfit_qgroupcat_ssfr(Mrcut=Mrcut, clobber=True) 
    subs.plot(masses, fit_line_param[0].item() * (masses - 10.5) + fit_line_param[1].item(), 
            lw=4, ls='--', c=pretty_colors[3]) 

    subs.set_xlim([9.5, 12.0]) 
    subs.set_ylim([-15.0, -10.0]) 
    subs.set_xlabel('log(M)') 
    subs.set_ylabel('log(sSFR)') 
    
    fig_name = ''.join(['figure/tinker/',
        'ssfr_mass_q_groupcat_mrcut', str(Mrcut), '.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear()

if __name__=="__main__": 
    plot_sfms_groupcat(Mrcut=18, writeout=True)
