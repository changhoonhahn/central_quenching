'''

Plot the quiescent fraction of SF assigned or Evolved
CenQue class object

'''
import numpy as np
import matplotlib.pyplot as plt

#----- Local -----
from plots import Plots
#from cenque import CenQue
# quiescent fraction
from quiescent_fraction import cq_fq
from quiescent_fraction import get_fq

from util import cenque_utility as util

from group_catalog.group_catalog import central_catalog


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
        """ 
        Parameterized queiscent fraction  
        """

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
        """ Set up axes
        """
        self.subs.set_xlim([9.0, 12.0])
        self.subs.set_ylim([0.0, 1.0])
        
        self.subs.set_xlabel(r'Mass $\mathtt{M_*}$') 
        self.subs.set_ylabel(r'Quiescent Fraction $\mathtt{f_Q}$', fontsize=20) 

        self.subs.legend(loc='upper left', frameon=False)

        return None
    
    def save_fig(self, file_name): 
        ''' 
        save figure to file 
        '''
        self.fig.savefig(file_name, bbox_inches='tight')
"""
def plot_fqobs_snapshot_evol( 
        nsnaps = [1,2,3,4,5,6,7,8,9,10,11,12], sf_prop={'name': 'average'}, fq_prop={'name': 'wetzelsmooth'}, cenque_type='sf_assigned', tau_prop={'name': 'instant'}, **kwargs): 
    ''' Plot the observed quiescent fraction of snapshots

    Parameters
    ----------
    nsnaps : list of snapshots
    fqtype : type of queiscent fraction  
    '''
    #prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[nsnaps]

    mass_bin = np.arange(9.0, 12.0, 0.25)       # mass bin 
    mass_low = mass_bin[:-1]
    mass_high = mass_bin[1:]
    mass_mid = 0.5 * (mass_low + mass_high)
       
    fig, subs = plt.subplots(1,2, figsize=[10, 5]) 
    subs = subs.ravel() 

    if 'evol_from' in cenque_type: 

        start_nsnap = cenque_type.split('evol_from')[-1]

        snap = CenQue(n_snap = start_nsnap) 
        snap.cenque_type = 'sf_assigned'
        snap.tau_prop = tau_prop
        snap.readin()

        fq_mass, fq = [], [] 

        for i_m in xrange(len(mass_mid)): 

            mass_bin_index = np.where(
                    (snap.mass > mass_low[i_m]) & 
                    (snap.mass <= mass_high[i_m])
                    )
            ngal = len(mass_bin_index[0])
            
            if ngal == 0: 
                continue 

            sfqs = sfq_classify( 
                    snap.mass[mass_bin_index], 
                    snap.sfr[mass_bin_index], 
                    snap.zsnap 
                    ) 

            ngal_q = np.float(np.sum(sfqs == 'quiescent')) 

            fq_mass.append(mass_mid[i_m]) 
            fq.append(ngal_q/ngal) 
        
        subs[0].plot(fq_mass, fq, color='black', ls='--', lw=4) 
    
    for ii_nsnap, i_nsnap in enumerate(nsnaps): 

        snap = CenQue(n_snap = i_nsnap ) 
        snap.cenque_type = cenque_type
        snap.fq_prop = fq_prop
        snap.sf_prop = sf_prop
        snap.tau_prop = tau_prop
        snap.readin() 

        fq_mass, fq = [], [] 
        for i_m in xrange(len(mass_mid)): 

            mass_bin_index = np.where(
                    (snap.mass > mass_low[i_m]) & 
                    (snap.mass <= mass_high[i_m])
                    )
            ngal = len(mass_bin_index[0])
            
            if ngal == 0: 
                continue 

            sfqs = sfq_classify( 
                    snap.mass[mass_bin_index], 
                    snap.sfr[mass_bin_index], 
                    snap.zsnap 
                    ) 
            ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
            fq_mass.append(mass_mid[i_m]) 
            fq.append(ngal_q/ngal) 
        
        subs[0].plot(fq_mass, fq, color=pretty_colors[i_nsnap-1], lw=4) 
        
        # plot fq(Mass) 
        fq_mass = get_fq(mass_mid, zbin[ii_nsnap], lit=fq_prop['name']) 
        subs[0].plot(
                mass_mid, 
                fq_mass, 
                color=pretty_colors[i_nsnap-1], 
                lw=2, 
                ls='--' 
                ) 
    
    subs[0].set_title('Snapshots')
    subs[0].set_xlim([9.0, 12.0])
    subs[0].set_ylim([0.0, 1.0])
    subs[0].set_xlabel('Mass') 
    for i_z, z in enumerate(zbin): 

        # plot fq(Mass) 
        fq_mass = get_fq(mass_mid, z, lit=fq_prop['name']) 
        subs[1].plot(
                mass_mid, 
                fq_mass, 
                color=pretty_colors[i_z], 
                lw=4, 
                label=str(z) 
                ) 
    
    subs[1].set_title(fq_prop['name']) 

    subs[1].set_xlim([9.0, 12.0])
    subs[1].set_ylim([0.0, 1.0])

    subs[1].set_xlabel('Mass') 

    subs[0].set_ylabel('Quiescent Fraction') 
    #plt.legend(loc='upper left')
    #plt.show()
    
    if tau_prop['name'] in ('instant', 'constant', 'satellite', 'long'): 
        tau_str = ''.join(['_', tau_prop['name'], 'tau'])
    elif tau_prop['name'] in ('line'):
        tau_str = ''.join([
            '_', tau_prop['name'], 'tau', 
            '_Mfid', str(tau_prop['fid_mass']), 
            '_slope', str(round(tau_prop['slope'], 4)), 
            '_yint', str(round(tau_prop['yint'],4))
            ])
    fig.savefig(
            ''.join([
                '/home/users/hahn/research/pro/tinker/central_quenching/figure/', 
                'fq_', 
                tau_str, 
                '.png']), 
            bbox_inches='tight'
            )
    fig.clear()
    plt.close()

def plot_fqobs_snapshot(nsnap, fq_prop={'name': 'wetzelsmooth'}, cenque_type='sf_assigned', tau_prop={'name': 'instant'}, **kwargs):
    ''' 
    Plot the quiescent fraction of specified CenQue class object
    along with a parameterized fQ for SINGLE given snapshot redshift

    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    nsnap : (int) snapshot #
    fq_type : type of queiscent fraction  

    '''
    #prettyplot()        # make pretty 
    pretty_colors = prettycolors() 
    
    # Corresponding redshift for the snapshot 
    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[nsnap]

    mass_bin = np.arange(9.0, 12.0, 0.2)   # mass bin 
    mass_low = mass_bin[:-1]
    mass_high = mass_bin[1:]
    mass_mid = 0.5 * (mass_low + mass_high)
       
    fig = plt.figure(figsize=[5, 5]) 
    subs = fig.add_subplot(111)

    snap = CenQue(n_snap = nsnap) 
    snap.cenque_type = cenque_type 
    snap.fq_prop = fq_prop
    snap.tau_prop = tau_prop
    snap.readin() 
        
    # classify CenQue class object with SF properties using 
    # sfq_classify function 
    sfqs = sfq_classify( 
            snap.mass, 
            snap.sfr, 
            snap.zsnap 
            ) 

    mass_fq, fq = [], [] 

    for i_m in xrange(len(mass_mid)): 

        mass_bin_index = np.where(
                (snap.mass > mass_low[i_m]) & 
                (snap.mass <= mass_high[i_m])
                )
        
        ngal = np.float(len(mass_bin_index[0]))

        if not ngal > 0.: 
            continue 
        
        ngal_q = np.float(np.sum(sfqs[mass_bin_index] == 'quiescent')) 

        mass_fq.append(mass_mid[i_m]) 
        fq.append(ngal_q/ngal) 
    
    subs.plot(
            mass_fq, fq, 
            color = pretty_colors[3], 
            lw = 4, 
            label='Snapshot '+str(snap.nsnap)
            ) 
    
    # parameterized fq 
    fq_mass = get_fq(
            mass_mid, zbin, 
            lit = fq_prop['name']
            )
    subs.plot(
            mass_mid, fq_mass, 
            color = 'black',  
            lw = 4, 
            ls = '--', 
            label = fq_prop['name']+'; z = '+str(zbin) 
            ) 

    subs.set_xlim([9.0, 12.0])
    subs.set_ylim([0.0, 1.0])
    subs.set_xlabel('Mass') 
    subs.set_ylabel('Quiescent Fraction') 
    subs.legend(loc='upper left') 

    plt.show()

def plot_fq_geha_groupcat(Mrcut=18): 
    ''' Plot comparison of Modified Tinker Catalog fq from Geha and SDSS group catalog 

    Notes
    -----
    * Mrcut should only be = 18 because the modified tinker catalog is constructed from Mrcut=-18 with z <= 0.06. 
    * SDSS Group Catalog fQ is calculated based on SFR-M* cut. 
    * Geha Modified Tinker catalog fQ is calculated based on EW Halpha and Dn4000. 
    * Should they match?  

    '''
    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    if Mrcut != 18: 
        print 'Mrcut should only be = 18 because the modified tinker catalog is constructed from Mrcut=-18 with z <= 0.06.'
    
    # load literature data 
    # modified tinker
    mod_tink_file = ''.join(['dat/central_quenching/literature/modified_tinker_fq.dat']) 
    mod_tink_mass, mod_tink_fq = np.loadtxt(mod_tink_file, unpack=True, usecols=[0,1])   
    
    # read group catalog 
    central = cq_group.central_catalog(Mrcut=Mrcut) 
    central_z = central.z    # redshifts 
    median_z = np.median(central_z)
   
    fq_mass, fq = [], [] 

    mass_bins = util.simple_mass_bin() 
    for i_m in range(mass_bins.nbins): 

        mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
        mass_bin_mid = round(mass_bins.mass_mid[i_m], 2)
        mass_bin_high = round(mass_bins.mass_high[i_m], 2) 

        # boolean list for mass range
        mass_bin_bool = (central.mass > mass_bin_low) & (central.mass <= mass_bin_high)
        
        if np.sum(mass_bin_bool) < 10: 
            continue 

        sfqs = sfq_classify( central.mass[mass_bin_bool], central.sfr[mass_bin_bool], median_z ) 
        ngal = np.float(len(sfqs))
        ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
        fq_mass.append(mass_bin_mid) 
        fq.append(ngal_q/ngal) 

    fig = plt.figure(figsize=[8, 8]) 
    sub = fig.add_subplot(111)

    sub.plot(mod_tink_mass, mod_tink_fq, 
        color='black', lw=6, ls='--', label='Modified Tinker Group' ) 
    
    sub.plot(fq_mass, fq, 
            color='black', lw=6, label=r'SDSS Group Catalog; z = '+str(median_z))

    sub.set_xlim([9.0, 12.0])
    sub.set_ylim([0.0, 1.0])

    sub.set_xlabel('Mass') 
    sub.set_ylabel('Quiescent Fraction') 
    sub.legend(loc='upper left') 

    fig_name = ''.join(['figure/fq_evol_groupcat', str(Mrcut), '_geha_comp.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear() 

def plot_fq_evol(): 
    ''' Plot fq evolution comparison
    '''
    prettyplot()    
    pretty_colors = prettycolors() 

    # snapshot redshifts
    zbin, tcosmic = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2,3])
    ofinterest = np.where(zbin < 1.6)
    zbin = zbin[ofinterest]
    tcosmic = tcosmic[ofinterest]

    mass_bin = np.arange(9.0, 12.0, 0.2)   # mass bin 
    mass_low = mass_bin[:-1]
    mass_high = mass_bin[1:]
    mass_mid = 0.5 * (mass_low + mass_high)
    
    # load modified tinker group catalog data from Marla
    mod_tink_file = ''.join(['dat/central_quenching/literature/modified_tinker_fq.dat']) 
    mod_tink_mass, mod_tink_fq = np.loadtxt(
            mod_tink_file, 
            unpack=True, 
            usecols=[0,1]
            )   
    gc = central_catalog(Mrcut=18, clobber=False)
    groupcat = CenQue()
    groupcat.sfr = gc.sfr
    groupcat.mass = gc.mass
    groupcat.zsnap = np.median(gc.z)
    gc_mass, gc_fq = cq_fq(groupcat)

    #fq_types = ['cosmosinterp', 'wetzel', 'wetzelsmooth'] 
    fq_types = ['wetzelsmooth']
    
    fig = plt.figure(1, figsize=[len(fq_types)*5, 5]) 
    subs = [fig.add_subplot(1, len(fq_types), i) for i in xrange(len(fq_types))]

    for i_fq, fq_type in enumerate(fq_types): 

        fq_s = [] 
        for i_z, z in enumerate(zbin): 

            # plot fq(Mass) 
            fq_mass = get_fq(mass_mid, z, lit=fq_type) 
            subs[i_fq].plot(
                    mass_mid, 
                    fq_mass, 
                    color=pretty_colors[i_z], 
                    lw=4, 
                    label='z = '+str(z) 
                    ) 
            fq_s.append(fq_mass)

        subs[i_fq].plot(
                mod_tink_mass, 
                mod_tink_fq, 
                color='black', 
                lw=6, 
                label='Modified Tinker Group' 
                ) 

        subs[i_fq].plot(
                gc_mass, 
                gc_fq, 
                color='red', 
                lw=4, 
                ls='--', 
                label=r'Group Catalog $\mathtt{M_r = 18}$' 
                ) 

        subs[i_fq].set_title(fq_type) 

        subs[i_fq].set_xlim([9.0, 12.0])
        subs[i_fq].set_ylim([0.0, 1.0])

        subs[i_fq].set_xlabel('Mass') 
    subs[0].set_ylabel('Quiescent Fraction') 
    #subs[0].legend(loc='upper left') 

    fig2 = plt.figure(2, figsize=[14,10]) 
    subs2 = fig2.add_subplot(1,1,1)
    for i_mass, mass in enumerate(mass_mid): 
        subs2.plot(
                tcosmic, 
                [fq[i_mass] for fq in fq_s], 
                lw=4, 
                c=pretty_colors[i_mass], 
                label = r'$\mathtt{log\;M_* = '+str(mass)+'}$'
                )
        subs2.scatter(
                tcosmic, 
                [fq[i_mass] for fq in fq_s], 
                s=12, 
                c=pretty_colors[i_mass]
                )
    subs2.set_xlim([4.0, 14.0])
    subs2.set_xlabel(r'$\mathtt{t_{cosmic}}$', fontsize=20)
    subs2.set_ylabel(r'Quiescent Fraction')
    #subs2.set_yscale('log')
    #subs2.set_ylim([0.01, 1.0])
    subs2.set_ylim([0.0, 1.0])
    subs2.legend(loc='upper left', frameon=False)

    fig_name = ''.join(['figure/fq_evol_comp_', '_'.join(fq_types), '.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear() 

    fig_name2 = ''.join(['figure/fq_tcosmic_evol.png'])
    fig2.savefig(fig_name2, bbox_inches='tight')

if __name__=="__main__":
    plot_fq_evol()
    #plot_fqobs_snapshot_evol( nsnaps = [1,2,3,4,5,6,7,8,9,10,11,12], 
    #        cenque_type = 'evol_from13', 
    #        fq_prop = {'name': 'wetzelsmooth'}, 
    #        sf_prop = {'name': 'average'}, 
    #        tau_prop={'name': 'line', 'fid_mass': 10.75, 'slope': -0.57, 'yint': 0.5}
    #        )
    #
    #for i_nsnap in [12,11,10,9,8,7,6,5,4,3,2,1]:
    #    plot_fqobs_snapshot(i_nsnap, fq_type='wetzelsmooth', cenque_type='evol_from13', tau_prop={'name': 'satellite'})
"""
