import numpy as np
import matplotlib.pyplot as plt

#----- Local -----
from cenque import CenQue
from quiescent_fraction import get_fq
from util import cenque_utility as util
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
import bovy_plot as bovy

# quiescent fraction -----------------------------------
def plot_snapshot_fqobs_evol(
        nsnaps = [1,2,3,4,5,6,7,8,9,10,11,12],
        fq_type = 'wetzelsmooth', 
        **kwargs
        ): 
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
    print mass_mid
       
    fig, subs = plt.subplots(1,2, figsize=[10, 5]) 
    subs = subs.ravel() 

    snap = CenQue(n_snap = 13) 
    snap.cenque_type = 'sf_assigned'
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

        sfqs = util.sfq_classify( 
                snap.mass[mass_bin_index], 
                snap.sfr[mass_bin_index], 
                snap.zsnap 
                ) 

        ngal_q = np.float(np.sum(sfqs == 'quiescent')) 

        fq_mass.append(mass_mid[i_m]) 
        fq.append(ngal_q/ngal) 
    
    subs[0].plot(fq_mass, fq, color='black', ls='--', lw=4, label='z= '+str(snap.zsnap)) 
    
    for i_nsnap in nsnaps: 
        snap = CenQue(n_snap = i_nsnap ) 
        snap.cenque_type = 'evol_from13'
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

            sfqs = util.sfq_classify( 
                    snap.mass[mass_bin_index], 
                    snap.sfr[mass_bin_index], 
                    snap.zsnap 
                    ) 
            ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
            fq_mass.append(mass_mid[i_m]) 
            fq.append(ngal_q/ngal) 
        
        subs[0].plot(fq_mass, fq, color=pretty_colors[i_nsnap-1], lw=4, label='z= '+str(snap.zsnap)) 
    
    subs[0].set_title('Snapshots')
    subs[0].set_xlim([9.0, 12.0])
    subs[0].set_ylim([0.0, 1.0])
    subs[0].set_xlabel('Mass') 
    for i_z, z in enumerate(zbin): 

        # plot fq(Mass) 
        fq_mass = get_fq(mass_mid, z, lit=fq_type) 
        subs[1].plot(
                mass_mid, 
                fq_mass, 
                color=pretty_colors[i_z], 
                lw=4, 
                label='z = '+str(z) 
                ) 
    
    subs[1].set_title(fq_type) 

    subs[1].set_xlim([9.0, 12.0])
    subs[1].set_ylim([0.0, 1.0])

    subs[1].set_xlabel('Mass') 

    subs[0].set_ylabel('Quiescent Fraction') 
    plt.legend()
    plt.show()
    """    
    # file name ----------------------------------------------------------------
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

    fig_file = ''.join(['../figure/', 
        'fq_obs_snapshots_', tau_str, '_', fqing_str, file_type_str, '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear()
    plt.close(fig)
    """

def plot_snapshot_fqobs(nsnap, fq_type='wetzel', **kwargs): 
    ''' Plot the observed quiescent fraction of a snapshot with parameterized fQ

    Parameters
    ----------
    nsnap : (int) snapshot #
    fqtype : type of queiscent fraction  

    '''
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 

    # snapshot redshifts
    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[nsnap]

    mass_bins = util.simple_mass_bin()                    # mass bin 
       
    fig = plt.figure(figsize=[5, 5]) 
    subs = fig.add_subplot(111)

    snap = cq.CenQue() 
    snap.readin(nsnap=nsnap, file_type='evol from 13', **kwargs) 

    fq_mass, fq = [], [] 

    for i_m in range(mass_bins.nbins): 
        mass_bin_low = round(mass_bins.mass_low[i_m], 2)        # for floating point errors
        mass_bin_mid = round(mass_bins.mass_mid[i_m], 2)
        mass_bin_high = round(mass_bins.mass_high[i_m], 2) 

        # boolean list for mass range
        mass_bin_bool = (snap.mass > mass_bin_low) & (snap.mass <= mass_bin_high)
        
        if np.sum(mass_bin_bool) == 0: 
            continue 

        sfqs = util.sfq_classify( snap.mass[mass_bin_bool], snap.sfr[mass_bin_bool], snap.zsnap ) 
        ngal = np.float(len(sfqs))
        ngal_q = np.float(np.sum(sfqs == 'quiescent')) 
        fq_mass.append(mass_bin_mid) 
        fq.append(ngal_q/ngal) 
    
    subs.plot(fq_mass, fq, color=pretty_colors[3], lw=4, label='Snapshot '+str(snap.nsnap)) 
    
    # parameterized fq 
    fq_mass = [util.get_fq(mass_bins.mass_mid[i], zbin, lit=fq_type) for i in range(len(mass_bins.mass_mid))]
    subs.plot(mass_bins.mass_mid, fq_mass, 
            color='black', lw=4, ls='--', label='Wetzel; z = '+str(zbin) ) 

    subs.set_xlim([9.0, 12.0])
    subs.set_ylim([0.0, 1.0])
    subs.set_xlabel('Mass') 
    subs.set_ylabel('Quiescent Fraction') 
    subs.legend(loc='upper left') 

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
        'fq_obs_snapshot', str(nsnap), '_', tau_str, '_', fqing_str, file_type_str, '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear()
    plt.close(fig)

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

        sfqs = util.sfq_classify( central.mass[mass_bin_bool], central.sfr[mass_bin_bool], median_z ) 
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

    # snapshot redshifts
    zbin = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2])
    zbin = zbin[zbin < 1.0]
    #zbin = [0.1*np.float(i) for i in range(1,10)]   # zbin 

    mass_bin = util.simple_mass_bin()                    # mass bin 
       
    prettyplot()        # make pretty 
    pretty_colors = prettycolors() 
    
    # load literature data 
    # modified tinker
    mod_tink_file = ''.join(['dat/central_quenching/literature/modified_tinker_fq.dat']) 
    mod_tink_mass, mod_tink_fq = np.loadtxt(mod_tink_file, unpack=True, usecols=[0,1])   

    fq_types = ['cosmosinterp', 'wetzel', 'wetzelsmooth'] 
    
    fig, subs = plt.subplots(1, len(fq_types), figsize=[len(fq_types)*5, 5]) 
    subs = subs.ravel() 

    for i_fq, fq_type in enumerate(fq_types): 
        for i_z, z in enumerate(zbin): 

            # plot fq(Mass) 
            fq_mass = [util.get_fq(mass_bin.mass_mid[i], z, lit=fq_type) 
                    for i in range(len(mass_bin.mass_mid))]
            subs[i_fq].plot(mass_bin.mass_mid, fq_mass, 
                    color=pretty_colors[i_z], lw=4, label='z = '+str(z) ) 
        
        subs[i_fq].plot(mod_tink_mass, mod_tink_fq, 
            color='black', lw=6, label='Modified Tinker Group' ) 

        subs[i_fq].set_title(fq_type) 

        subs[i_fq].set_xlim([9.0, 12.0])
        subs[i_fq].set_ylim([0.0, 1.0])

        subs[i_fq].set_xlabel('Mass') 

    subs[0].set_ylabel('Quiescent Fraction') 
    subs[0].legend(loc='upper left') 
    fig_name = ''.join(['figure/fq_evol_comp_', '_'.join(fq_types), '.png'])
    fig.savefig(fig_name, bbox_inches='tight')
    fig.clear() 

if __name__=="__main__":
    plot_snapshot_fqobs_evol(nsnaps = [2,3,4,5,6,7,8,9,10,11,12])

