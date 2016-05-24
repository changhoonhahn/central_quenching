'''

Module specifically to test the central_subhalo module.
The CentralSubhalos and Subhalos class objects in particular.  

'''
import h5py
import numpy as np
import os.path
from gal_prop import SMF
from sham_hack import LFClass

# local ---- 
from central_subhalo import Subhalos
from central_subhalo import CentralSubhalos
from util.util import get_z_nsnap
from util.util import intersection_index

# Plotting
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors 

def SubhaloSMF(type, scatter=0.0, source='li-drory-march', 
        nsnap_ancestor=20 ): 
    ''' Test the Subhalos/CentralSubhalos imported from TreePM by 
    comparing their measured SMFs to the analytic SMFs 

    Parameters
    ----------
    type : string 
        String that specifies whether we're interested in CentralSubhalos 
        or all Subhalos.
    scatter : float
        Float that specifies the scatter in the SMHM relation, which 
        affects the SHAM masses
    source : string
        String that specifies which SMF is used for SHAM
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    for i_snap in range(1, nsnap_ancestor+1):
        smf = SMF()
        if type == 'all': 
            subh = Subhalos()
        elif type == 'central': 
            subh = CentralSubhalos()
        else: 
            raise ValueError

        subh.Read(i_snap, scatter=scatter, source=source, nsnap_ancestor=nsnap_ancestor)
        
        subh_mass, subh_phi = smf.Obj(subh, dlogm=0.1, LF=False)
        sub.plot(subh_mass, subh_phi, lw=4, c=pretty_colors[i_snap % 19], alpha=0.5) 

        analytic_mass, analytic_phi = smf.analytic(get_z_nsnap(i_snap), source=source) 
        sub.plot(analytic_mass, analytic_phi, lw=4, ls='--', c=pretty_colors[i_snap % 19], 
                label=r"$ z = "+str(round(get_z_nsnap(i_snap),2))+"$") 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    #sub.legend(loc='upper right')
    if type == 'all': 
        subhalo_str = 'subhalo' 
    elif type == 'central': 
        subhalo_str = 'central_subhalo' 
    else: 
        raise ValueError
    plt.show()
    fig_file = ''.join(['figure/test/'
        'SubhaloSMF', 
        '.', subhalo_str, 
        '.scatter', str(round(scatter,2)), 
        '.ancestor', str(nsnap_ancestor),
        '.', source, 
        '.png']) 
    fig.savefig(fig_file, bbox_inches='tight')
    plt.close()

def SubhaloLF(scatter=0.0, source='cool_ages', 
        nsnap_ancestor=20 ): 
    ''' Test the Subhalos/CentralSubhalos imported from TreePM by 
    comparing their measured SMFs to the analytic SMFs 

    Parameters
    ----------
    scatter : float
        Float that specifies the scatter in the SMHM relation, which 
        affects the SHAM masses
    source : string
        String that specifies which SMF is used for SHAM
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    for i_snap in [11, 7, 4, 1]:
        smf = SMF()
        subh = Subhalos()

        subh.Read(i_snap, scatter=scatter, source=source, nsnap_ancestor=nsnap_ancestor)
        subh.mag_r *= -1.
        mass, phi = smf.Obj(subh, dlogm=0.1, m_arr=np.arange(-24.0, -16., 0.1), LF=True)

        sub.plot(mass, phi, lw=2, ls='-', c=pretty_colors[i_snap % 19])

        analytic_mass, analytic_phi = smf.analytic(get_z_nsnap(i_snap), source=source) 
        sub.plot(analytic_mass, analytic_phi, lw=4, ls='--', c=pretty_colors[i_snap % 19], 
                label=r"$ z = "+str(round(get_z_nsnap(i_snap),2))+"$") 

    sub.set_yscale('log')
    sub.set_ylim([10**-7, 10**-1])
    sub.set_xlim([-24.5, -17.8])
    sub.legend(loc='upper right')
    fig_file = ''.join(['figure/test/'
        'SubhaloLF', 
        '.scatter', str(round(scatter,2)), 
        '.ancestor', str(nsnap_ancestor),
        '.', source, 
        '.png']) 
    fig.savefig(fig_file, bbox_inches='tight')
    plt.close()

def Subhalo_MhaloMag(scatter=0.0, source='cool_ages', 
        nsnap_ancestor=20 ): 
    ''' Test the Subhalos/CentralSubhalos imported from TreePM by 
    comparing their measured SMFs to the analytic SMFs 

    Parameters
    ----------
    scatter : float
        Float that specifies the scatter in the SMHM relation, which 
        affects the SHAM masses
    source : string
        String that specifies which SMF is used for SHAM
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    for i_snap in [7, 4, 1]:

        mu_mag = []
        sig_mag = [] 

        subh = Subhalos()

        subh.Read(i_snap, scatter=scatter, source=source, nsnap_ancestor=nsnap_ancestor)
        subh.mag_r *= -1.
        
        m_halo = getattr(subh, 'halo.m.max') 
        for m_low in np.arange(10., 15.5, 0.1): 
            m_bin = np.where((m_halo > m_low) & (m_halo <= m_low+0.1))

            mu_mag.append(np.mean(subh.mag_r[m_bin]))
            sig_mag.append(np.std(subh.mag_r[m_bin]))
    
        #sub.errorbar(np.arange(10.05, 15.55, 0.1), mu_mag, yerr=sig_mag, c=pretty_colors[i_snap])
        sub.fill_between(np.arange(10.05, 15.55, 0.1), 
                np.array(mu_mag) - np.array(sig_mag), 
                np.array(mu_mag) + np.array(sig_mag), 
                color=pretty_colors[i_snap], alpha=0.75, 
                label=r"$ z = "+str(round(get_z_nsnap(i_snap),2))+"$") 
    sub.set_xlim([10, 15.5])
    sub.set_ylim([-17.8, -23.])
    sub.set_ylabel('$\mathtt{M}_\mathtt{halo}^\mathtt{peak}$', fontsize=25)
    sub.set_ylabel('Magnitude', fontsize=25)
    sub.legend(loc='lower right')
    fig_file = ''.join(['figure/test/'
        'Subhalo_Mhalo_Mag', 
        '.scatter', str(round(scatter,2)), 
        '.ancestor', str(nsnap_ancestor),
        '.', source, 
        '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()

def DescendantSubhaloSMF(type, nsnap_ancestor=20, scatter=0.0, source='li-drory-march', 
        nomass=False): 
    ''' SMF for 'descendant' all/central subhalos that have ancestors. If nomass is 
    specified, then we highlight the portion of the SMF from galaxies whose 
    ancestors have M* = 0 at snapshot nsnap_ancestor. 
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # Central subhalo at nsnap_ancestor
    if type == 'central': 
        subh = CentralSubhalos()
    elif type == 'all':
        subh = Subhalos()

    # read in ancestor snapshot  
    subh.Read(nsnap_ancestor, scatter=scatter, source=source, nsnap_ancestor=nsnap_ancestor)
    ancestor_index = subh.index
    ancestor_mass = subh.mass
    
    # galaxy has M* = 0 at nsnap_ancestor
    nomass = np.where(subh.mass == 0.0)
    nomass_anc_index = subh.index[nomass]
    nomass_anc_mass = subh.mass[nomass]
    # galaxy has M* > 0 at nsnap_ancestor
    massive = np.where(subh.mass > 0.0)
    massive_anc_index = subh.index[massive]
    massive_anc_mass = subh.mass[massive]

    for i_snap in [1, 5, 10, 15]:
        if i_snap == 1: 
            massive_label = 'Galaxies whose Ancestors have M* > 0'
            nomass_label = 'Galaxies whose Ancestors have M* = 0'
            label = 'Total'
        elif i_snap >= nsnap_ancestor: 
            continue 
        else: 
            massive_label = None 
            nomass_label = None
            label = None

        smf = SMF()
        # total central subhalo SMF
        if type == 'central': 
            subh = CentralSubhalos()
        elif type == 'all': 
            subh = Subhalos()

        subh.Read(i_snap, scatter=scatter, source=source, nsnap_ancestor=nsnap_ancestor)
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        sub.plot(subh_mass, subh_phi, lw=3, c=pretty_colors[i_snap], label=label) 
            
        anc_ind = getattr(subh, 'ancestor'+str(nsnap_ancestor))
        if not nomass: 
            # SMF of central subhalos who's ancestors are massive centrals at snapshot 20
            if i_snap == 1: 
                has_ancestor, has_descendant = intersection_index(anc_ind, massive_anc_index)
            else: 
                has_ancestor, has_d = intersection_index(anc_ind, massive_anc_index)
            mass, phi = smf.smf(subh.mass[has_ancestor])
            sub.plot(mass, phi, lw=2, ls='--', c=pretty_colors[i_snap], label=massive_label) 
        
        if nomass:
            # SMF of central subhalos who's ancestors have M*=0 centrals at snapshot 20
            if i_snap == 1: 
                has_ancestor, has_descendant = intersection_index(anc_ind, nomass_anc_index)
            else: 
                has_ancestor, has_d = intersection_index(anc_ind, nomass_anc_index)
            mass, phi = smf.smf(subh.mass[has_ancestor])
            sub.plot(mass, phi, lw=2, ls='--', c=pretty_colors[i_snap], label=nomass_label) 

        del subh

    anc_mass, anc_phi = smf.smf(ancestor_mass)
    sub.plot(anc_mass, anc_phi, lw=2, c='k', label='Initial Redshift') 
    #print 'With descendants at snapshot 1', len(ancestor_mass[has_descendant])
    #anc_massive = np.where(ancestor_mass[has_descendant] > 0.)
    #print 'With mass greater than 0', len(ancestor_mass[has_descendant[anc_massive]])
    #anc_mass, anc_phi = smf.smf(ancestor_mass[has_descendant])
    #sub.plot(anc_mass, anc_phi, ls='--', lw=4, c='k') 
    #anc_mass, anc_phi = smf.smf(ancestor_mass[has_descendant[anc_massive]])
    #sub.plot(anc_mass, anc_phi, ls='--', lw=2, c='green') 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    sub.legend(loc='upper right')
    if nomass: 
        nomass_str = '.ancestor_Mstar0'
    else: 
        nomass_str = ''
    
    fig_file = ''.join([
        'figure/test/'
        'DescendantSubhaloSMF', 
        '.', type, '_subhalo',  
        '.scatter', str(round(scatter, 2)), 
        '.ancestor', str(nsnap_ancestor),
        '.', source,
        nomass_str, 
        '.png'])

    fig.savefig(fig_file, bbox_inches='tight')
    plt.close()
    return None 

def OrphanSubhaloSMF(type, nsnap_ancestor=20, scatter=0.0, source='li-drory-march'): 
    ''' SMF for orphan central subhalos that do not have ancestors at snapshot 
    nsnap_ancestor.
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # Central subhalo at nsnap_ancestor
    if type == 'central': 
        subh = CentralSubhalos()
    elif type == 'all': 
        subh = Subhalos()
    subh.Read(nsnap_ancestor, scatter=scatter, source=source, nsnap_ancestor=nsnap_ancestor)
    ancestor_index = subh.index
    ancestor_mass = subh.mass

    for i_snap in [1, 5, 10, 15]:
        if i_snap >= nsnap_ancestor: 
            continue 

        smf = SMF()
        if type == 'central': 
            subh = CentralSubhalos()
        elif type == 'all': 
            subh = Subhalos()
        subh.Read(i_snap, scatter=scatter, source=source, nsnap_ancestor=nsnap_ancestor)
    
        # total central subhalo SMF
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        sub.plot(subh_mass, subh_phi, lw=4, c=pretty_colors[i_snap], alpha=0.5) 

        # SMF of central subhalos who's ancestors are centrals at snapshot 20
        if i_snap == 1: 
            label = 'Galaxies that do not have Ancestors'
        else: 
            label = None
        orphan = np.invert(np.in1d(getattr(subh, 'ancestor'+str(nsnap_ancestor)), ancestor_index))
        orph_mass, orph_phi = smf.smf(subh.mass[orphan])
        sub.plot(orph_mass, orph_phi, lw=2, ls='--', c=pretty_colors[i_snap], label=label) 

    anc_mass, anc_phi = smf.smf(ancestor_mass)
    sub.plot(anc_mass, anc_phi,  lw=4, c='gray', label='Total') 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    sub.legend(loc='upper right')

    fig_file = ''.join([
        'figure/test/'
        'OrphanSubhaloSMF', 
        '.', type, '_subhalo',  
        '.scatter', str(round(scatter, 2)), 
        '.ancestor', str(nsnap_ancestor),
        '.', source,
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight')
    plt.close()
    return None

def test_ancestors_without_descendant(nsnap_ancestor = 20, scatter = 0.0, source='li-drory-march'): 
    '''
    What happens to ancestors that do not have descendants at nsnap = 1


    Notes
    -----
    * More than 50% of the subhalos at nsnap=20 do not have descendants at nsnap = 1. 
        What happens to them? 
    * Do they become satellites? --> Some become satellites, others 
        with no mass subhalos don't stay in the catalog at all
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # Central subhalo at nsnap_ancestor (ancesotr)
    anc = CentralSubhalos()
    anc.read(nsnap_ancestor, scatter=scatter, source=source)
    
    # Centrals subhalo at nsnap = 1 (descendant)
    des = CentralSubhalos()
    des.read(1, scatter=scatter, source=source)
    
    no_descendants = np.invert(np.in1d(anc.index, des.ancestor20))
    massive_nodescendants = np.where(anc.mass[no_descendants] > 0.0)
    print 'N_SH no descendants ', len(anc.mass[no_descendants])
    print 'N_SH total ', len(anc.mass)
    print np.float(len(anc.mass[no_descendants]))/np.float(len(anc.mass))

    no_des_index = anc.index[no_descendants][massive_nodescendants][:5]
    print anc.mass[no_descendants][massive_nodescendants][:5]
    for isnap in range(1, nsnap_ancestor)[::-1]: 
        i_des = CentralSubhalos()
        i_des.read(isnap, scatter=scatter, source=source)
    
        #print isnap, np.in1d(no_des_index, i_des.ancestor20)

        #if not np.in1d(no_des_index, i_des.ancestor20)[0]: 
        des_sh = Subhalos()
        des_sh.read(isnap, scatter=scatter, source=source)
        #if np.sum(np.in1d(des_sh.ancestor20, no_des_index)) != len(no_des_index): 
        #    raise ValueError
        print des_sh.ilk[np.in1d(des_sh.ancestor20, no_des_index)][np.argsort(des_sh.ancestor20[np.in1d(des_sh.ancestor20, no_des_index)])]
        print 'M* ', des_sh.mass[np.in1d(des_sh.ancestor20, no_des_index)][np.argsort(des_sh.ancestor20[np.in1d(des_sh.ancestor20, no_des_index)])]
        print 'M_halo ', getattr(des_sh, 'halo.m')[np.in1d(des_sh.ancestor20, no_des_index)][np.argsort(des_sh.ancestor20[np.in1d(des_sh.ancestor20, no_des_index)])]
        #print des_sh.mass[np.in1d(des_sh.index, no_des_index)]
        #np.in1d(i_des.ancestor20, no_des_index)

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([6.0, 12.0])
    sub.legend(loc='upper right')

"""
    #test_orphan_subhalo_smf(scatter=0.0, source='li-march')
    #test_descendant_subhalo_smf(scatter=0.0, source='li-march')
    #test_nomass_descendant_subhalo_smf(scatter=0.0, source='li-march')
    #test_ancestors_without_descendant(scatter=0.0, source='li-march')
"""





if __name__=='__main__': 
    #SubhaloSMF('all', scatter=0.2, source='li-march', 
    #        nsnap_ancestor=15)
    #SubhaloLF(scatter=0.0, source='cool_ages', nsnap_ancestor=11)
    for scat in [0.0]: #[0.2, 0.3, 0.4]: 
        SubhaloLF(scatter=scat, source='cool_ages', nsnap_ancestor=11)
        Subhalo_MhaloMag(scatter=scat, source='cool_ages', nsnap_ancestor=11)

    #for scat in [0.0, 0.2]: 
    #    for type in ['all', 'central']: 
    #        #SubhaloSMF(type, scatter=scat, source='li-march', nsnap_ancestor=20)
    #        #SubhaloSMF(type, scatter=scat, source='li-march', nsnap_ancestor=15)
    #        #SubhaloSMF(type, scatter=scat, source='li-march', nsnap_ancestor=10)
    #        #DescendantSubhaloSMF(type, nsnap_ancestor=20, scatter=scat, source='li-march')
    #        #DescendantSubhaloSMF(type, nsnap_ancestor=15, scatter=scat, source='li-march')
    #        #DescendantSubhaloSMF(type, nsnap_ancestor=10, scatter=scat, source='li-march')
    #        OrphanSubhaloSMF(type, nsnap_ancestor=20, scatter=scat, source='li-march')
    #        OrphanSubhaloSMF(type, nsnap_ancestor=15, scatter=scat, source='li-march')
    #        OrphanSubhaloSMF(type, nsnap_ancestor=10, scatter=scat, source='li-march')
