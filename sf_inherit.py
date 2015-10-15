"""

Evolve star forming properties of ancestor CenQue objects within the 
Lineage object ancestor for the descendant CenQue objects

"""
import time
import numpy as np

from ssfr import Ssfr
from cenque import CenQue
from lineage import Lineage
from quiescent_fraction import get_fq
from util.cenque_utility import get_zsnap
from plotting.plot_cenque import PlotCenque
from sfms.fitting import get_param_sfr_mstar_z
from util.cenque_utility import get_q_ssfr_mean
from util.tau_quenching import get_quenching_efold

def sf_inherit(
        nsnap_ancestor = 13, 
        nsnap_descendant = 1, 
        tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
        quiet = True, 
        qaplot = False):
    """ For Linage object class, evolve star forming properties of ancestors based on 
    quenching timescale 
    """

    # read in the lineage (< 0.05 seconds)
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.readin()
    
    ancestor = bloodline.ancestor_cq 
    descendant = bloodline.descendant_cq
    
    if not np.array_equal(ancestor.snap_index, descendant.ancestor13_index): 
        raise ValueError
    
    # initialize SF properties of descendant 
    n_descendant = len(descendant.snap_index)
    descendant.sfr    = np.array([-999. for i in xrange(n_descendant)])
    descendant.ssfr   = np.array([-999. for i in xrange(n_descendant)])
    descendant.q_ssfr = np.array([-999. for i in xrange(n_descendant)])
    
    
    q_ancestors = np.where( ancestor.gal_type == 'quiescent' )[0]
    sf_ancestors = np.where( ancestor.gal_type == 'star-forming' )[0]

    # quiescent evolution (ssfr reamins constant)
    descendant.ssfr[q_ancestors] = ancestor.ssfr[q_ancestors]
    descendant.sfr[q_ancestors] = descendant.mass[q_ancestors] + descendant.ssfr[q_ancestors]
    
    # star-forming evolution 
    # quenching probability (just quiescent fraction)
    P_q = get_fq(
            descendant.mass[sf_ancestors], descendant.zsnap, 
            lit = ancestor.fq_prop['name'])
    np.random.seed()
    randoms = np.random.uniform(0, 1, len(sf_ancestors)) 
    
    # which SF galaxies are quenching or not 
    is_qing = np.where(P_q > randoms) 
    is_notqing = np.where(P_q <= randoms)
    
    # star forming decendants that aren't quenching 
    # SFRs all evolve by 0.76 * Delta z
    descendant.sfr[sf_ancestors[is_notqing]] = \
            ancestor.sfr[sf_ancestors[is_notqing]] + \
            0.76 * (descendant.zsnap - ancestor.zsnap) 
    descendant.ssfr[sf_ancestors[is_notqing]] = \
            descendant.sfr[sf_ancestors[is_notqing]] - \
            descendant.mass[sf_ancestors[is_notqing]]

    # star forming decendants that ARE quenching
    # determine time at which the decendant start quenching based on some PDF 
    # test case is uniform 
    t_descendant = descendant.t_cosmic
    t_ancestor = ancestor.t_cosmic
    
    np.random.seed()
    # cosmic time/redshift at which the SF galaxy is quenched
    t_q = np.random.uniform(t_descendant, t_ancestor, len(is_qing[0]))
    z_q = get_zsnap(t_q)

    # quenching e-fold timescale
    tau_q = get_quenching_efold(
            descendant.mass[sf_ancestors[is_qing]], 
            tau_param = tau_prop
            )

    # q_ssfr_final (final quenched SSFR. This is specified due to the fact that 
    # there's a lower bound on the SSFR)
    q_ssfr_mean = get_q_ssfr_mean(descendant.mass[sf_ancestors[is_qing]])
    descendant.q_ssfr[sf_ancestors[is_qing]] = \
            0.18 * np.random.randn(len(is_qing[0])) + q_ssfr_mean 
    
    # apply SFR quenching evolution
    descendant.sfr[sf_ancestors[is_qing]] = \
            ancestor.sfr[sf_ancestors[is_qing]] + \
            0.76 * (descendant.zsnap - z_q) + \
            np.log10( np.exp( (t_ancestor - t_q) / tau_q ) )
    descendant.ssfr[sf_ancestors[is_qing]] = \
            descendant.sfr[sf_ancestors[is_qing]] - \
            descendant.mass[sf_ancestors[is_qing]]
    
    overquenched = np.where(
            descendant.q_ssfr[sf_ancestors[is_qing]] > descendant.ssfr[sf_ancestors[is_qing]]
            )
    descendant.ssfr[sf_ancestors[is_qing[0][overquenched]]] =\
            descendant.q_ssfr[sf_ancestors[is_qing[0][overquenched]]]

    if qaplot: 
        cqplot = PlotCenque()
        cqplot.cenque_ssfr_dist(ancestor, line_color='b')
        cqplot.cenque_ssfr_dist(descendant, line_color='r')
        cqplot.groupcat_ssfr_dist(Mrcut=18)
        plt.show()
    
    return descendant

def ssfr_lineage_evol(
        start_nsnap = 13, 
        final_nsnap = 1, 
        tau_prop = {'name': 'instant'}, 
        **kwargs
        ): 
    """ Calculate sSFR distribution of evolved Lineage
    """

    evolved_cq = evolve_lineage(
            nsnap_ancestor = start_nsnap, 
            nsnap_descendant = final_nsnap, 
            tau_prop = tau_prop
            )
    
    evolved_ssfr = Ssfr()

    return evolved_ssfr.cenque(evolved_cq)

def rho_ssfr_lineage_evol(
        start_nsnap = 13, 
        final_nsnap = 1, 
        tau_prop = {'name': 'instant'}, 
        Mrcut=18, 
        **kwargs
        ): 
    """ Compare sSFR distribution of evolved CenQue and
    SDSS Group Catalog in the green valley
    """

    if Mrcut == 18: 
        z_med = 0.03
    elif Mrcut == 19: 
        z_med = 0.05
    elif Mrcut == 20: 
        z_med = 0.08

    evol_cq_ssfr_bin, evol_cq_ssfr_hist = ssfr_lineage_evol(
            start_nsnap = start_nsnap, 
            final_nsnap = final_nsnap, 
            tau_prop = tau_prop
            )

    group_ssfr = Ssfr()
    group_ssfr_bin, group_ssfr_hist = group_ssfr.groupcat(Mrcut=Mrcut)
    
    l2_ssfr = 0.0
    for i_massbin, massbin in enumerate(group_ssfr.mass_bins): 

        if not np.array_equal(evol_cq_ssfr_bin[i_massbin], group_ssfr_bin[i_massbin]):
            raise ValueError()

        # sSFR comparison range

        q_ssfr_massbin = np.mean(get_q_ssfr_mean(massbin)) 

        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

        sf_ssfr_massbin = sfr_mstar_z(massbin[1], z_med) - massbin[1]

        green_range = np.where(
                (evol_cq_ssfr_bin[i_massbin] > q_ssfr_massbin) &
                (evol_cq_ssfr_bin[i_massbin] < sf_ssfr_massbin)
                )

        #print np.sum((evol_cq_ssfr_hist[i_massbin][green_range] - group_ssfr_hist[i_massbin][green_range])**2)

        l2_ssfr += np.sum((evol_cq_ssfr_hist[i_massbin][green_range] - group_ssfr_hist[i_massbin][green_range])**2)

    return l2_ssfr

if __name__=="__main__":
    print rho_ssfr_lineage_evol( 
            start_nsnap = 13, 
            final_nsnap = 1, 
            tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6},
            Mrcut=18
            ) 
    
    print rho_ssfr_lineage_evol( 
            start_nsnap = 13, 
            final_nsnap = 1, 
            tau_prop = {'name': 'satellite'},
            Mrcut=18
            ) 

    print rho_ssfr_lineage_evol( 
            start_nsnap = 13, 
            final_nsnap = 1, 
            tau_prop = {'name': 'instant'},
            Mrcut=18
            ) 

