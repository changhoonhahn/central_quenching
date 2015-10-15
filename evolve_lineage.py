"""

Evolve Lineage CenQue object ancestor and descendants

"""
import time
import numpy as np
    
from cenque import CenQue
from lineage import Lineage
from quiescent_fraction import get_fq
from util.cenque_utility import get_zsnap
from plotting.plot_cenque import PlotCenque
from util.cenque_utility import get_q_ssfr_mean
from util.tau_quenching import get_quenching_efold

def evolve_lineage(
        nsnap_ancestor = 13, 
        nsnap_descendant = 1, 
        tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
        quiet = True):
    """ Evolve lienage 
    """
    start_time = time.time()

    # read in the lineage (< 0.05 seconds)
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor, nsnap_descendant = nsnap_descendant)
    bloodline.readin()

    print 'Reading in the lineage takes', time.time() - start_time, ' seconds'
    start_time = time.time()
    
    ancestor = bloodline.ancestor_cq 
    descendant = bloodline.descendant_cq
    
    if not np.array_equal(ancestor.snap_index, descendant.ancestor13_index): 
        raise ValueError
    
    # initialize SF properties of descendant 
    n_descendant = len(descendant.snap_index)
    descendant.sfr    = np.array([-999. for i in xrange(n_descendant)])
    descendant.ssfr   = np.array([-999. for i in xrange(n_descendant)])
    descendant.q_ssfr = np.array([-999. for i in xrange(n_descendant)])

    #qaplot = PlotCenque()
    #qaplot.cenque_ssfr_dist(ancestor, line_color='b')
    
    q_ancestors = np.where( ancestor.gal_type == 'quiescent' )[0]
    sf_ancestors = np.where( ancestor.gal_type == 'star-forming' )[0]

    # quiescent evolution (ssfr reamins constant)
    descendant.ssfr[q_ancestors] = ancestor.ssfr[q_ancestors]
    descendant.sfr[q_ancestors] = descendant.mass[q_ancestors] + descendant.ssfr[q_ancestors]
    
    print 'Reading in the lineage takes', time.time() - start_time, ' seconds'
    start_time = time.time()
    
    # star-forming evolution 
    # quenching probability (just quiescent fraction)
    P_q = get_fq(
            descendant.mass[sf_ancestors], descendant.zsnap, 
            lit = ancestor.fq_prop['name'])
    np.random.seed()
    randoms = np.random.uniform(0, 1, len(sf_ancestors)) 
    print 'Reading in the lineage takes', time.time() - start_time, ' seconds'
    start_time = time.time()
    
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
    print 'Reading in the lineage takes', time.time() - start_time, ' seconds'
    start_time = time.time()

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

    start_time = time.time()
    #q_ssfr_final
    q_ssfr_mean = get_q_ssfr_mean(descendant.mass[sf_ancestors[is_qing]])
    descendant.q_ssfr[sf_ancestors[is_qing]] = \
            0.18 * np.random.randn(len(is_qing[0])) + q_ssfr_mean 

    print 'Calculating tau takes', time.time() - start_time, ' seconds'
    start_time = time.time()

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

    #qaplot.cenque_ssfr_dist(descendant, line_color='r')
    #qaplot.groupcat_ssfr_dist(Mrcut=18)
    #plt.show()
    
    print time.time() - start_time, ' seconds'
    return None

if __name__=="__main__":
    evolve_lineage(
            nsnap_ancestor = 13, 
            nsnap_descendant = 1, 
            tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
            )

    evolve_lineage(
            nsnap_ancestor = 13, 
            nsnap_descendant = 1, 
            tau_prop = {'name': 'satellite'}, 
            )

