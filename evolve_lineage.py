"""

Evolve Lineage CenQue object ancestor and descendants

"""
import time
import numpy as np
    
from cenque import CenQue
from lineage import Lineage
from quiescent_fraction import get_fq
from plotting.plot_cenque import PlotCenque

def evolve_lineage(nsnap_ancestor = 13, nsnap_descendant = 1, quiet = True):
    """ Evolve lienage 
    """
    start_time = time.time()

    # read in the lineage
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor, nsnap_descendant = nsnap_descendant)
    bloodline.readin()
    
    ancestor = bloodline.ancestor_cq 
    descendant = bloodline.descendant_cq

    n_ancestor = len(ancestor.snap_index)
    n_descendant = len(descendant.snap_index)
    #print n_ancestor, n_descendant
    descendant.sfr  = np.array([-999. for i in xrange(n_descendant)])
    descendant.ssfr = np.array([-999. for i in xrange(n_descendant)])

    qaplot = PlotCenque()
    qaplot.cenque_ssfr_dist(ancestor, line_color='k')
    
    q_ancestors = np.where( ancestor.gal_type == 'quiescent' ) 
    sf_ancestors = np.where( ancestor.gal_type == 'star-forming' ) 

    print 'quiescent galaxy ', len(q_ancestors[0])
    print 'star forming galaxy ', len(sf_ancestors[0])

    print ancestor.snap_index
    print descendant.ancestor13_index

    # quiescent evolution (ssfr reamins constant)
    descendant.ssfr[q_ancestors] = ancestor.ssfr[q_ancestors]
    descendant.sfr[q_ancestors] = descendant.mass[q_ancestors] + descendant.ssfr[q_ancestors]
    
    # star-forming evolution 
    # quenching probabiliyt 
    P_q = get_fq(descendant.mass[sf_ancestors], descendant.zsnap, lit = ancestor.fq_prop['name'])
    print P_q
    
    descendant.sfr[sf_ancestors] = ancestor.sfr[sf_ancestors] + 0.76 * (ancestor.zsnap - descendant.zsnap) 
    descendant.ssfr[sf_ancestors] = descendant.sfr[sf_ancestors] - descendant.mass[sf_ancestors]
    
    qaplot.cenque_ssfr_dist(descendant, line_color='r', line_style='--')
    plt.show()
    
    #print time.time() - start_time, ' seconds'

if __name__=="__main__":
    evolve_lineage(nsnap_ancestor = 13, nsnap_descendant = 1)

