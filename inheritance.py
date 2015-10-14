"""

Test parent children match over multiple snapshot steps

Author(s): ChangHoon Hahn

"""
import time
import random 
import numpy as np
import warnings

from cenque import CenQue
from util import cenque_utility as util
from assign_sfr import assign_sfr
from quiescent_fraction import get_fq
from quiescent_fraction import get_fq_nsnap
from util.gal_classify import sfr_cut 
from util.gal_classify import sfq_classify
from sfms.fitting import get_param_sfr_mstar_z
from util.tau_quenching import get_quenching_efold

def evolve_cq( 
        nsnap_0 = 13, 
        nsnap_f = 1,
        quiet = True, 
        **kwargs
        ): 
    """ Evolve SF properties and stellar mass properties of galaxies 
    in the CenQue class object over cosmic time until final_nsnap.
    
    ----------------------------------------------------------------
    Parameters
    ----------------------------------------------------------------
    nsnap_0 : first snapshot #
    nsnap_f : final snapshot #
    quiet : yak or not  

    ----------------------------------------------------------------
    Notes
    ----------------------------------------------------------------
    * Quenching of SF galaxies now involve quenching-efold *and* 
    overall SF MS SFR decrease

    """
    
    if not quiet: 
        evo_start_time = time.time()

    original_cq = CenQue() 
    original_cq.import_treepm(nsnap_0) 

    parent_cq = original_cq

    for i_snap in range(nsnap_f, nsnap_0)[::-1]:    
    
        if not quiet:
            print ''
            print '---------------------------------------'
            print 'Evolving to ', str(i_snap) 
    
        child_cq = CenQue() 
        child_cq.import_treepm(i_snap) 

        child_cq = evolve_onestep(parent_cq, child_cq, nsnap_0 = nsnap_0, quiet=quiet)
        
        parent_cq = child_cq
    
    if not quiet: 
        print 'Total Evolution takes ', time.time() - evo_start_time
    has_ancestor = np.where(child_cq.ancestor13_index >= 0)
    ancestor, descendant = parent_children_match(original_cq.snap_index, child_cq.ancestor13_index[has_ancestor])

    print original_cq.snap_index[ancestor]
    print child_cq.ancestor13_index[has_ancestor[0][descendant]]
    print original_cq.mass[ancestor]
    print child_cq.mass[has_ancestor[0][descendant]]
    return child_cq 

def evolve_onestep(parent_cq, child_cq, nsnap_0 = 13, quiet=False):
    """ Evolve CenQue class object by one time step 
    """
    evo_start = time.time()

    mass_bins = child_cq.mass_bins
    
    # remove galaxies below the min and max mass
    within_massbin = np.where((child_cq.mass > min(mass_bins.mass_low))) 
    child_cq.sample_trim(within_massbin)                   

    n_child = len(within_massbin[0])     
    n_parent = len(parent_cq.mass)
    if not quiet: 
        print n_parent, ' parents ', n_child, ' children'
        
    setattr(child_cq, 'ancestor'+str(nsnap_0)+'_index', np.array([-999 for i in xrange(n_child)]))

    # parent children match which assigns indices into dictionaries and then get dictionary values
    parents, children = parent_children_match(parent_cq.snap_index, child_cq.parent)

    # Deal with parent to children inheritance (children inherit the parents' attributes)
    # This step needs to be thought out more. 
    if not quiet: 
        inheritance_time = time.time()

    try: 
        ancestor_index = getattr(parent_cq, 'ancestor'+str(nsnap_0)+'_index')[parents]
        getattr(child_cq, 'ancestor'+str(nsnap_0)+'_index')[children] = ancestor_index
        if not quiet: 
            print ''
            print child_cq.nsnap, '---------' 
            print 'Parent with ancestors ', len(np.where(parent_cq.ancestor13_index >= 0)[0]) 
            print 'Children with ancestors ', len(np.where(child_cq.ancestor13_index >= 0)[0])
            print len(np.where(parent_cq.ancestor13_index >= 0)[0]) -\
                    len(np.where(child_cq.ancestor13_index >= 0)[0]), ' ancestors lost'
    except AttributeError: 
        ancestor_index = parent_cq.snap_index[parents]
        getattr(child_cq, 'ancestor'+str(nsnap_0)+'_index')[children] = ancestor_index

    return child_cq

def parent_children_match(parent_snap_index, child_parent_snap_index):
    """ Match snapshot index of parents to the parent snapshot index of
    childrens by matching (index, parent), (index, child) key pairs 
    Takes approximately < 1 second
    """
    parent_sort_index = np.argsort(parent_snap_index)
    child_sort_index = np.argsort(child_parent_snap_index)
    sorted_parent_snap_index = parent_snap_index[parent_sort_index]
    sorted_child_parent_snap_index = child_parent_snap_index[child_sort_index]

    parent_child_in1d = np.in1d(sorted_parent_snap_index, sorted_child_parent_snap_index)
    child_parent_in1d = np.in1d(sorted_child_parent_snap_index, sorted_parent_snap_index)
    
    parents = parent_sort_index[parent_child_in1d]
    children = child_sort_index[child_parent_in1d]
    #print 'Parent Child match takes ', time.time() - parent_child_start_time, ' seconds'

    return parents, children 

if __name__=='__main__': 
    evolve_cq( 
        nsnap_0 = 13, 
        nsnap_f = 1
        )
