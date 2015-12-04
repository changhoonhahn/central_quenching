import time
import numpy as np

from ssfr import Ssfr
from cenque import CenQue
from lineage import Lineage

import sfr_evol
import mass_evol

from plotting.plot_fq import PlotFq
from quiescent_fraction import cq_fq
from quiescent_fraction import get_fq
from plotting.plot_tau import plot_tau
from plotting.plot_ssfr import PlotSSFR
from util.cenque_utility import get_zsnap
from sfms.fitting import get_param_sfr_mstar_z
from util.cenque_utility import get_q_ssfr_mean
from util.tau_quenching import get_quenching_efold
from group_catalog.group_catalog import central_catalog

from sf_inherit import sf_inherit

def rho_fq_ssfr_descendant(
        nsnap_descendant = 1, 
        nsnap_ancestor = 20, 
        pq_prop = {'slope': 0.0, 'yint':0.0}, 
        tau_prop = {'name': 'instant'}, 
        fq_prop = {'name': 'wetzelsmooth'},
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

    bloodline = sf_inherit(
            [nsnap_descendant], 
            nsnap_ancestor = nsnap_ancestor, 
            pq_prop = pq_prop,
            tau_prop = tau_prop
            )

    descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
    descendant_ssfr = Ssfr()
    descendant_ssfr_bin, descendant_ssfr_hist = descendant_ssfr.cenque(descendant)

    group_ssfr = Ssfr()
    group_ssfr_bin, group_ssfr_hist = group_ssfr.groupcat(Mrcut=Mrcut)

    descendant_masses, descendant_f_q = cq_fq(descendant)
            
    param_f_q = get_fq(
            descendant_masses, 
            descendant.zsnap, 
            lit = fq_prop['name']
            )  
    mass_range = np.where(
            (descendant_masses > 9.5) & 
            (descendant_masses < 12.0)
            )
    l2_f_q = np.sum((descendant_f_q[mass_range] - param_f_q[mass_range])**2/param_f_q[mass_range]) 
            
    l2_ssfr = 0.0
    n_length = 0
    for i_massbin, massbin in enumerate(group_ssfr.mass_bins): 

        if not np.array_equal(descendant_ssfr_bin[i_massbin], group_ssfr_bin[i_massbin]):
            raise ValueError()

        # sSFR comparison range

        q_ssfr_massbin = np.min(get_q_ssfr_mean(massbin)) 

        sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

        sf_ssfr_massbin = sfr_mstar_z(massbin[1], z_med) - massbin[1]

        green_range = np.where(
                (descendant_ssfr_bin[i_massbin] > q_ssfr_massbin) &
                (descendant_ssfr_bin[i_massbin] < sf_ssfr_massbin)
                )

        n_length += len(green_range[0])

        #print np.sum((evol_cq_ssfr_hist[i_massbin][green_range] - group_ssfr_hist[i_massbin][green_range])**2)
        #print "----------------------------------------"
        #print massbin
        #print q_ssfr_massbin, sf_ssfr_massbin
        #print np.sum( 
        #        (descendant_ssfr_hist[i_massbin][green_range] - group_ssfr_hist[i_massbin][green_range])**2/group_ssfr_hist[i_massbin][green_range]
        #        )

        l2_ssfr += np.sum( 
                (descendant_ssfr_hist[i_massbin][green_range] - group_ssfr_hist[i_massbin][green_range])**2/group_ssfr_hist[i_massbin][green_range]
                )
    #l2_ssfr /= np.float(n_length)
    return l2_f_q+l2_ssfr
