'''

sSFR distribution of galaxy catalogs

'''

import warnings
import numpy as np

# -- Local --
from cenque import CenQue
from group_catalog.group_catalog import central_catalog

class Ssfr(object): 

    def __init__(self, **kwargs): 
        """ Class that describes the sSFR distribution of a CenQue object
        """
        self.kwargs = kwargs

        self.mass_bins = self.get_mass_bins()
        
        self.ssfr_hist = None
        self.ssfr_bin_low = None 
        self.ssfr_bin_mid = None 
        self.ssfr_bin_high = None 

    def cenque(self, cenque): 
        """ Calculate sSFR distribution from CenQue object

        ssfr_bin_edges are preset 
        """

        self.ssfr_hist = [] 
        self.ssfr_bin_low = []
        self.ssfr_bin_mid = []
        self.ssfr_bin_high = []

        for mass_bin in self.mass_bins: 

            mass_limit = np.where(
                    (cenque.mass >= mass_bin[0]) & 
                    (cenque.mass < mass_bin[1])
                    )

            ngal_bin = len(mass_limit[0])

            # calculate SSFR distribution  
            ssfr_dist, ssfr_bin_edges = np.histogram(
                    cenque.ssfr[mass_limit], 
                    range = [-13.0, -7], 
                    bins = 40, 
                    normed = True
                    )

            self.ssfr_hist.append(ssfr_dist)

            self.ssfr_bin_low.append(ssfr_bin_edges[:-1])
            self.ssfr_bin_high.append(ssfr_bin_edges[1:])
            self.ssfr_bin_mid.append(0.5 * (ssfr_bin_edges[:-1] + ssfr_bin_edges[1:]))
        
        return [self.ssfr_bin_mid, self.ssfr_hist]

    def groupcat(self, Mrcut=18): 
        """ Calculate sSFR distribution from SDSS Group Catalog
        """
        
        groupcat = central_catalog(Mrcut=Mrcut)
        
        if self.ssfr_hist != None: 
            warnings.warn("Overwriting sSFR distribution")

        self.ssfr_hist = [] 
        self.ssfr_bin_low = []
        self.ssfr_bin_mid = []
        self.ssfr_bin_high = []

        for mass_bin in self.mass_bins: 

            mass_limit = np.where(
                    (groupcat.mass >= mass_bin[0]) & 
                    (groupcat.mass < mass_bin[1])
                    )

            ngal_bin = len(mass_limit[0])

            # calculate SSFR distribution  
            ssfr_dist, ssfr_bin_edges = np.histogram(
                    groupcat.ssfr[mass_limit], 
                    range = [-13.0, -7], 
                    bins = 40, 
                    normed = True
                    )

            self.ssfr_hist.append(ssfr_dist)

            self.ssfr_bin_low.append(ssfr_bin_edges[:-1])
            self.ssfr_bin_high.append(ssfr_bin_edges[1:])
            self.ssfr_bin_mid.append(0.5 * (ssfr_bin_edges[:-1] + ssfr_bin_edges[1:]))
        
        return [self.ssfr_bin_mid, self.ssfr_hist]

    def get_mass_bins(self): 
        """ Mass bins where sSFR distribution will be calculated 

        Currently preset
        """

        mass_bins = [ 
                [9.7, 10.1], [10.1, 10.5], [10.5, 10.9], [10.9, 11.3]
                ]

        return mass_bins


def ssfr_cq_evol(
        start_nsnap = 13, 
        final_nsnap = 1, 
        sf_prop = {'name': 'average'}, 
        fq_prop = {'name': 'wetzelsmooth'}, 
        tau_prop = {'name': 'instant'}, 
        mass_evol = 'sham', 
        **kwargs
        ): 
    """ Calculate sSFR distribution of evolved CenQue 
    """

    evolved_cq = CenQue(n_snap = final_nsnap, cenque_type = 'evol_from'+str(start_nsnap))
    evolved_cq.sf_prop = sf_prop
    evolved_cq.fq_prop = fq_prop
    evolved_cq.tau_prop = tau_prop
    evolved_cq.mass_evol = mass_evol
    
    if not os.path.isfile(evolved_cq.file()):
        start_cq = CenQue(n_snap = start_nsnap, cenque_type = 'sf_assigned')
        start_cq.readin()

        evolved_cq = evolve_cq(
                start_cq, 
                final_nsnap = final_nsnap, 
                sf_prop = sf_prop, 
                fq_prop = fq_prop, 
                tau_prop = tau_prop, 
                mass_evol = mass_evol, 
                **kwargs
                )
    else: 
        evolved_cq.readin()

    evolved_cq_ssfr = Ssfr()

    return evolved_cq_ssfr.cenque(evolved_cq)
