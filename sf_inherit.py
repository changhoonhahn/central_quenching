"""

Evolve star forming properties of ancestor CenQue objects within the 
Lineage object ancestor for the descendant CenQue objects

"""
import time
import numpy as np
from scipy import interpolate

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

def sf_inherit(nsnap_descendants, 
        nsnap_ancestor = 13, 
        pq_prop = {'slope': 0.05, 'yint': 0.0}, 
        tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
        sfrevol_prop = {'name': 'notperiodic'}, 
        sfrevol_massdep = False,
        massevol_prop = {'name': 'sham'}, 
        scatter = 0.0, 
        quiet = True, 
        qaplot = False,
        qaplotname = None, 
        clobber = False):
    """ 
    Evolve star forming properties of ancestor CenQue object in Lineage Class 
    to descendant CenQue object
    """
    if sfrevol_massdep: 
        print "SFR(M*(t), t)"
    else: 
        print "SFR(M*(t0), t)"
    
    # make sure that snapshot = 1 is included among imported descendants
    # and the first element of the list
    if 1 not in nsnap_descendants: 
        nsnap_descendants = [1] + nsnap_descendants
    elif nsnap_descendants[0] != 1: 
        nsnap_descendants.pop(1)
        nsnap_descendants = [1] + nsnap_descendants
    
    # spline between z and t_cosmic
    z, t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
    z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 
    
    logsfr_mstar_z, sig_logsfr_mstar_z = get_param_sfr_mstar_z()

    # read in the lineage (< 0.05 seconds for one snapshot)
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.readin(nsnap_descendants, scatter = scatter, clobber = clobber)

    ancestor = bloodline.ancestor_cq    # ancestor CenQue object
    t_ancestor = ancestor.t_cosmic
    z_ancestor = ancestor.zsnap
    
    q_ancestors = np.where(ancestor.gal_type == 'quiescent')[0]
    sf_ancestors = np.where(ancestor.gal_type == 'star-forming')[0]
        
    if qaplot: 
        # plot ancestor SSFR distribution 
        cqplot = PlotSSFR()
        cqplot.cenque_ssfr_dist(ancestor, line_color='b')

        fqplot = PlotFq()
        fqplot.cenque_fq(ancestor, line_color='k', label = 't = ' + str(round(ancestor.t_cosmic,2)) )
        fqplot.param_fq(line_color = 'k', label = None)
    
    for nsnap_descendant in nsnap_descendants:      # Descendants 
    
        # descendent CenQue object
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
        t_descendant = descendant.t_cosmic 
        z_descendant = descendant.zsnap

        if not np.array_equal(ancestor.snap_index, getattr(descendant, 'ancestor'+str(nsnap_ancestor)+'_index')): 
            # check that lineage is properly tracked
            raise ValueError
    
        # initialize SF properties of descendant 
        n_descendant = len(descendant.snap_index)
        descendant.sfr    = np.array([-999. for i in xrange(n_descendant)])
        descendant.ssfr   = np.array([-999. for i in xrange(n_descendant)])
        descendant.q_ssfr = np.array([-999. for i in xrange(n_descendant)])
    
        # --------------------------------------------------------------------------------
        # quiescent evolution (ssfr reamins constant)
        # --------------------------------------------------------------------------------
        descendant.ssfr[q_ancestors] = ancestor.ssfr[q_ancestors]
        
        def logsfr_quiescent(logmass, t_input):
            return ancestor.ssfr[q_ancestors] + logmass
    
        if massevol_prop['name'] == 'integrated': 
            descendant.mass[q_ancestors], descendant.sfr[q_ancestors] = mass_evol.integrated_rk4(
                    logsfr_quiescent, 
                    ancestor.mass[q_ancestors], 
                    t_ancestor, 
                    t_descendant, 
                    f_retain=massevol_prop['f_retain'], 
                    delt = massevol_prop['t_step']
                    )
        elif massevol_prop['name'] == 'sham': 
            descendant.sfr[q_ancestors] = logsfr_quiescent(
                    descendant.mass[q_ancestors], 
                    t_descendant
                    )

        if nsnap_descendant == 1: 

            # star-forming evolution 
            # quenching probability (just quiescent fraction)
            # P_Q = ( f_Q(zf) - f_Q(z0) ) / (1 - f_Q(z0))
            P_q_offset = pq_prop['slope'] * (descendant.mass[sf_ancestors] - 9.5) + pq_prop['yint']
            P_q = (
                    get_fq(
                        descendant.mass[sf_ancestors], 
                        descendant.zsnap, 
                        lit = ancestor.fq_prop['name']
                        ) - 
                    get_fq(
                        ancestor.mass[sf_ancestors], 
                        ancestor.zsnap, 
                        lit = ancestor.fq_prop['name']
                        )
                    ) / (1.0 -  get_fq( ancestor.mass[sf_ancestors], ancestor.zsnap, lit = ancestor.fq_prop['name']) ) + P_q_offset

            np.random.seed()
            randoms = np.random.uniform(0, 1, len(sf_ancestors)) 
        
            # determine which SF galaxies are quenching or not 
            is_qing = np.where(P_q > randoms) 
            is_notqing = np.where(P_q <= randoms)
            
            # Initialize SFR evolution parameters
            sfrevol_param = sfr_evol.get_sfrevol_param(len(descendant.mass), sf_ancestors, **sfrevol_prop)

            np.random.seed()
            # cosmic time/redshift at which the SF galaxy is quenched 
            # is sampled uniformly between t_ancestor and t_nsnap=1
            t_q = np.random.uniform(t_descendant, t_ancestor, len(sf_ancestors))
            z_q = get_zsnap(t_q)
            t_q[is_notqing] = 999.0     # galaxies that will never quench
            z_q[is_notqing] = -999.0
        
            # quenching e-fold timescale
            tau_q = get_quenching_efold(
                    descendant.mass[sf_ancestors], 
                    tau_param = tau_prop
                    )
        
            # final quenched SSFR. This is specified due to the fact that 
            # there's a lower bound on the SSFR
            q_ssfr_mean = get_q_ssfr_mean(descendant.mass[sf_ancestors])
            final_q_ssfr = 0.18 * np.random.randn(len(sf_ancestors)) + q_ssfr_mean 
        
        # star forming decendants that will quench by Snapshot 1
        descendant.q_ssfr[sf_ancestors] = final_q_ssfr
        q_started = np.where(t_q <= t_descendant)   # SF galaxies that have started quenching
        q_notstarted = np.where(t_q > t_descendant) # SF galaxies taht have NOT started quenching
        
        # ---------------------------------------------------------------------------------------------------
        # star forming evolution 
        # ---------------------------------------------------------------------------------------------------
        # SFR evolution parameters
        sfrevol_param_sf = []
        for i_p in xrange(len(sfrevol_param)):
            sfrevol_param_sf.append(sfrevol_param[i_p][sf_ancestors])

        def logsfr_m_t(logmass, t_input): 

            if sfrevol_massdep:  
                # SFR evolution based on solving an ODE of SFR
                avglogsfr = logsfr_mstar_z(logmass, z_ancestor)
            else: 
                # SFR evolution based on just time evolution 
                avglogsfr = ancestor.avg_sfr[sf_ancestors]

            # log(SFR)_SFMS evolutionfrom t0 to tQ
            logsfr_sfms = sfr_evol.logsfr_sfms_evol(
                    z_ancestor, 
                    z_of_t(t_input),
                    z_q = z_q
                    )
            
            # log(SFR)_duty evolution from t0 to tQ
            logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                    t_ancestor, 
                    t_input, 
                    t_q = t_q, 
                    delta_sfr=ancestor.delta_sfr[sf_ancestors],
                    sfrevol_param=sfrevol_param_sf, 
                    **sfrevol_prop
                    )

            logsfr_quench = sfr_evol.logsfr_quenching(
                    t_q, 
                    t_input, 
                    tau=tau_q
                    )

            return avglogsfr + logsfr_sfms + logsfr_sfduty + logsfr_quench
         
        print 'SHAM Masses : ', descendant.mass[sf_ancestors]

        if massevol_prop['name'] == 'integrated': 

            start_time = time.time()
            descendant.mass[sf_ancestors], descendant.sfr[sf_ancestors] = \
                    mass_evol.integrated_rk4(
                            logsfr_m_t, 
                            ancestor.mass[sf_ancestors], 
                            t_ancestor, 
                            t_descendant, 
                            f_retain = massevol_prop['f_retain'], 
                            delt = massevol_prop['t_step']
                            )
            print 'integration time ', (time.time() - start_time)
        
            print 'Integrated Masses : ', descendant.mass[sf_ancestors]
        
            #print 'blah', descendant.sfr[sf_ancestors]
            #print 'dlah', logsfr_m_t(descendant.mass[sf_ancestors], t_descendant)

        elif massevol_prop['name'] == 'sham': 
            if sfrevol_massdep:
                raise ValueError("You can't both have SHAM masses SFR(M*(t),t) dependence!")

            descendant.sfr[sf_ancestors] = logsfr_m_t(
                    descendant.mass[sf_ancestors], 
                    t_descendant
                    )

        descendant.ssfr[sf_ancestors] = descendant.sfr[sf_ancestors] - descendant.mass[sf_ancestors]
        
        # ---------------------------------------------------------------------------------------------------
        # Account for over quenching  
        # ---------------------------------------------------------------------------------------------------
        overquenched = np.where(
                descendant.q_ssfr[sf_ancestors] > descendant.ssfr[sf_ancestors]
                )
        descendant.ssfr[sf_ancestors[overquenched]] = descendant.q_ssfr[sf_ancestors[overquenched]]

        setattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant), descendant)

        if qaplot: 
            if nsnap_descendant != 1: 
                cqplot.cenque_ssfr_dist(descendant, lw=2, line_color=cqplot.pretty_colors[nsnap_descendant], label=None)

                fqplot.cenque_fq(
                        descendant, 
                        lw=2, 
                        label = 't=' + str(round(descendant.t_cosmic, 1))+',z='+str(round(descendant.zsnap, 2))  
                        )
                fqplot.param_fq(lw=2, line_color = fqplot.pretty_colors[descendant.nsnap], label = None)
            else:
                cqplot.cenque_ssfr_dist(descendant, line_color='r')
                fqplot.cenque_fq(
                        descendant, 
                        label = 't=' + str(round(descendant.t_cosmic, 1))+',z='+str(round(descendant.zsnap, 2))  
                        )
                fqplot.param_fq(
                        line_color = fqplot.pretty_colors[descendant.nsnap], 
                        label = None
                        )

    if qaplot: 
        cqplot.groupcat_ssfr_dist(Mrcut=18)
        cqplot.set_axes()

        gc = central_catalog(Mrcut=18, clobber=False)
        groupcat = CenQue()
        groupcat.sfr = gc.sfr
        groupcat.mass = gc.mass
        groupcat.zsnap = np.median(gc.z)
        gc_mass, gc_fq = cq_fq(groupcat)
        
        fqplot.subs.plot(
                gc_mass, 
                gc_fq, 
                color='red', 
                lw=5, 
                label=r'Group Catalog $\mathtt{M_r = 18}$' 
                ) 

        fqplot.set_axes()
        tauplot = plot_tau([tau_prop])

        if qaplotname is not None: 
            if '.png' not in qaplotname: 
                qaplotname += '.png'
            cqplot.fig.savefig('figure/cq_'+qaplotname, bbox_inches = 'tight')
            fqplot.fig.savefig('figure/fq_'+qaplotname, bbox_inches = 'tight')
            tauplot.savefig('figure/tau_'+qaplotname, bbox_inches = 'tight')
        else: 
            plt.show()
    
    return bloodline 

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
