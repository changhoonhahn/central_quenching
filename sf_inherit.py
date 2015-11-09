"""

Evolve star forming properties of ancestor CenQue objects within the 
Lineage object ancestor for the descendant CenQue objects

"""
import time
import numpy as np

from ssfr import Ssfr
from cenque import CenQue
from lineage import Lineage
# SFR evolution 
import sfr_evol

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
        scatter = 0.0, 
        quiet = True, 
        qaplot = False,
        qaplotname = None, 
        clobber = False):
    """ 
    Evolve star forming properties of ancestor CenQue object in Lineage Class 
    to descendant CenQue object
    """
    
    # make sure that snapshot = 1 is included among imported descendants
    # and the first element of the list
    if 1 not in nsnap_descendants: 
        nsnap_descendants = [1] + nsnap_descendants
    elif nsnap_descendants[0] != 1: 
        nsnap_descendants.pop(1)
        nsnap_descendants = [1] + nsnap_descendants

    # read in the lineage (< 0.05 seconds for one snapshot)
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.readin(nsnap_descendants, scatter = scatter, clobber = clobber)

    ancestor = bloodline.ancestor_cq    # ancestor CenQue object
    t_ancestor = ancestor.t_cosmic
    
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

        if not np.array_equal(ancestor.snap_index, getattr(descendant, 'ancestor'+str(nsnap_ancestor)+'_index')): 
            # check that lineage is properly tracked
            raise ValueError
    
        # initialize SF properties of descendant 
        n_descendant = len(descendant.snap_index)
        descendant.sfr    = np.array([-999. for i in xrange(n_descendant)])
        descendant.ssfr   = np.array([-999. for i in xrange(n_descendant)])
        descendant.q_ssfr = np.array([-999. for i in xrange(n_descendant)])
    
        # quiescent evolution (ssfr reamins constant)
        descendant.ssfr[q_ancestors] = ancestor.ssfr[q_ancestors]
        descendant.sfr[q_ancestors] = descendant.mass[q_ancestors] + descendant.ssfr[q_ancestors]
    
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

            ### Initialize SFR(t_cosmic) function
            ### 
            ### 
            ### 
            ### 
            ### 
            ### 

            np.random.seed()
            # cosmic time/redshift at which the SF galaxy is quenched 
            # is sampled uniformly between t_ancestor and t_nsnap=1
            t_q = np.random.uniform(t_descendant, t_ancestor, len(is_qing[0]))
            z_q = get_zsnap(t_q)
        
            # quenching e-fold timescale
            tau_q = get_quenching_efold(
                    descendant.mass[sf_ancestors[is_qing]], 
                    tau_param = tau_prop
                    )
        
            # final quenched SSFR. This is specified due to the fact that 
            # there's a lower bound on the SSFR
            q_ssfr_mean = get_q_ssfr_mean(descendant.mass[sf_ancestors[is_qing]])
            final_q_ssfr = 0.18 * np.random.randn(len(is_qing[0])) + q_ssfr_mean 
        
        # star forming decendants that never quench 
        # SFRs all evolve by 0.76 * Delta z

        sfrevol_param_nq = []  
        for i_p in xrange(len(sfrevol_param)):
            sfrevol_param_nq.append(sfrevol_param[i_p][sf_ancestors[is_notqing]])

        def logsfr_nq(logmass, t_input): 
            # this is a simplification. logsfr should depend on the new mass at each time step 
            # but that complicates the SF duty cycle 
            # avglogsfr = logsfr_mstar_z(logmass, t_ancestor)    
            avglogsfr = ancestor.avg_sfr[sf_ancestors[is_notqing]]

            logsfr_sfms = sfr_evol.logsfr_sfms_evol(t_ancestor, t_input)

            logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                    t_ancestor, 
                    t_input, 
                    delta_sfr=ancestor.delta_sfr[sf_ancestors[is_notqing]],
                    sfrevol_param=sfrevol_param_nq, 
                    **sfrevol_prop
                    )

            return avglogsfr + logsfr_sfms + logsfr_sfduty
        
        descendant.sfr[sf_ancestors[is_notqing]] = \
                sfr_evol.sfr_evol(
                        t_cosmic = t_descendant - t_ancestor, 
                        indices = sf_ancestors[is_notqing], 
                        sfrevol_param = sfrevol_param, 
                        ancestor_sfr = ancestor.sfr, 
                        ancestor_delta_sfr = ancestor.delta_sfr, 
                        **sfrevol_prop
                        ) + \
                                0.76 * (descendant.zsnap - ancestor.zsnap) 
        print ''
        print 'Not Quenching'
        print len(descendant.sfr[sf_ancestors[is_notqing]])
        print len(logsfr_nq(ancestor.mass, t_descendant))
        print np.min(descendant.sfr[sf_ancestors[is_notqing]] - logsfr_nq(ancestor.mass, t_descendant))
        print np.max(descendant.sfr[sf_ancestors[is_notqing]] - logsfr_nq(ancestor.mass, t_descendant))
        print ''
        print ''


        descendant.ssfr[sf_ancestors[is_notqing]] = \
                descendant.sfr[sf_ancestors[is_notqing]] - \
                descendant.mass[sf_ancestors[is_notqing]]

        # star forming decendants that ARE quenching
        descendant.q_ssfr[sf_ancestors[is_qing]] = final_q_ssfr
         
        q_started = np.where(t_q <= t_descendant) 
        q_notstarted = np.where(t_q > t_descendant)
        
        sfrevol_param_qs = []  
        for i_p in xrange(len(sfrevol_param)):
            sfrevol_param_qs.append(sfrevol_param[i_p][sf_ancestors[is_qing[0][q_started]]])

        def logsfr_qs(logmass, t_input): 
            # this is a simplification. logsfr should depend on the new mass at each time step 
            # but that complicates the SF duty cycle 
            # avglogsfr = logsfr_mstar_z(logmass, t_ancestor)    
            avglogsfr = ancestor.avg_sfr[sf_ancestors[is_qing[0][q_started]]]

            logsfr_sfms = sfr_evol.logsfr_sfms_evol(t_ancestor, t_q[q_started])

            logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                    t_ancestor, 
                    t_q[q_started], 
                    delta_sfr=ancestor.delta_sfr[sf_ancestors[is_qing[0][q_started]]],
                    sfrevol_param=sfrevol_param_qs, 
                    **sfrevol_prop
                    )
            
            logsfr_quench = sfr_evol.logsfr_quenching(
                    t_q[q_started], 
                    t_input, 
                    tau=tau_q[q_started]
                    )

            return avglogsfr + logsfr_sfms + logsfr_sfduty + logsfr_quench
        
        descendant.sfr[sf_ancestors[is_qing[0][q_started]]] = \
                sfr_evol.sfr_evol(t_cosmic = t_q[q_started] - t_ancestor, 
                        indices = sf_ancestors[is_qing[0][q_started]],
                        sfrevol_param = sfrevol_param, 
                        ancestor_sfr = ancestor.sfr, 
                        ancestor_delta_sfr = ancestor.delta_sfr, 
                        **sfrevol_prop) + \
                                0.76 * (z_q[q_started] - ancestor.zsnap) + \
                                np.log10( np.exp( (t_q[q_started] - t_descendant) / tau_q[q_started] ) )

        print ''
        print 'Quenching Started'
        print len(descendant.sfr[sf_ancestors[is_qing[0][q_started]]])
        print len(logsfr_qs(ancestor.mass, t_descendant))
        print np.min(descendant.sfr[sf_ancestors[is_qing[0][q_started]]]-logsfr_qs(ancestor.mass, t_descendant))
        print np.max(descendant.sfr[sf_ancestors[is_qing[0][q_started]]]-logsfr_qs(ancestor.mass, t_descendant))
        print ''
        print ''

        if len(q_notstarted[0]) > 0:

            sfrevol_param_qns = []  
            for i_p in xrange(len(sfrevol_param)):
                sfrevol_param_qns.append(sfrevol_param[i_p][sf_ancestors[is_qing[0][q_notstarted]]])

            def logsfr_qns(logmass, t_input): 
                # this is a simplification. logsfr should depend on the new mass at each time step 
                # but that complicates the SF duty cycle 
                # avglogsfr = logsfr_mstar_z(logmass, t_ancestor)    
                avglogsfr = ancestor.avg_sfr[sf_ancestors[is_qing[0][q_notstarted]]]

                logsfr_sfms = sfr_evol.logsfr_sfms_evol(t_ancestor, t_descendant)

                logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                        t_ancestor, 
                        t_descendant, 
                        delta_sfr=ancestor.delta_sfr[sf_ancestors[is_qing[0][q_notstarted]]],
                        sfrevol_param=sfrevol_param_qns, 
                        **sfrevol_prop
                        )

                return avglogsfr + logsfr_sfms + logsfr_sfduty

            descendant.sfr[sf_ancestors[is_qing[0][q_notstarted]]] = \
                    sfr_evol.sfr_evol(
                            t_cosmic = t_descendant - t_ancestor, 
                            indices = sf_ancestors[is_qing[0][q_notstarted]],
                            sfrevol_param = sfrevol_param, 
                            ancestor_sfr = ancestor.sfr, 
                            ancestor_delta_sfr = ancestor.delta_sfr, 
                            **sfrevol_prop
                            ) + \
                    0.76 * (descendant.zsnap - ancestor.zsnap)
            print ''
            print 'Quenching Not Started'
            print len(descendant.sfr[sf_ancestors[is_qing[0][q_notstarted]]])
            print len(logsfr_qns(ancestor.mass, t_descendant))
            print np.min(descendant.sfr[sf_ancestors[is_qing[0][q_notstarted]]] - logsfr_qns(ancestor.mass, t_descendant))
            print np.max(descendant.sfr[sf_ancestors[is_qing[0][q_notstarted]]] - logsfr_qns(ancestor.mass, t_descendant))
            print ''
            print ''

        descendant.ssfr[sf_ancestors[is_qing]] = descendant.sfr[sf_ancestors[is_qing]] - descendant.mass[sf_ancestors[is_qing]]
        
        overquenched = np.where(
                descendant.q_ssfr[sf_ancestors[is_qing]] > descendant.ssfr[sf_ancestors[is_qing]]
                )
        descendant.ssfr[sf_ancestors[is_qing[0][overquenched]]] =\
                descendant.q_ssfr[sf_ancestors[is_qing[0][overquenched]]]

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

if __name__=="__main__":
    pass
    #print rho_fq_ssfr_descendant(
    #            nsnap_descendant = 1, 
    #            nsnap_ancestor = 20, 
    #            pq_prop = {'slope': 0.1, 'yint': 0.0054}, 
    #            tau_prop = {'name': 'line', 'fid_mass': 11.2, 'slope': -0.68, 'yint': 0.6},
    #            Mrcut=18
    #            )
    
    #sf_inherit(
    #        [1, 2, 3, 6, 9, 12, 15, 18], 
    #        nsnap_ancestor = 20, 
    #        pq_prop = {'slope': 0.0, 'yint': 0.0}, 
    #        tau_prop = {'name': 'instant'}, 
    #        qaplot = True,
    #        scatter = 0.2 
    #        )

    #tau_prop = {'name': 'line', 'fid_mass': 11.2, 'slope': -0.68, 'yint': 0.6},

    #print 'terrible fit both f_q and ssfr' 
    #print rho_fq_ssfr_descendant(
    #            nsnap_descendant = 1, 
    #            nsnap_ancestor = 20, 
    #            pq_prop = {'slope': 0.25, 'yint': 0.25}, 
    #            tau_prop = {'name': 'line', 'fid_mass': 10.5, 'slope': -0.9, 'yint': 1.0},
    #            Mrcut=18)
    
    #print 'decent fit both f_q and ssfr' 
    #print rho_fq_ssfr_descendant(
    #            nsnap_descendant = 1, 
    #            nsnap_ancestor = 20, 
    #            pq_prop = {'slope': 0.0, 'yint': 0.0}, 
    #            tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6},
    #            Mrcut=18)
    
    #sf_inherit(
    #        [1, 2, 3, 6, 9, 12, 15, 18], 
    #        nsnap_ancestor = 20, 
    #        pq_prop = {'slope': 0.05, 'yint': 0.0}, 
    #        tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6},
    #        qaplot=True,
    #        qaplotname='decentfit_example'
    #        )
    
    #sf_inherit(
    #        [1, 2, 3, 6, 9, 12, 15, 18], 
    #        nsnap_ancestor = 20, 
    #        pq_prop = {'slope': 0.05, 'yint': 0.1}, 
    #        tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6},
    #        qaplot=True,
    #        qaplotname='decentfit_example'
    #        )

    """
    fqs, ssfrs = [], [] 
    for yint in [0.0, 0.05, 0.1]: 
        for slope in [ 0.0, 0.05, 0.1]: 
            #print slope, yint 
            fq, ssfr = rho_fq_ssfr_descendant(
                nsnap_descendant = 1, 
                nsnap_ancestor = 20, 
                pq_prop = {'slope': slope, 'yint': yint}, 
                tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6},
                Mrcut=18)

            fqs.append(fq)
            ssfrs.append(ssfr)


    mu_fq = np.mean(np.array(fqs))
    std_fq = np.std(np.array(fqs))
    print mu_fq, std_fq
    mu_ssfr = np.mean(np.array(ssfrs))
    std_ssfr = np.std(np.array(ssfrs))
    print mu_ssfr, std_ssfr

    def rho_ssfr_lineage_evol(
            start_nsnap = 13, 
            final_nsnap = 1, 
            tau_prop = {'name': 'instant'}, 
            Mrcut=18, 
            **kwargs
            ): 
         Compare sSFR distribution of evolved CenQue and
        SDSS Group Catalog in the green valley

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


    def ssfr_descendant(
            nsnap_descendant, 
            nsnap_ancestor = 13, 
            pq_prop = {'slope': 0.0, 'yint': 0.0},
            tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
            **kwargs
            ):
        ''' SSFR distribution of specified descendant snapshot 
        '''

        bloodline = sf_inherit(
                [nsnap_descendant], 
                nsnap_ancestor = nsnap_ancestor, 
                pq_prop = pq_prop, 
                tau_prop = tau_prop
                )

        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))

        ssfr = Ssfr()

        return ssfr.cenque(descendant)

    def fq_descendant(
            nsnap_descendant, 
            nsnap_ancestor = 13, 
            tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
            **kwargs
            ):
        ''' SSFR distribution of specified descendant snapshot 
        '''

        bloodline = sf_inherit(
                [nsnap_descendant], 
                nsnap_ancestor = nsnap_ancestor, 
                tau_prop = tau_prop
                )

        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
        
        return cq_fq(descendant)
    """
