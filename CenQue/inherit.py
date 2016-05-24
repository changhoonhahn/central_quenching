import os 
import time
import pickle 
import numpy as np
from scipy import interpolate

from lineage import Lineage

import sfr_evol

from gal_prop import Fq
from gal_prop import dFqdt
import util.util as Util
from util.util import get_nsnap_t
from util.util import intersection_index

import matplotlib.pyplot as plt 

# spline between z and t_cosmic
z_snap, t_snap = np.loadtxt(Util.snapshottable(), unpack=True, usecols=[2, 3]) 
z_of_t = interpolate.interp1d(list(reversed(t_snap)), list(reversed(z_snap)), kind='cubic') 

class Inherit(object): 
    '''
    Evolve star formation properties of 'ancestor' CentralGalaxyPopulation class at 
    redshift specified by nsnap_ancestor to descendant CentralGalaxlyPopulation object
    at redshift specified by nsnap_descendant. Both ancestor and descendant objects
    are attributes in the Lineage class, which contains all the 'lineage' information 
    i.e. all the halo tracking information. 
    '''
    def __init__(self, nsnap_descendants, nsnap_ancestor=20, 
            subhalo_prop=None, sfr_prop=None, evol_prop=None, 
            gv=None, tau=None, fudge=None, quiet=True):
        ''' 
        Parameters
        ----------
        nsnap_ancestor : int
            Snapshot number of ancestor CGPop object. The ancestor object is constructed 
            with from subhalo catalog with subhalo_prop properties. They are also assigned
            SFRs with sfr_prop properties. 
        subhalo_prop : dict
            Dictionary that describes the subhalo properties. The key 'scatter' corresponds
            to the M*-M_halo relation. The key 'soruce' describes the source SMF used for 
            the SHAM masses.
        sfr_prop : dict
            Dictionary that describes the SFR properties assigned to the ancestor CenQue object. 
            The key 'fq' describes the quiescent fraction used for the ancestor while the key
            'sfr' describes the properties of the SFR assignment. 
        evol_prop : dict
            Dictionary that consists of dictionaries which each describe paramter choices in 
            the model. 
            - evol_prop['pq'] dictates the quenching properties. 
            - evol_prop['tau'] dictates the quenching timescale. 
            - evol_prop['sfr'] dictates the SFR evolution. 
            - evol_prop['mass'] dictates the mass evolution. 
        '''
        self.quiet = quiet 
        
        self.nsnap_ancestor = nsnap_ancestor     # snapshot of ancestor
        if not isinstance(nsnap_descendants, list): 
            self.nsnap_descendants = [nsnap_descendants]
        else:  
            self.nsnap_descendants = nsnap_descendants
        self.tsnap_descendants = [Util.get_t_nsnap(nsnap) for nsnap in self.nsnap_descendants] 
        self.subhalo_prop = subhalo_prop 
        self.sfr_prop = sfr_prop 
        self.evol_prop = evol_prop
        self.fudge = fudge
        self.tau = tau
        self.gv = gv
        self._InitProps()
        self._UnpackProps() # unpack dictionaries 
        
        # Read in Lineage object and import necessary GalPop objects
        self._Read_Lineage()
        
        # Evolve quiescent population 
        self._QuiescentEvol() 

        # Evolve star-forming population 
        self._StarformingEvol()
        
        # Account for overquenching 
        self._Overquenching()

    def __call__(self):  
        '''
        '''
        self.descendant_dict['ancestor'] = self.ancestor
        for key in self.descendant_dict.keys(): 
            cgpop_obj = self.descendant_dict[key]
            #self.descendant_dict[key] = cgpop_obj._Jettison() 

        return self.descendant_dict

    def _StarformingEvol(self): 
        ''' Evolve the star-forming galaxy population. This involves evolving their 
        stellar mass and star formation rates while quenching some of them to match 
        the quiescent fraction, fQ(M*, z), simultaneously. 
        '''
        for nsnap_d in self.descendant_dict.keys(): 
            allwill = self.descendant_dict[str(nsnap_d)].will.copy()
            allsucc = self.descendant_dict[str(nsnap_d)].succession.copy()
            anc_sfr_class = self.ancestor.sfr_class.copy()
            
            q_ancestors = np.where(anc_sfr_class[allwill] == 'quiescent')[0]  # Q ancestors
            sf_ancestors = np.where(anc_sfr_class[allwill] == 'star-forming')[0]  # SF ancestors

            # star-formation duty cycle parameters (takes a long time ~15 sec)
            if not self.quiet: 
                dutycycle_time = time.time()
            self.dutycycle_prop = sfr_evol.dutycycle_param(
                    len(sf_ancestors), 
                    dutycycle_prop=self.evol_prop['sfr']['dutycycle'])
            self.dutycycle_prop['delta_sfr'] = self.ancestor.delta_sfr[allwill[sf_ancestors]].copy()
            if not self.quiet: 
                print 'SFR dutycycle properties take ', time.time() - dutycycle_time, ' to generate'
        
            logsfr, t_quench = self._Evol_MshamSFR(nsnap_d, allwill, sf_ancestors, q_ancestors)

            self.descendant_dict[str(nsnap_d)].sfr[allsucc[sf_ancestors]] = logsfr
            self.descendant_dict[str(nsnap_d)].ssfr[allsucc[sf_ancestors]] = \
                    self.descendant_dict[str(nsnap_d)].sfr[allsucc[sf_ancestors]] - \
                    self.descendant_dict[str(nsnap_d)].mass[allsucc[sf_ancestors]]

            is_qing = np.where(t_quench < self.descendant_dict[str(nsnap_d)].t_cosmic) 
            is_notqing = np.where(t_quench >= self.descendant_dict[str(nsnap_d)].t_cosmic) 

            self.descendant_dict[str(nsnap_d)].sfr_class[is_qing] = 'quiescent'
            self.descendant_dict[str(nsnap_d)].sfr_class[is_notqing] = 'star-forming'

        return None

    def _Evol_MshamSFR(self, nsnap_d, allwill, sf_ancestor, q_ancestor): 
        ''' Evolve SFR of the SF galaxies while quenching some of them at the same time. 

        Notes
        -----
        * SF galaxies are quenched based on the following prescription:
            The number of galaxies that start quenching between t0 and t0+tstep 
            is determined by first guessing, 

            N_quenching = N_sf(t0) * (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep
        '''
        tQ0 = self.ancestor.tQ[allwill[sf_ancestor]].copy()
        tQ = self.ancestor.tQ[allwill[sf_ancestor]].copy()
        Mq = self.ancestor.MQ[allwill[sf_ancestor]].copy()
        
        # earliest starting time (cosmic time of ancestor snapshot) 
        t0 = self.ancestor.tsnap_genesis[allwill[sf_ancestor]].copy() 
        t00 = t0.min()              

        i_low = np.abs(t_snap-np.max(self.tsnap_descendants)).argmin()-1
        i_high = np.abs(t_snap-t00).argmin()+1
        t_evol = t_snap[i_low:i_high][::-1]

        qf = Fq()
        Fq_anal = qf.model
        M_bins = np.arange(5.5, 13., 0.5) 
        M_mid = 0.5 * (M_bins[:-1] + M_bins[1:]) 

        q_Mshams = self.MshamEvol[allwill[q_ancestor],:].copy() 
        sf_Mshams = self.MshamEvol[allwill[sf_ancestor],:].copy() 

        for i_t, tt in enumerate(t_evol[:-1]): 
            t_step = t_evol[i_t+1] - tt 
            if not self.quiet: 
                print 't_cosmic = ', tt
                t_one_tstep = time.time()
            # M_sham of the quiescent galaxies at the closest snapshot. 
            closest_t0_snap = np.abs(t_snap-tt).argmin() - 1
            Msham_q0 = q_Mshams[:, closest_t0_snap]
            Msham_sf0 = sf_Mshams[:, closest_t0_snap]

            within = np.where((Msham_sf0 > 0.) & (t0 <= tt))    # tsnap_genesis <= t
            sf_within = np.where((Msham_sf0 > 0.) & (t0 <= tt) & (tQ == 999.))[0] # Not quenching SF galaxies 
            Nsf_0 = len(sf_within)          
            M_t0_sample = np.concatenate(
                    [Msham_sf0[within], Msham_q0[np.where(Msham_q0 > 0.)]]
                    )
            P_sf = np.random.uniform(0., 1., Nsf_0)
            
            # Initial estimate of quenching probability given by dFq/dt 
            # P_q = (1/( 1 - fQ(t0) )) * dFq/dt(t0 + 0.5 tstep) * tstep
            Ng0, dum = np.histogram(M_t0_sample, bins=M_bins)
            Nsf0, dum = np.histogram(Msham_sf0[sf_within], bins=M_bins)
            P_q_arr = Ng0/Nsf0 * \
                    dFqdt(M_mid, tt + 0.5 * t_step, lit=self.fq_prop['name']) * t_step
            P_q_arr[np.where(Nsf0 == 0)] = 0.
            if not self.quiet: 
                print 'M* : ', M_mid[-6:]
                print 'P_q : ', P_q_arr[-6:]
            Pq_M = interpolate.interp1d(M_mid, P_q_arr, kind='linear') 
            P_q = Pq_M(Msham_sf0[sf_within]) 
            q_ing0 = np.where(P_q > P_sf)
            Nqing0 = len(q_ing0[0])

            # M_sham of the quiescent galaxies at the closest snapshot. 
            closest_t1_snap = np.abs(t_snap - t_evol[i_t+1]).argmin() - 1 
            Msham_q1 = q_Mshams[:, closest_t1_snap]
            Msham_sf1 = sf_Mshams[:, closest_t1_snap]
            M_t1_sample = np.concatenate([Msham_sf1[within], Msham_q1])
            
            # dPQ correction to the quenching to account for change in Ng(M*,t)
            Ngp, dum = np.histogram(M_t1_sample, bins=M_bins)
            dPq = (Ngp - Ng0)/Nsf0.astype('float') * \
                    Fq_anal(M_mid, z_of_t(tt + t_step), lit=self.fq_prop['name'])
            dPq[np.where(dPq < 0.)] = 0.
            dPq[np.where(Nsf0 == 0)] = 0.
            if not self.quiet: 
                print 'dPq : ', dPq[-6:]
            dPq_M = interpolate.interp1d(M_mid, dPq, kind='linear') 
            P_q += dPq_M(Msham_sf0[sf_within])
            # fudge factor for P_Q
            fudge_factor = self.fudge_prop['slope'] * (Msham_sf0[sf_within] - self.fudge_prop['fidmass']) + self.fudge_prop['offset']
            fudge_factor[np.where(fudge_factor < 1.)] = 1.
            P_q *= fudge_factor 
            
            q_ing = np.where(P_sf < P_q)
            if not self.quiet: 
                print 'Initial guess ', Nqing0, ' final: ', len(q_ing[0]), ' SF gal out of ', Nsf_0, '  start quenching'
                print time.time() - t_one_tstep
                print '----------------'
            # assign them quenching times then actually evolve their stellar masses  
            tQ[sf_within[q_ing]] = np.random.uniform(low=tt, high=tt+t_step, size=len(q_ing[0]))

        quenching = np.where((Mq == -999.) & (tQ < 999.) & (tQ > 0.))[0] 
        Nquenching = len(quenching)
        closest_tQ_index = np.abs(
                np.tile(t_snap, (Nquenching, 1)) - np.tile(tQ[quenching].reshape(Nquenching, 1), (1, len(t_snap)))
                ).argmin(axis=1) - 1
        Mq[quenching] = sf_Mshams[quenching, closest_tQ_index]

        z0 = z_of_t(t0)     # initial redshift
        gv0 = np.where(tQ0 == 0.)[0]    # galaxies that were initial quenching 
        
        t_f = self.descendant_dict[str(nsnap_d)].t_cosmic
        z_f = self.descendant_dict[str(nsnap_d)].zsnap

        closest_tf_snap = np.abs(t_snap - t_f).argmin() - 1
        Msham_f = sf_Mshams[:, closest_tf_snap].copy()

        q_ed = np.where(t_f > tQ)
        tmp_M = Msham_f.copy()
        tmp_M[q_ed] = Mq[q_ed].copy()

        # log(SFR)_SFMS evolution from t0
        logsfr_sfms = sfr_evol.AverageLogSFR_sfms(tmp_M, z0, sfms_prop=self.sfms_prop) + \
                sfr_evol.DeltaLogSFR_sfms(z0, z_f, sfms_prop=self.sfms_prop) 

        # log(SFR)_duty cycle evolution from t0 to tQ
        logsfr_sfduty = sfr_evol.DeltaLogSFR_dutycycle(
                t0, 
                t_f, 
                t_q=tQ, 
                dutycycle_prop=self.dutycycle_prop
                )

        logsfr_quench = sfr_evol.DeltaLogSFR_quenching(
                tQ, 
                t_f, 
                M_q=Mq,
                tau_prop=self.tau_prop)

        logSFR = logsfr_sfms + logsfr_quench + logsfr_sfduty 
        
        # correct for the galaxies that were originally in the green valley. 
        # ultimately doesn't make a difference in the outcome, but to be meticulous 
        logSFR[gv0] = (self.ancestor.sfr[allwill[sf_ancestor]])[gv0] + \
                sfr_evol.DeltaLogSFR_sfms(z0[gv0], z_f, sfms_prop=self.sfms_prop) + \
                sfr_evol.DeltaLogSFR_quenching(
                        np.repeat(self.ancestor.t_cosmic, len(gv0)), 
                        t_f, M_q=Mq[gv0], tau_prop=self.tau_prop)

        return logSFR, tQ
    
    def _Overquenching(self): 
        ''' Account for overquenching 

        Assign final quenched SSFR to star-forming galaxies. This is specified 
        due to the fact that there's a lower bound on the SSFR. These values are effectively 
        hardcoded in order to reproduce the quiescent peak of the SSFR 
        distribution, which is more of a lower bound. 
        '''
        for n_d in self.descendant_dict.keys(): 

            des = self.descendant_dict[str(n_d)]
            sf_succ = des.succession[des.sf_ancestor]

            avg_q_ssfr = sfr_evol.AverageLogSSFR_q_peak(des.mass[sf_succ])
            sigma_q_ssfr = sfr_evol.ScatterLogSSFR_q_peak(des.mass[sf_succ])

            min_q_ssfr = sigma_q_ssfr * np.random.randn(len(sf_succ)) + avg_q_ssfr 

            self.descendant_dict[str(n_d)].min_ssfr[sf_succ] = min_q_ssfr.copy()

            # Deal with over quenched galaxies 
            overquenched = np.where(
                    self.descendant_dict[str(n_d)].min_ssfr[sf_succ] > self.descendant_dict[str(n_d)].ssfr[sf_succ]
                    )
            if len(overquenched[0]) > 0: 
                self.descendant_dict[str(n_d)].ssfr[sf_succ[overquenched]] = \
                        self.descendant_dict[str(n_d)].min_ssfr[sf_succ[overquenched]] 
                self.descendant_dict[str(n_d)].sfr[sf_succ[overquenched]] = \
                    self.descendant_dict[str(n_d)].ssfr[sf_succ[overquenched]] + \
                    self.descendant_dict[str(n_d)].mass[sf_succ[overquenched]]
            del n_d
            del des 
        return None 

    def _QuiescentEvol(self): 
        ''' Evolve quiescent ancestor galaxy population to descendant. 
        Both ancestor and descendant have to be CGPop objects. Galaxies 
        that are quiescent at z_init remain quiescent until z_final. 
        There are no mechanisms to induce star forming. These galaxies are 
        evolved by keeping the SSFR constant. M*(zf) is the SHAM M* at
        redshift zf. 
        '''     
        for n_d in self.descendant_dict.keys():
            des = self.descendant_dict[str(n_d)]

            # quiescent ancestors
            a_index = des.will[des.q_ancestor]
            
            # matching descendants 
            d_index = des.succession[des.q_ancestor]

            self.descendant_dict[str(n_d)].sfr_class[d_index] = 'quiescent'
            self.descendant_dict[str(n_d)].ssfr[d_index] = self.ancestor.ssfr[a_index].copy()
            self.descendant_dict[str(n_d)].sfr[d_index] = \
                    self.descendant_dict[str(n_d)].ssfr[d_index] +\
                    self.descendant_dict[str(n_d)].mass[d_index]
            del des
            del n_d
        return None  

    def logSFR_M_t(self, logmass, t_input, t_init=None, t_q=None, M_q=None): 
        ''' log(SFR) as a function of M* and t_cosmic.

        Notes
        -----
        * kwarg indices is for dutycycle_prop 
        '''
        # SFR evolution based on solving an ODE of SFR
        logsfr_time = time.time()
        z_init = z_of_t(t_init) # initial redshift
        
        # update quenched M* ( this is the stellar mass roughly at the time of quenching)
        #just_quenched = np.where(
        #        (t_input > t_q) & 
        #        (M_q == -999.))
        #if len(just_quenched[0]) > 0: 
        #    M_q[just_quenched] = logmass[just_quenched].copy()

        # average SFR of SFMS at M* and z_init
        quenched = np.where(t_input > t_q)
        tmp_M = logmass.copy()
        tmp_M[quenched] = M_q[quenched].copy()
        avglogsfr = sfr_evol.AverageLogSFR_sfms(tmp_M, z_init, sfms_prop=self.sfms_prop)
        del tmp_M 

        # log(SFR)_SFMS evolutionfrom t0
        logsfr_sfms = sfr_evol.DeltaLogSFR_sfms(
                z_init, 
                z_of_t(t_input),
                sfms_prop=self.sfms_prop)

        # log(SFR)_duty cycle evolution from t0 to tQ
        logsfr_sfduty = sfr_evol.DeltaLogSFR_dutycycle(
                t_init, 
                t_input, 
                t_q=t_q, 
                dutycycle_prop=self.dutycycle_prop, 
                indices=None)

        logsfr_quench = sfr_evol.DeltaLogSFR_quenching(
                t_q, 
                t_input, 
                M_q=M_q,
                tau_prop=self.tau_prop)

        print logsfr_quench

        logsfr_tot = avglogsfr + logsfr_sfms + logsfr_sfduty + logsfr_quench
         
        return logsfr_tot

    def _Read_Lineage(self): 
        ''' Read in Lineage object and import all the GalPop objects within 
        lineage to Inherit object. 
        '''
        # read in the lineage (< 0.05 seconds for one snapshot)
        if not self.quiet: 
            read_time = time.time()
        bloodline = Lineage(
                nsnap_ancestor=self.nsnap_ancestor, 
                subhalo_prop=self.subhalo_prop, quiet=self.quiet)
        bloodline.Read(self.nsnap_descendants, quiet=self.quiet)

        #DEFUNCT
        #if 'subhalogrowth' in self.sfr_prop.keys():  
        #    # depending on whether SFR assign includes subhalo growth AM
        #    sfr_prop['subhalogrowth']['nsnap_descendant'] = nsnap_descendant

        # assign SFR to ancestor object
        bloodline.AssignSFR_ancestor(sfr_prop=self.sfr_prop, quiet=self.quiet)
        if not self.quiet: 
            print 'Lineage Read Time = ', time.time() - read_time 

        self.ancestor = bloodline.ancestor    # ancestor object
        self.ancestor.tQ[np.where(self.ancestor.tQ != 999.)] = 0.
        self.MshamEvol = self.ancestor.Msham_evol[0]

        self._ImportDescendants(bloodline)

        return None 
        
    def _ImportDescendants(self, bloodline): 
        ''' Import specified descendants from bloodline into 
        a descendant dictionary 
        '''
        anc_sfr_class = self.ancestor.sfr_class.copy()
        anc_snap_index = self.ancestor.snap_index

        self.descendant_dict = {} 
        for n_d in self.nsnap_descendants: 
            des = getattr(bloodline, 'descendant_snapshot'+str(n_d))
            des._clean_initialize()   # initialize SF properties
            des.sfr_prop = self.ancestor.sfr_prop 
            # match indices up with each other 
            succession, will = intersection_index(
                    getattr(des, 'ancestor'+str(self.nsnap_ancestor)), 
                    anc_snap_index
                    )
            if len(succession) != len(des.mass): 
                raise ValueError('Something wrong with the lineage')
            des.succession = succession
            des.will = will
    
            q_ancestor = np.where(anc_sfr_class[will] == 'quiescent')[0]      # Q ancestors
            sf_ancestor = np.where(anc_sfr_class[will] == 'star-forming')[0]  # SF ancestors
            if not self.quiet: 
                print "nsnap_descendant = ", n_d
                print "Ancestors: Nq = ", len(q_ancestor), ', Nsf = ', len(sf_ancestor)

            des.q_ancestor = q_ancestor
            des.sf_ancestor = sf_ancestor

            self.descendant_dict[str(n_d)] = des

            del des
            del n_d
        return None

    def _InitProps(self): 
        ''' Initial the dictionaries that specify properties
        '''
        if self.subhalo_prop is None:   # subhalo
            self.subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}
        if self.sfr_prop is None:   # sfr
            fq_prop = {'name': 'wetzel'}     # quiescent fraction
            if self.gv is None:         # green valley
                gv_prop = {'fidmass': 10.5, 'slope': 0.0, 'offset': 0.0}  
            elif isinstance(self.gv, list): 
                gv_prop = {'fidmass': 10.5, 'slope': self.gv[0], 'offset': self.gv[1]}  
            elif isinstance(self.gv, dict): 
                gv_prop = self.gv
            sfms_prop = {'name': 'linear', 'zslope': 1.1}   # SFMS 
            self.sfr_prop = {'fq': fq_prop, 'sfms': sfms_prop, 'gv': gv_prop}

        if self.evol_prop is None: 
            mass_prop = {'name': 'sham'} 
            dutycycle_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3}

            if self.tau == 'satellite': 
                tau_dict = {'name': 'satellite'} 
            elif isinstance(self.tau, list): 
                tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': self.tau[0], 'yint': self.tau[1]}
            elif isinstance(self.tau, dict): 
                tau_dict = self.tau 
            del self.tau 

            if self.fudge is None:
                fudge_dict = {'slope': 0.0, 'fidmass': 10.5, 'offset': 1.0}
            elif isinstance(self.fudge, list): 
                fudge_dict = {'fidmass': 10.5, 'slope': self.fudge[0], 'offset': self.fudge[1]}
            elif isinstance(self.fudge, dict): 
                fudge_dict = fudge 
            del self.fudge
            self.evol_prop = {
                    'sfr': {'dutycycle': dutycycle_prop},
                    'mass': mass_prop,
                    'tau': tau_dict, 
                    'fudge': fudge_dict
                    }
        return None 

    def _UnpackProps(self): 
        ''' Unpack dictionaries that specify properties  
        '''
        self.sfms_prop = self.sfr_prop['sfms']        # SFMS properties
        self.fq_prop = self.sfr_prop['fq']            # fq_prop 
        self.tau_prop = self.evol_prop['tau']         # Quenching timescale properties
        self.fudge_prop = self.evol_prop['fudge']
        self.dutycycle_prop = self.evol_prop['sfr']['dutycycle']   # SF dutycycle properties 
        return None

    def _Jettison(self): 
        ''' Inherit class is usually very massive. Jettison needless data 
        '''
        return None 


class ABCInherit(object):  
    def __init__(self, t, abcrun=None, prior_name=None): 
        ''' A wrapper class for Inherit results with ABC best-fit parameters
        '''
        if abcrun is None: 
            raise ValueError('specify abcrun') 
        if prior_name is None: 
            raise ValueError('specify prior_name') 

        self.t = t
        self.abcrun= abcrun
        self.prior_name = prior_name 
    
        # Read in the theta values of the ABC run with 'prior_name' priors 
        # at time step t
        theta_file = ''.join([Util.code_dir(), 
            'dat/pmc_abc/', 
            'CenQue_theta_t', str(self.t), '_', self.abcrun, '.dat']) 
        w_file = ''.join([Util.code_dir(), 
            'dat/pmc_abc/', 
            'CenQue_w_t', str(t), '_', abcrun, '.dat']) 
        self.theta = np.loadtxt(theta_file)      # theta values 
        self.w = np.loadtxt(w_file)              # w values 

        self.med_theta = [np.median(self.theta[:,i]) for i in range(len(self.theta[0]))]

        self.sim_kwargs = self._ABCrun_kwargs() 
    
    def _ABCrun_kwargs(self): 
        '''
        '''
        gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = self.med_theta

        sfinherit_kwargs = self._ReadABCrun()
        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {
                'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset
                }
        sim_kwargs['evol_prop']['fudge'] = {
                'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset
                }
        sim_kwargs['evol_prop']['tau'] = {
                'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset
                }
        
        return sim_kwargs

    def PickleFile(self, snapshots): 
        ''' Name of the pickle file where the inherited snapshots are saved  
        '''
        if isinstance(snapshots, float): 
            snapshots = [snapshots] 
        elif isinstance(snapshots, list): 
            pass
        elif isinstance(snapshots, np.ndarray): 
            pass 
        elif isinstance(snapshots, str): 
            if snapshots == 'all': 
                snapshots = range(1, self.sim_kwargs['nsnap_ancestor']) 
            else: 
                raise ValueError
        else: 
            raise ValueError
        self.snapshots = snapshots

        if snapshots == range(1, self.sim_kwargs['nsnap_ancestor']): 
            snap_str = '_all' 
        else: 
            snap_str = '_'.join([str(snap) for snap in self.snapshots])
        
        abcinh_file = ''.join([Util.code_dir(), 'dat/pmc_abc/pickle/', 
            'Inherit', 
            '.t', str(self.t), 
            '.snapshots', snap_str, 
            '.ancestor', str(self.sim_kwargs['nsnap_ancestor']), 
            '.abcrun_', self.abcrun, '.p']) 
        return abcinh_file  

    def PickleDump(self, snapshots, clobber=False):  
        ''' Run Inherit for the snapshots and save to pickle file 
        '''
        pickle_file = self.PickleFile(snapshots) 

        if os.path.isfile(pickle_file) and not clobber: 
            raise ValueError("Already exists. Specify clobber=True, if you want to overwrite") 

        inh = Inherit(self.snapshots, 
                nsnap_ancestor=self.sim_kwargs['nsnap_ancestor'],
                subhalo_prop=self.sim_kwargs['subhalo_prop'], 
                sfr_prop=self.sim_kwargs['sfr_prop'], 
                evol_prop=self.sim_kwargs['evol_prop'])

        inh_dict = inh() 
        pickle.dump(inh_dict, open(pickle_file, 'wb'))
        return inh_dict

    def PickleLoad(self, snapshots, clobber=False): 
        ''' Load the precomputed Inherit saved to pickel file 
        '''
        pickle_file = self.PickleFile(snapshots) 
        if not os.path.isfile(pickle_file) or clobber: 
            inh = Inherit(self.snapshots, 
                    nsnap_ancestor=self.sim_kwargs['nsnap_ancestor'],
                    subhalo_prop=self.sim_kwargs['subhalo_prop'], 
                    sfr_prop=self.sim_kwargs['sfr_prop'], 
                    evol_prop=self.sim_kwargs['evol_prop'])

            inh_dict = inh() 
            pickle.dump(inh_dict, open(pickle_file, 'wb'))
        else: 
            inh_dict = pickle.load(open(pickle_file, 'rb'))

        return inh_dict

    def _ReadABCrun(self): 
        ''' Read file that specifies all the choices of parameters and ABC settings.
        '''
        abcrun_file = ''.join([Util.code_dir(), 
            'dat/pmc_abc/run/', 'abcrun_', self.abcrun, '.p']) 
        return pickle.load(open(abcrun_file, 'rb'))
