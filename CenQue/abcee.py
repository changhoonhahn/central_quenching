'''

Module for ABC-PMC inference for the Central Quenching project

'''
import time
import pickle
import numpy as np
from datetime import datetime 

# --- Local ---
from util.util import code_dir
from gal_prop import Fq 
from observations import GroupCat
from sf_inherit import InheritSF

import corner
import abcpmc
from abcpmc import mpi_util

from plotting.plots import QAplot

import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors

def DataSummary(Mrcut=18, observables=['ssfr']): 
    ''' Summary statistics of the data. In our case that is the 
    SSFR distribution of the SDSS group catalog.
    '''
    obvs = []
    if 'ssfr' in observables:
        # Group Catalog object
        groupcat = GroupCat(Mrcut=Mrcut, position='central')
        # SSFR distribution of group catalog
        bins, dist = groupcat.Ssfr()   
        obvs.append([np.array(bins), np.array(dist)])
    if 'fqz03' in observables: 
        qfrac = Fq()
        M_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
        M_mid = 0.5 * (M_bin[:-1] + M_bin[1:])
        fq_model = qfrac.model(M_mid, 0.3412, lit='wetzel')
        obvs.append([M_mid, fq_model])
    
    if len(observables) == 1: 
        obvs = obvs[0]
    return obvs

def SimSummary(observables=['ssfr'], **sim_kwargs):
    ''' Summary statistics of the simulation. In our case, the simulation is InheritSF 
    and the summary statistic is the SSFR distribution. 
    '''
    obvs = []
    if 'ssfr' in observables:
        bloodline = InheritSF(
                1,
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        descendant = getattr(bloodline, 'descendant_snapshot1') 
        
        bins, dist = descendant.Ssfr()
        obvs.append([np.array(bins), np.array(dist)])
    if 'fqz03' in observables: 
        bloodline = InheritSF(
                6,
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        descendant = getattr(bloodline, 'descendant_snapshot6') 

        qfrac = Fq()
        M_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
        M_mid = 0.5 * (M_bin[:-1] + M_bin[1:]) 

        sfq = qfrac.Classify(descendant.mass, descendant.sfr, descendant.zsnap, 
                sfms_prop=sim_kwargs['sfr_prop']['sfms'])

        ngal, dum = np.histogram(descendant.mass, bins=M_bin)
        ngal_q, dum = np.histogram(descendant.mass[sfq == 'quiescent'], bins=M_bin)

        obvs.append([M_mid, ngal_q.astype('float')/ngal.astype('float')])
    
    if len(observables) == 1: 
        obvs = obvs[0]
    return obvs
    
def PriorRange(name): 
    ''' Return predetermined prior ranges based on their assigned name.

    theta = [gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset]
    '''
    if name =='try0':   # guess 
        prior_min = [0.0, 0.0, -3., 1., -1.5, 0.01]
        prior_max = [1.0, 1.0, 0.0, 3., 0.0, 1.0]
    else: 
        raise NotImplementedError

    return [prior_min, prior_max] 

def RhoSSFR(simum, datum):
    ''' L2 Norm between the data SSFR distribution and the simulation SSFR distribution. 
    
    Parameters
    ----------
    datum : list
        List [bins, dist] where bins is a np.ndarray of the SSFR bin values and 
        dist is a np.ndarray that specifies the SSFR distribution of the group catalog.
    model : list 
        List [bins, dist] where bins is a np.ndarray of the SSFR bin values and 
        dist is a np.ndarray that specifies the SSFR distribution of the simulation.

    '''
    datum_dist = datum[1]
    simum_dist = simum[1]
    drho = np.sum((datum_dist - simum_dist)**2)
    return drho

def RhoSSFR_Fq(simum, datum): 
    ''' Multivariate distance metric. First dimension is 
    L2 Norm between the data SSFR distribution and the simulation SSFR distribution. 
    While the second dimension is L2 Norm between the model fQ(z=0.34) and the 
    simulated fQ at snapshot 6. 
    
    Parameters
    ----------
    datum : list
        List [bins, dist] where bins is a np.ndarray of the SSFR bin values and 
        dist is a np.ndarray that specifies the SSFR distribution of the group catalog.
    model : list 
        List [bins, dist] where bins is a np.ndarray of the SSFR bin values and 
        dist is a np.ndarray that specifies the SSFR distribution of the simulation.

    '''
    
    datum_ssfr = datum[0][1]
    simum_ssfr = simum[0][1]
    rho_ssfr = np.sum((datum_ssfr - simum_ssfr)**2)
    
    datum_fq = datum[1][1]
    simum_fq = simum[1][1]
    rho_fq = np.sum((datum_fq - simum_fq)**2)
    return [rho_ssfr, rho_fq]

def MakeABCrun(abcrun=None, nsnap_start=15, subhalo=None, fq=None, sfms=None, dutycycle=None, mass_evol=None, 
        Niter=None, Npart=None, prior_name=None, eps_val=None, restart=False): 
    '''Pickle file that specifies all the choices of parameters and ABC settings 
    for record keeping purposes.
    '''
    if abcrun is None: 
        abcrun = ''.join(['RENAME_ME_', str(datetime.today().date())])
    restart_str = ''
    if restart: 
        restart_str = '.restart'
    file_name = ''.join([code_dir(), 'dat/pmc_abc/run/', 'abcrun_', abcrun, restart_str, '.txt']) 
    pickle_name = file_name.replace('.txt', '.p')

    f = open(file_name, 'w')

    # ABC properities
    f.write('# ABC Run Properties \n') 
    f.write('N_iter = '+str(Niter)+'\n')
    f.write('N_particles = '+str(Npart)+'\n')
    f.write('Prior Name = '+prior_name+'\n')
    f.write('Epsilon_0 = '+str(eps_val)+'\n')
    if restart: 
        f.write('Restarting')
    f.write('\n')

    # Subhalo properties
    f.write('# Subhalo Properties \n') 
    if subhalo is None:    
        subhalo = {'scatter': 0.2, 'source': 'li-march'}
    f.write(''.join(['scatter = ', str(subhalo['scatter']), '\n']))
    f.write(''.join(['source = ', subhalo['source'], '\n']))
    f.write('\n')
    
    # Quiescent fraction properties
    f.write('# Quiescent Fraction Properties \n') 
    if fq is None:          # quiescent fraction
        fq = {'name': 'wetzel'}
    f.write(''.join(['name = ', fq['name'], '\n']))
    f.write('\n')

    # Starforming Main Sequence
    f.write('# Starforming Main Sequence Properties \n') 
    if sfms is None:   # SFMS  
        sfms = {'name': 'linear', 'zslope': 1.14}
    f.write(''.join(['name = ', sfms['name'], '\n']))
    f.write(''.join(['zslope = ', str(sfms['zslope']), '\n']))
    f.write('\n')
    
    # SF dutycycle properties
    f.write('# SF Dutycycle Properties \n') 
    if dutycycle is None:   
        dutycycle = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3}
    f.write(''.join(['name = ', dutycycle['name'], '\n']))
    f.write(''.join(['frequency range = ', str(dutycycle['freq_range'][0]), ',', str(dutycycle['freq_range'][1]), '\n']))
    f.write(''.join(['sigma = ', str(dutycycle['sigma']), '\n']))
    f.write('\n')
    
    # Mass Evolution properties
    f.write('# Mass Evolution Properties \n') 
    if mass_evol is None:
        mass_evol = {'name': 'sham'} 
    if mass_evol['name'] == 'sham': 
        f.write(''.join(['name = ', mass_evol['name'], '\n']))
    f.close()

    evol_dict = { 
            'type': 'simult',
            'sfr': {'dutycycle': dutycycle},
            'mass': mass_evol
            }
    kwargs = {
            'nsnap_ancestor': nsnap_start, 
            'subhalo_prop': subhalo, 
            'sfr_prop': {'fq': fq, 'sfms': sfms},
            'evol_prop': evol_dict
            }

    pickle.dump(kwargs, open(pickle_name, 'wb'))
    return kwargs, abcrun 

def ReadABCrun(abcrun, restart=False): 
    ''' Read text file that specifies all the choices of parameters and ABC settings.
    '''
    if abcrun is None: 
        abcrun = ''.join(['RENAME_ME_', str(datetime.today().date())])
    restart_str = ''
    if restart: 
        restart_str = '.restart'
    pickle_name = ''.join([code_dir(), 'dat/pmc_abc/run/', 'abcrun_', abcrun, restart_str, '.p']) 
    return pickle.load(open(pickle_name, 'rb')), abcrun


def ABC(T, eps_input, Npart=1000, prior_name='try0', observables=['ssfr'], abcrun=None, 
        restart=False, t_restart=None, eps_restart=None):
    ''' ABC-PMC implementation. 

    Parameters
    ----------
    T : (int) 
        Number of iterations

    eps_input : (float)
        Starting epsilon threshold value 

    N_part : (int)
        Number of particles

    prior_name : (string)
        String that specifies what prior to use.

    abcrun : (string)
        String that specifies abc run information 
    '''
    if isinstance(eps_input, list):
        if len(eps_input) != len(observables): 
            raise ValueError
    if len(observables) > 1 and isinstance(eps_input, float): 
        raise ValueError

    # output abc run details
    sfinherit_kwargs, abcrun_flag = MakeABCrun(
            abcrun=abcrun, Niter=T, Npart=Npart, prior_name=prior_name, 
            eps_val=eps_input, restart=restart) 

    # Data 
    data_sum = DataSummary(observables=observables)
    # Priors
    prior_min, prior_max = PriorRange(prior_name)
    prior = abcpmc.TophatPrior(prior_min, prior_max)    # ABCPMC prior object

    ## Read in ABC run file
    #sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun, restart=restart)

    def Simz(tt):       # Simulator (forward model) 
        gv_slope = tt[0]
        gv_offset = tt[1]
        fudge_slope = tt[2]
        fudge_offset = tt[3]
        tau_slope = tt[4]
        tau_offset = tt[5]

        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
        sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
        sim_kwargs['evol_prop']['tau'] = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

        return SimSummary(observables=observables, **sim_kwargs)

    theta_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'CenQue_theta_t', str(pewl), '_', abcrun_flag, '.dat']) 
    w_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'CenQue_w_t', str(pewl), '_', abcrun_flag, '.dat']) 
    dist_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'CenQue_dist_t', str(pewl), '_', abcrun_flag, '.dat']) 
    eps_file = ''.join([code_dir(), 
        'dat/pmc_abc/epsilon_', abcrun_flag, '.dat'])

    if observables == ['ssfr']:
        distfn = RhoSSFR
    elif observables == ['ssfr', 'fqz03']: 
        distfn = RhoSSFR_Fq
   
    if restart:
        if t_restart is None: 
            raise ValueError

        last_thetas = np.loadtxt(theta_file(t_restart))
        last_ws = np.loadtxt(w_file(t_restart))
        last_dist = np.loadtxt(dist_file(t_restart))

        init_pool = abcpmc.PoolSpec(t_restart, None, None, last_thetas, last_dist, last_ws)
    else: 
        init_pool = None 

    eps = abcpmc.ConstEps(T, eps_input)
    try:
        mpi_pool = mpi_util.MpiPool()
        abcpmc_sampler = abcpmc.Sampler(
                N=Npart,                # N_particles
                Y=data_sum,             # data
                postfn=Simz,            # simulator 
                dist=distfn,           # distance function  
                pool=mpi_pool)  
    except AttributeError: 
        abcpmc_sampler = abcpmc.Sampler(
                N=Npart,                # N_particles
                Y=data_sum,             # data
                postfn=Simz,            # simulator 
                dist=distfn)           # distance function  
    abcpmc_sampler.particle_proposal_cls = abcpmc.ParticleProposal

    pools = []
    if init_pool is None: 
        f = open(eps_file, "w")
        f.close()
    eps_str = ''
    for pool in abcpmc_sampler.sample(prior, eps, pool=init_pool):
        print '----------------------------------------'
        print 'eps ', pool.eps
        new_eps_str = '\t'+str(pool.eps)+'\n'
        if eps_str != new_eps_str:  # if eps is different, open fiel and append 
            f = open(eps_file, "a")
            eps_str = new_eps_str
            f.write(eps_str)
            f.close()

        print("T:{0},ratio: {1:>.4f}".format(pool.t, pool.ratio))
        print eps(pool.t)

        # write theta, weights, and distances to file 
        np.savetxt(theta_file(pool.t), pool.thetas, 
            header='gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset')
        np.savetxt(w_file(pool.t), pool.ws)
        np.savetxt(dist_file(pool.t), pool.dists)
    
        # update epsilon based on median thresholding 
        if len(observables) == 1: 
            eps.eps = np.median(pool.dists)
        else:
            #print pool.dists
            print np.median(np.atleast_2d(pool.dists), axis = 0)
            eps.eps = np.median(np.atleast_2d(pool.dists), axis = 0)
        print '----------------------------------------'
        pools.append(pool)

    return pools 


# ABC Corner Plot 
class PlotABC(object): 
    def __init__(self, t, abcrun=None, prior_name='try0'): 
        '''
        '''
        if abcrun is None: 
            raise ValueError('specify ABC run string') 

        self.t = t 
        self.abcrun = abcrun

        theta_file = ''.join([
            code_dir(), 'dat/pmc_abc/', 'CenQue_theta_t', str(t), '_', abcrun, '.dat']) 
        w_file = ''.join([
            code_dir(), 'dat/pmc_abc/', 'CenQue_w_t', str(t), '_', abcrun, '.dat']) 
        self.theta = np.loadtxt(theta_file)      # theta values 
        self.w = np.loadtxt(w_file)              # w values 
        
        self.med_theta = [np.median(self.theta[:,i]) for i in range(len(self.theta[0]))]
        
        # prior range 
        prior_min, prior_max = PriorRange(prior_name)
        self.prior_range = [(prior_min[i], prior_max[i]) for i in range(len(prior_min))]

        self.descendant = None
    
    def Corner(self): 
        ''' Wrapper to generate the corner plot of the ABC particle pool. The plot parameter 
        ranges are the prior ranges. The median of the parameters are marked in the plots.
        '''
        # list of parameters
        params = [
                'slope_gv', 
                'offset_gv', 
                'slope_fudge', 
                'offset_fudge', 
                'slope_tau', 
                'offset_tau'
                ]
        fig_name = ''.join([code_dir(), 
            'figure/', 'abc_step', str(self.t), '_', self.abcrun, '_weighted.png']) 
        self._thetas(self.theta, w=self.w, truths=self.med_theta, plot_range=self.prior_range, 
                parameters=params, 
                fig_name=fig_name)
        fig_name = ''.join([code_dir(), 
            'figure/', 'abc_step', str(self.t), '_', self.abcrun, '_unweighted.png']) 
        self._thetas(self.theta, truths=self.med_theta, plot_range=self.prior_range, 
                parameters=params, 
                fig_name=fig_name)
        return None

    def _thetas(self, theta, w=None, truths=None, plot_range=None, 
            parameters=['slope_gv', 'offset_gv', 'slope_fudge', 'offset_fudge', 'slope_tau', 'offset_tau'], 
            fig_name = None
            ):
        ''' Actually runs the corner plots for the ABC pool thetas. 
        '''
        par_labels = [] 
        for par in parameters: 
            if par == 'slope_gv': 
                par_labels.append(r'$\mathtt{A_{gv}}$')
            if par == 'offset_gv': 
                par_labels.append(r'$\delta_\mathtt{gv}$')
            if par == 'slope_fudge': 
                par_labels.append(r'$\mathtt{A_{fudge}}$')
            if par == 'offset_fudge': 
                par_labels.append(r'$\delta_\mathtt{fudge}$')
            if par == 'slope_tau': 
                par_labels.append(r'$\mathtt{A}_\tau$')
            if par == 'offset_tau': 
                par_labels.append(r'$\delta_\tau$')

        if w is not None: 
            # weighted theta
            fig = corner.corner(
                    theta,
                    weights=w.flatten(),
                    truths=truths,
                    truth_color='#ee6a50',
                    labels=par_labels,
                    label_kwargs={'fontsize': 25},
                    range=plot_range,
                    quantiles=[0.16,0.5,0.84],
                    show_titles=True,
                    title_args={"fontsize": 12},
                    plot_datapoints=True,
                    fill_contours=True,
                    levels=[0.68, 0.95],
                    color='b',
                    bins=20,
                    smooth=1.0)
        else: 
            # weighted theta
            fig = corner.corner(
                    theta,
                    truths=truths,
                    truth_color='#ee6a50',
                    labels=par_labels,
                    label_kwargs={'fontsize': 25},
                    range=plot_range,
                    quantiles=[0.16,0.5,0.84],
                    show_titles=True,
                    title_args={"fontsize": 12},
                    plot_datapoints=True,
                    fill_contours=True,
                    levels=[0.68, 0.95],
                    color='b',
                    bins=20,
                    smooth=1.0)

        if fig_name is not None: 
            plt.savefig(fig_name)
            plt.close()
        
        return None

    def Ssfr(self, subsample_thetas=False): 
        ''' Plot the SSFR distribution of the median paramater values of the 
        ABC particle pool.
        '''
        sfinherit_kwargs, abcrun_flag = ReadABCrun(self.abcrun)
        
        ssfr_plot = None
        if subsample_thetas: 
            N_theta = len(self.theta)
            i_subs = np.random.choice(N_theta, size=20)

            for i_sub in i_subs: 
                gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = self.theta[i_sub]

                sim_kwargs = sfinherit_kwargs.copy()
                sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
                sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
                sim_kwargs['evol_prop']['tau'] = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

                bloodline = InheritSF(
                        1,
                        nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                        subhalo_prop=sim_kwargs['subhalo_prop'], 
                        sfr_prop=sim_kwargs['sfr_prop'], 
                        evol_prop=sim_kwargs['evol_prop'])
                # descendant
                desc = getattr(bloodline, 'descendant_snapshot1') 

                ssfr_plot = desc.plotSsfr(line_color='red', lw=2, alpha=0.25, ssfr_plot=ssfr_plot)

        # plot the median 'best-fit' theta
        gv_slope = self.med_theta[0]
        gv_offset = self.med_theta[1]
        fudge_slope = self.med_theta[2]
        fudge_offset = self.med_theta[3]
        tau_slope = self.med_theta[4]
        tau_offset = self.med_theta[5]

        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
        sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
        sim_kwargs['evol_prop']['tau'] = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}

        if self.descendant is None: 
            bloodline = InheritSF(
                    1,
                    nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                    subhalo_prop=sim_kwargs['subhalo_prop'], 
                    sfr_prop=sim_kwargs['sfr_prop'], 
                    evol_prop=sim_kwargs['evol_prop'])
            # descendant
            descendant = getattr(bloodline, 'descendant_snapshot1') 
            self.descendant = descendant
        else: 
            descendant = self.descendant

        fig_name = ''.join([code_dir(), 
            'figure/SSFR.abc_step', str(self.t), '_', self.abcrun, '.png']) 
        descendant.plotSsfr(line_color='red', line_width=4, 
                sfms_prop=sim_kwargs['sfr_prop']['sfms'], z=descendant.zsnap, 
                groupcat=True, ssfr_plot=ssfr_plot, savefig=fig_name)
        plt.close()
        return None 

    def QAplot(self, nsnap_descendant=1):
        ''' Quality assurance plots
        '''
        sfinherit_kwargs, abcrun_flag = ReadABCrun(self.abcrun)

        gv_slope = self.med_theta[0]
        gv_offset = self.med_theta[1]
        fudge_slope = self.med_theta[2]
        fudge_offset = self.med_theta[3]
        tau_slope = self.med_theta[4]
        tau_offset = self.med_theta[5]
        
        # tau slopes and offsets of random particles 
        tau_slopes = self.theta[:,4]
        tau_offsets = self.theta[:,5]

        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {
                'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
        sim_kwargs['evol_prop']['fudge'] = {
                'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
        sim_kwargs['evol_prop']['tau'] = {
                'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
        
        if self.descendant is not None and self.descendant.nsnap == nsnap_descendant: 
            descendant = self.descendant
        else: 
            bloodline = InheritSF(
                    nsnap_descendant,
                    nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                    subhalo_prop=sim_kwargs['subhalo_prop'], 
                    sfr_prop=sim_kwargs['sfr_prop'], 
                    evol_prop=sim_kwargs['evol_prop'])
            # descendant
            descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
            self.descendant = descendant
        
        fig_name = ''.join([code_dir(), 
            'figure/QAPlot',
            '.nsnap', str(nsnap_descendant), 
             '.abc_step', str(self.t), '_', self.abcrun, '.png']) 
        QAplot(descendant, sim_kwargs, fig_name=fig_name, taus=[tau_slopes, tau_offsets])

        return None 




if __name__=="__main__": 
    for tf in [10]:
        ppp = PlotABC(tf, abcrun='run1')
        #ppp.Corner()
        #ppp.Ssfr()
        ppp.QAplot(nsnap_descendant=1)
        #ppp.QAplot(nsnap_descendant=6)
        #ppp.QAplot(nsnap_descendant=5)
        #ppp.QAplot(nsnap_descendant=7)
