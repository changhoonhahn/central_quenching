'''

Module for ABC-PMC inference for the Central Quenching project

'''
import time
import pickle
import numpy as np
from datetime import datetime 

# --- Local ---
from util.util import code_dir
from observations import GroupCat
from sf_inherit import InheritSF

import corner
import abcpmc
from abcpmc import mpi_util

import matplotlib.pyplot as plt 

def DataSummary(Mrcut=18): 
    ''' Summary statistics of the data. In our case that is the 
    SSFR distribution of the SDSS group catalog.
    '''
    # Group Catalog object
    groupcat = GroupCat(Mrcut=Mrcut, position='central')
    # SSFR distribution of group catalog
    bins, dist = groupcat.Ssfr()   
    return [np.array(bins), np.array(dist)]

def SimSummary(**sim_kwargs):
    ''' Summary statistics of the simulation. In our case, the simulation is InheritSF 
    and the summary statistic is the SSFR distribution. 
    '''
    bloodline = InheritSF(
            1,
            nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
            subhalo_prop=sim_kwargs['subhalo_prop'], 
            sfr_prop=sim_kwargs['sfr_prop'], 
            evol_prop=sim_kwargs['evol_prop'])
    descendant = getattr(bloodline, 'descendant_snapshot1') 
    
    bins, dist = descendant.Ssfr()
    return [np.array(bins), np.array(dist)]
    
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

def MakeABCrun(abcrun=None, nsnap_start=15, subhalo=None, fq=None, sfms=None, dutycycle=None, mass_evol=None, 
        Niter=None, Npart=None, prior_name=None, eps_val=None): 
    '''Pickle file that specifies all the choices of parameters and ABC settings 
    for record keeping purposes.
    '''
    if abcrun is None: 
        abcrun = ''.join(['RENAME_ME_', str(datetime.today().date())])
    file_name = ''.join([code_dir(), 'dat/pmc_abc/run/', 'abcrun_', abcrun, '.txt']) 
    pickle_name = file_name.replace('.txt', '.p')

    f = open(file_name, 'w')

    # ABC properities
    f.write('# ABC Run Properties \n') 
    f.write('N_iter = '+str(Niter)+'\n')
    f.write('N_particles = '+str(Npart)+'\n')
    f.write('Prior Name = '+prior_name+'\n')
    f.write('Epsilon_0 = '+str(eps_val)+'\n')
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
    return None 

def ReadABCrun(abcrun): 
    ''' Read text file that specifies all the choices of parameters and ABC settings.
    '''
    if abcrun is None: 
        abcrun = ''.join(['RENAME_ME_', str(datetime.today().date())])
    pickle_name = ''.join([code_dir(), 'dat/pmc_abc/run/', 'abcrun_', abcrun, '.p']) 
    return pickle.load(open(pickle_name, 'r')), abcrun

def ABC(T, eps_val, Npart=1000, prior_name='try0', abcrun=None):
    ''' ABC-PMC implementation. 

    Parameters
    ----------
    T : (int) 
        Number of iterations

    eps_val : (float)
        Starting epsilon threshold value 

    N_part : (int)
        Number of particles

    prior_name : (string)
        String that specifies what prior to use.

    abcrun : (string)
        String that specifies abc run information 
    '''
    # output abc run details
    MakeABCrun(abcrun=abcrun, Niter=T, Npart=Npart, prior_name=prior_name, eps_val=eps_val) 

    # Data 
    data_sum = DataSummary()
    # Priors
    prior_min, prior_max = PriorRange(prior_name)
    prior = abcpmc.TophatPrior(prior_min, prior_max)    # ABCPMC prior object

    # Read in ABC run file
    sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)

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

        bins, dists = SimSummary(**sim_kwargs)
        
        return [bins, dists]

    theta_file = lambda pewl: ''.join([code_dir(), 'dat/pmc_abc/', 'CenQue_theta_t', str(pewl), '_', abcrun_flag, '.dat']) 
    w_file = lambda pewl: ''.join([code_dir(), 'dat/pmc_abc/', 'CenQue_w_t', str(pewl), '_', abcrun_flag, '.dat']) 
    eps_file = ''.join([code_dir(), 'dat/pmc_abc/epsilon_', abcrun_flag, '.dat'])
    
    try:
        mpi_pool = mpi_util.MpiPool()
        abcpmc_sampler = abcpmc.Sampler(
                N=Npart,                # N_particles
                Y=data_sum,             # data
                postfn=Simz,            # simulator 
                dist=RhoSSFR,           # distance function  
                pool=mpi_pool)  
    except AttributeError: 
        abcpmc_sampler = abcpmc.Sampler(
                N=Npart,                # N_particles
                Y=data_sum,             # data
                postfn=Simz,            # simulator 
                dist=RhoSSFR)           # distance function  
    abcpmc_sampler.particle_proposal_cls = abcpmc.ParticleProposal

    eps = abcpmc.ConstEps(T, eps_val)
    pools = []
    f = open(eps_file, "w")
    f.close()
    eps_str = ''
    for pool in abcpmc_sampler.sample(prior, eps):
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
        print pool.__dict__.keys()

        # write theta and w to file 
        np.savetxt(theta_file(pool.t), pool.thetas, 
	    header='gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset')
        np.savetxt(w_file(pool.t), pool.ws)
        print pool.dists
        #eps.eps = np.median(np.atleast_2d(pool.dists), axis = 0)
        eps.eps = np.median(pool.dists)
        print '----------------------------------------'
        pools.append(pool)

    abcpmc_sampler.close()
    return pools

# ABC Corner Plot 
class PlotABC(object): 
    def __init__(self, t, abcrun=None, prior_name='try0'): 
        '''
        '''
        if abcrun is None: 
            raise ValueError('specify ABC run string') 

        theta_file = ''.join([
            code_dir(), 'dat/pmc_abc/', 'CenQue_theta_t', str(t), '_', abcrun, '.dat']) 
        w_file = ''.join([
            code_dir(), 'dat/pmc_abc/', 'CenQue_w_t', str(t), '_', abcrun, '.dat']) 

        theta = np.loadtxt(theta_file)
        w = np.loadtxt(w_file)

        med_theta = [np.median(theta[:,i]) for i in range(len(theta[0]))]
    
        # prior range 
        prior_min, prior_max = PriorRange(prior_name)
        prior_range = [(prior_min[i], prior_max[i]) for i in range(len(prior_min))]
        
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
            'figure/', 'abc_step', str(t), '_', abcrun, '_weighted.png']) 
        self._thetas(theta, w=w, truths=med_theta, plot_range=prior_range, 
                parameters=params, 
                fig_name=fig_name)
        fig_name = ''.join([code_dir(), 
            'figure/', 'abc_step', str(t), '_', abcrun, '_unweighted.png']) 
        self._thetas(theta, truths=med_theta, plot_range=prior_range, 
                parameters=params, 
                fig_name=fig_name)


    def _thetas(self, theta, w=None, truths=None, plot_range=None, 
            parameters=['slope_gv', 'offset_gv', 'slope_fudge', 'offset_fudge', 'slope_tau', 'offset_tau'], 
            fig_name = None
            ):
        ''' Corner plots of input theta values 
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




if __name__=="__main__": 
    pass
