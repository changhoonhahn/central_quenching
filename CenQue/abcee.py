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
from inherit import Inherit

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
    if 'fqz_multi' in observables: 
        qfrac = Fq()
        M_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
        M_mid = 0.5 * (M_bin[:-1] + M_bin[1:])
        
        fq_out = [M_mid]
        for zz in [0.0502, 0.1581, 0.3412, 1.0833]: 
            fq_model = qfrac.model(M_mid, zz, lit='wetzel')
            fq_out += [fq_model]

        obvs.append(fq_out)
    
    if len(observables) == 1: 
        obvs = obvs[0]
    return obvs

def SimSummary(observables=['ssfr'], **sim_kwargs):
    ''' Summary statistics of the simulation. In our case, the simulation is Inherit 
    and the summary statistic is the SSFR distribution. 
    '''
    obvs = []

    if 'fqz03' in observables and 'fqz_multi' in observables: 
        raise ValueError

    nsnap_ds = [1]
    if 'fqz03' in observables: 
        nsnap_ds += [6]
    if 'fqz_multi' in observables: 
        nsnap_ds += [3, 6] 

    inh = Inherit(nsnap_ds, 
            nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
            subhalo_prop=sim_kwargs['subhalo_prop'], 
            sfr_prop=sim_kwargs['sfr_prop'], 
            evol_prop=sim_kwargs['evol_prop'])
    des_dict = inh() 

    # SSFR 
    des1 = des_dict['1']
    if 'ssfr' in observables: 
        bins, dist = des1.Ssfr()
        obvs.append([np.array(bins), np.array(dist)])

    if 'fqz03' in observables:  # fQ(nsnap = 6) 
        des6 = des_dict['6']

        qfrac = Fq()
        M_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
        M_mid = 0.5 * (M_bin[:-1] + M_bin[1:]) 

        sfq = qfrac.Classify(des6.mass, des6.sfr, des6.zsnap, 
                sfms_prop=sim_kwargs['sfr_prop']['sfms'])

        ngal, dum = np.histogram(des6.mass, bins=M_bin)
        ngal_q, dum = np.histogram(des6.mass[sfq == 'quiescent'], bins=M_bin)

        obvs.append([M_mid, ngal_q.astype('float')/ngal.astype('float')])

    if 'fqz_multi' in observables:  # fQ at multiple snapshots 
        qfrac = Fq()
        M_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
        M_mid = 0.5 * (M_bin[:-1] + M_bin[1:]) 

        fq_list = [M_mid]

        for nd in nsnap_ds: 
            des_tmp = des_dict[str(nd)]

            sfq = qfrac.Classify(des_tmp.mass, des_tmp.sfr, des_tmp.zsnap, 
                    sfms_prop=sim_kwargs['sfr_prop']['sfms'])

            ngal, dum = np.histogram(des_tmp.mass, bins=M_bin)
            ngal_q, dum = np.histogram(des_tmp.mass[sfq == 'quiescent'], bins=M_bin)
            fq_tmp = ngal_q.astype('float')/ngal.astype('float')

            fq_list += [fq_tmp]
        # ancesotr 
        sfq = qfrac.Classify(inh.ancestor.mass, inh.ancestor.sfr, inh.ancestor.zsnap, 
                sfms_prop=sim_kwargs['sfr_prop']['sfms'])
        ngal, dum = np.histogram(inh.ancestor.mass, bins=M_bin)
        ngal_q, dum = np.histogram(inh.ancestor.mass[sfq == 'quiescent'], bins=M_bin)
        fq_tmp = ngal_q.astype('float')/ngal.astype('float')
        fq_list += [fq_tmp]

        obvs.append(fq_list)
    
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
    elif name == 'updated': # updated after try0 
        # These priors are overall wider in order to 
        # explore more of the parameter space. 
        prior_min = [0.0, -0.4, -5., 0.5, -1.5, 0.01]
        prior_max = [2.0, 0.6, 0.0, 2.5, 0.5, 1.5]
    elif name == 'satellite':   # theta does not include taus
        prior_min = [0.0, -0.4, -5., 0.5]
        prior_max = [2.0, 0.6, 0.0, 2.5]
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
    # rho_SSFR  
    datum_ssfr = datum[0][1]
    simum_ssfr = simum[0][1]
    rho_ssfr = np.sum((datum_ssfr - simum_ssfr)**2)
    
    # rho_fq
    n_fq = len(datum[1])-1 
    rho_fq = 0.
    for i_fq in range(n_fq):
        datum_fq = datum[1][i_fq+1]
        simum_fq = simum[1][i_fq+1]
        rho_fq += np.sum((datum_fq - simum_fq)**2)
    return [rho_ssfr, rho_fq]

def RhoFq(simum, datum):
    ''' Distance metric that compares the L2 Norm of quiescent fraction at different 
    redhsifts    
    Parameters
    ----------
    datum : list
        List [bins, dist] where bins is a np.ndarray of the SSFR bin values and 
        dist is a np.ndarray that specifies the SSFR distribution of the group catalog.
    model : list 
        List [bins, dist] where bins is a np.ndarray of the SSFR bin values and 
        dist is a np.ndarray that specifies the SSFR distribution of the simulation.

    '''
    # rho_fq
    n_fq = len(datum)-1 
    rho_fq = 0.
    for i_fq in range(n_fq):
        datum_fq = datum[i_fq+1]
        simum_fq = simum[i_fq+1]
        rho_fq += np.sum((datum_fq - simum_fq)**2)
    return rho_fq

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
        subhalo_dict = {'scatter': 0.2, 'source': 'li-march'}
    elif isinstance(subhalo, str):  
        subhalo_dict = {'scatter': 0.2, 'source': subhalo}
    f.write(''.join(['scatter = ', str(subhalo_dict['scatter']), '\n']))
    f.write(''.join(['source = ', subhalo_dict['source'], '\n']))
    f.write('\n')
    
    # Quiescent fraction properties
    f.write('# Quiescent Fraction Properties \n') 
    if fq is None:          # quiescent fraction
        #fq = {'name': 'wetzel'}
        fq = {'name': 'cosmos_tinker'} 
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
        dutycycle = {'name': 'notperiodic'}
        #dutycycle = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3}
    f.write(''.join(['name = ', dutycycle['name'], '\n']))
    if dutycycle['name'] == 'newamp_squarewave': 
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
            'subhalo_prop': subhalo_dict, 
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
        restart=False, t_restart=None, eps_restart=None, **sim_kwargs):
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
            eps_val=eps_input, restart=restart, **sim_kwargs) 

    # Data 
    data_sum = DataSummary(observables=observables)
    # Priors
    prior_min, prior_max = PriorRange(prior_name)
    prior = abcpmc.TophatPrior(prior_min, prior_max)    # ABCPMC prior object

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
    elif observables == ['ssfr', 'fqz_multi']: 
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
        new_eps_str = str(pool.eps)+'\t'+str(pool.ratio)+'\n'
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


def SatelliteABC(T, eps_input, Npart=1000, prior_name='try0', observables=['fqz_multi'], abcrun=None, 
        restart=False, t_restart=None, eps_restart=None, **sim_kwargs):
    ''' ABC-PMC  

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
            eps_val=eps_input, restart=restart, **sim_kwargs) 

    # Data 
    data_sum = DataSummary(observables=observables)
    # Priors
    prior_min, prior_max = PriorRange(prior_name)
    prior = abcpmc.TophatPrior(prior_min, prior_max)    # ABCPMC prior object

    def Simz(tt):       # Simulator (forward model) 
        gv_slope = tt[0]
        gv_offset = tt[1]
        fudge_slope = tt[2]
        fudge_offset = tt[3]

        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
        sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
        sim_kwargs['evol_prop']['tau'] = {'name': 'satellite'}

        return SimSummary(observables=observables, **sim_kwargs)

    theta_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'CenQue_theta_t', str(pewl), '_', abcrun_flag, '.satellite.dat']) 
    w_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'CenQue_w_t', str(pewl), '_', abcrun_flag, '.satellite.dat']) 
    dist_file = lambda pewl: ''.join([code_dir(), 
        'dat/pmc_abc/', 'CenQue_dist_t', str(pewl), '_', abcrun_flag, '.satellite.dat']) 
    eps_file = ''.join([code_dir(), 
        'dat/pmc_abc/epsilon_', abcrun_flag, '.satellite.dat'])

    distfn = RhoFq
   
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
        new_eps_str = str(pool.eps)+'\t'+str(pool.ratio)+'\n'
        if eps_str != new_eps_str:  # if eps is different, open fiel and append 
            f = open(eps_file, "a")
            eps_str = new_eps_str
            f.write(eps_str)
            f.close()

        print("T:{0},ratio: {1:>.4f}".format(pool.t, pool.ratio))
        print eps(pool.t)

        # write theta, weights, and distances to file 
        np.savetxt(theta_file(pool.t), pool.thetas, 
            header='gv_slope, gv_offset, fudge_slope, fudge_offset')
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
    
        if prior_name != 'satellite': 
            theta_file = ''.join([
                code_dir(), 'dat/pmc_abc/', 'CenQue_theta_t', str(t), '_', abcrun, '.dat']) 
            w_file = ''.join([
                code_dir(), 'dat/pmc_abc/', 'CenQue_w_t', str(t), '_', abcrun, '.dat']) 
            self.satellite_run = False 
        else: 
            theta_file = ''.join([
                code_dir(), 'dat/pmc_abc/', 'CenQue_theta_t', str(t), '_', abcrun, '.satellite.dat']) 
            w_file = ''.join([
                code_dir(), 'dat/pmc_abc/', 'CenQue_w_t', str(t), '_', abcrun, '.satellite.dat']) 
            self.satellite_run = True 
        self.theta = np.loadtxt(theta_file)      # theta values 
        self.w = np.loadtxt(w_file)              # w values 
        
        self.med_theta = [np.median(self.theta[:,i]) for i in range(len(self.theta[0]))]
        
        # prior range 
        prior_min, prior_max = PriorRange(prior_name)
        self.prior_range = [(prior_min[i], prior_max[i]) for i in range(len(prior_min))]

        self.descendant = None
    
    def Corner(self, filename=None): 
        ''' Wrapper to generate the corner plot of the ABC particle pool. The plot parameter 
        ranges are the prior ranges. The median of the parameters are marked in the plots.
        '''
        # list of parameters
        if not self.satellite_run:
            params = [
                    'slope_gv', 
                    'offset_gv', 
                    'slope_fudge', 
                    'offset_fudge', 
                    'slope_tau', 
                    'offset_tau'
                    ]
        else: 
            params = [
                    'slope_gv', 
                    'offset_gv', 
                    'slope_fudge', 
                    'offset_fudge'
                    ]

        if filename is None: 
            fig_name = ''.join([code_dir(), 
                'figure/', 'abc_step', str(self.t), '_', self.abcrun, '_weighted.png']) 
        else: 
            fig_name = filename

        if 'weighted' not in fig_name:
            fig_name = '_weighted.pn'.join(fig_name.split('.pn'))

        self._thetas(self.theta, w=self.w, truths=self.med_theta, plot_range=self.prior_range, 
                parameters=params, 
                fig_name=fig_name)

        fig_name = 'unweighted'.join(fig_name.split('weighted'))
        # fig_name = ''.join([code_dir(), 
        #    'figure/', 'abc_step', str(self.t), '_', self.abcrun, '_unweighted.png']) 
        self._thetas(self.theta, truths=self.med_theta, plot_range=self.prior_range, 
                parameters=params, 
                fig_name=fig_name)
        return fig_name

    def _thetas(self, theta, w=None, truths=None, plot_range=None, 
            parameters=['slope_gv', 'offset_gv', 'slope_fudge', 'offset_fudge', 'slope_tau', 'offset_tau'], 
            fig_name = None
            ):
        ''' Actually runs the corner plots for the ABC pool thetas. 
        '''
        prettyplot()
        par_labels = [] 
        for par in parameters: 
            if par == 'slope_gv': 
                par_labels.append(r'$\mathtt{A_{gv}}$')
            if par == 'offset_gv': 
                par_labels.append(r'$\delta_\mathtt{gv}$')
            if par == 'slope_fudge': 
                par_labels.append(r'$\mathtt{A_{P_Q}}$')
            if par == 'offset_fudge': 
                par_labels.append(r'$\delta_\mathtt{P_Q}$')
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
                    truth_color='k', #'#ee6a50',
                    labels=par_labels,
                    label_kwargs={'fontsize': 25},
                    range=plot_range,
                    quantiles=[0.16,0.5,0.84],
                    show_titles=True,
                    title_args={"fontsize": 12},
                    plot_datapoints=True,
                    fill_contours=True,
                    levels=[0.68, 0.95],
                    color='#ee6a50',
                    bins=20,
                    smooth=1.0)
        else: 
            # unweighted theta
            fig = corner.corner(
                    theta,
                    truths=truths,
                    truth_color='k',
                    labels=par_labels,
                    label_kwargs={'fontsize': 30},
                    range=plot_range,
                    quantiles=[0.16,0.5,0.84],
                    show_titles=True,
                    title_args={"fontsize": 15},
                    plot_datapoints=True,
                    fill_contours=True,
                    levels=[0.68, 0.95],
                    color='#ee6a50',
                    bins=20,
                    smooth=1.0)

        if fig_name is not None: 
            plt.savefig(fig_name)
            plt.close()
        
        return None

    def Ssfr(self, subsample_thetas=False, filename=None ): 
        ''' Plot the SSFR distribution of the median paramater values of the 
        ABC particle pool.
        '''
        sfinherit_kwargs, abcrun_flag = ReadABCrun(self.abcrun)
        
        ssfr_plot = None
        if subsample_thetas: 
            N_theta = len(self.theta)
            i_subs = np.random.choice(N_theta, size=20)

            for i_sub in i_subs: 
                if not self.satellite_run: 
                    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = self.theta[i_sub]
                else: 
                    gv_slope, gv_offset, fudge_slope, fudge_offset = self.theta[i_sub]

                sim_kwargs = sfinherit_kwargs.copy()
                sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
                sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
                if not self.satellite_run: 
                    sim_kwargs['evol_prop']['tau'] = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
                else: 
                    sim_kwargs['evol_prop']['tau'] = {'name': 'satellite'}

                inh = Inherit([1], 
                        nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                        subhalo_prop=sim_kwargs['subhalo_prop'], 
                        sfr_prop=sim_kwargs['sfr_prop'], 
                        evol_prop=sim_kwargs['evol_prop'])
                des_dict = inh() 
                # descendant
                desc = des_dict['1'] 

                ssfr_plot = desc.plotSsfr(line_color='red', lw=2, alpha=0.25, ssfr_plot=ssfr_plot)

        # plot the median 'best-fit' theta
        gv_slope = self.med_theta[0]
        gv_offset = self.med_theta[1]
        fudge_slope = self.med_theta[2]
        fudge_offset = self.med_theta[3]
        if not self.satellite_run: 
            tau_slope = self.med_theta[4]
            tau_offset = self.med_theta[5]

        sim_kwargs = sfinherit_kwargs.copy()
        sim_kwargs['sfr_prop']['gv'] = {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
        sim_kwargs['evol_prop']['fudge'] = {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
        if not self.satellite_run: 
            sim_kwargs['evol_prop']['tau'] = {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
        else: 
            sim_kwargs['evol_prop']['tau'] = {'name': 'satellite'}

        inh = Inherit([1], 
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        des_dict = inh() 
        descendant = des_dict['1'] 
        self.descendant = descendant
    
        if filename is None: 
            fig_name = ''.join([code_dir(), 
                'figure/SSFR.abc_step', str(self.t), '_', self.abcrun, '.png']) 
        else: 
            fig_name = filename 

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
        if not self.satellite_run: 
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
        if not self.satellite_run: 
            sim_kwargs['evol_prop']['tau'] = {
                    'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
        else: 
            sim_kwargs['evol_prop']['tau'] = {'name': 'satellite'}
        
        inh = Inherit(nsnap_descendant, 
                nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
                subhalo_prop=sim_kwargs['subhalo_prop'], 
                sfr_prop=sim_kwargs['sfr_prop'], 
                evol_prop=sim_kwargs['evol_prop'])
        des_dict = inh() 
        for nd in nsnap_descendant: 
            fig_name = ''.join([code_dir(), 
                'figure/QAPlot',
                '.nsnap', str(nd), 
                 '.abc_step', str(self.t), '_', self.abcrun, '.png']) 
            QAplot(des_dict[str(nd)], sim_kwargs, fig_name=fig_name)#, taus=[tau_slopes, tau_offsets])

        return None 




if __name__=="__main__": 
    for tf in [7]:
        ppp = PlotABC(tf, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated')
        ppp = PlotABC(tf, abcrun='RHOssfrfq', prior_name='updated')
        ppp.Corner()
        ppp.Ssfr()
        ppp.QAplot(nsnap_descendant=[1, 6])
    
    #for run in ['multifq_wideprior_nosmfevo', 'multifq_wideprior_extremesmfevo']:
    #    ppp = PlotABC(8, abcrun=run)
    #for tf in [0, 1, 2, 3, 4, 5]: 
    #    ppp = PlotABC(tf, abcrun=)
    #    ppp.Ssfr()
    #    ppp.QAplot(nsnap_descendant=[1, 6])
