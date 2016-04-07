'''

Module for ABC-PMC inference for the Central Quenching project

'''
import sys 
import time
import pickle
import numpy as np
from datetime import datetime 

# --- Local ---
from observations import GroupCat
from sf_inherit import InheritSF

import abcpmc
from abcpmc import mpi_util

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
    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    
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


def RhoSSFR(datum, model):
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
    model_dist = model[1]
    return np.sum((datum_dist - model_dist)**2)



def MakeABCrun(file_name=None, nsnap_start=15, subhalo=None, fq=None, sfms=None, dutycycle=None, mass_evol=None): 
    '''Pickle file that specifies all the choices of parameters and ABC settings 
    for record keeping purposes.
    '''
    if file_name is None: 
        file_name = ''.join(['abc/',
            'RENAME_ME_', str(datetime.now()).replace(' ','.'), 
            '.txt']) 
    pickle_name = file_name.replace('.txt', '.p')

    f = open(file_name, 'w')

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
    f.write(''.join(['frequency range = ', str(dutycycle['freq_range'][0]), ',', str(dutycycle['freq_range'][0]), '\n']))
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

def ReadABCrun(file_name): 
    ''' Read text file that specifies all the choices of parameters and ABC settings.
    '''
    if not isinstance(file_name, str):
        raise ValueError('ABC run file name has to be a string')
    if file_name.rsplit('.')[-1] == 'txt':
        pickle_name = ''.join(['.'.join(file_name.rsplit('.')[:-1]), '.p']) 
    elif file_name.rsplit('.')[-1] == 'p':
        pickle_name = file_name

    return pickle.load(open(pickle_name, 'r'))


def ABC(T, eps_val, N_part=1000, prior_name='first_try', abcrun=None):
    ''' ABC-PMC implementation. 

    Parameters
    ----------
    - T : Number of iterations 
    - eps_val : 
    - N_part : Number of particles
    - observables : list of observables. Options are 'nbar', 'gmf', 'xi'
    - data_dict : dictionary that specifies the observation keywords 
    '''
    # Priors
    prior_min, prior_max = PriorRange(prior_name)
    prior = abcpmc.TophatPrior(prior_min, prior_max)    # ABCPMC prior object

    # Read in ABC run file
    sfinherit_kwargs = ReadABCrun(abcrun)

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
        
        return SimSummary(**sim_kwargs)

    theta_file = lambda pewl: ''.join(['dat/pmc_abc/', 'CenQue_theta_t', str(pewl), '.dat']) 
    w_file = lambda pewl: ''.join(['dat/pmc_abc/', 'CenQue_theta_t', str(pewl), '.dat']) 

    #mpi_pool = mpi_util.MpiPool()
    abcpmc_sampler = abcpmc.Sampler(
            N=N_part,               # N_particles
            Y=DataSummary,          # data
            postfn=Simz,            # simulator 
            dist=RhoSSFR            # distance function  
            )
    #        pool=mpi_pool)  
    abcpmc_sampler.particle_proposal_cls = abcpmc.ParticleProposal

    #eps = abcpmc.ConstEps(T, [1.e13,1.e13])
    eps = abcpmc.ConstEps(T, eps_val)
    pools = []
    f = open("dat/pmc_abc/abc_tolerance.dat" , "w")
    f.close()
    eps_str = ''
    for pool in abcpmc_sampler.sample(prior, eps):
        #while pool.ratio > 0.01:
        new_eps_str = '\t'.join(eps(pool.t).astype('str'))+'\n'
        if eps_str != new_eps_str:  # if eps is different, open fiel and append 
            f = open("abc_tolerance.dat" , "a")
            eps_str = new_eps_str
            f.write(eps_str)
            f.close()

        print("T:{0},ratio: {1:>.4f}".format(pool.t, pool.ratio))
        print eps(pool.t)
        # plot theta
        #plot_thetas(pool.thetas, pool.ws , pool.t, 
        #        Mr=data_dict["Mr"], truths=data_hod, plot_range=prior_range, observables=observables)
        # write theta and w to file 
        np.savetxt(theta_file(pool.t), pool.thetas)
        np.savetxt(w_file(pool.t), pool.ws)
        eps.eps = np.median(np.atleast_2d(pool.dists), axis = 0)


        pools.append(pool)

    abcpmc_sampler.close()
    return pools



if __name__=="__main__": 
    MakeABCrun(file_name='test.txt')
    ABC(1, 10., N_part=10, prior_name='try0', abcrun='test.txt')
