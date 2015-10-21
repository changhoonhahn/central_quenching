import numpy as np
import abcpmc
import matplotlib.pyplot as plt
from interruptible_pool import InterruptiblePool
import time
plt.switch_backend("Agg")

from ssfr import rho_ssfr_cq_evol
from sf_inherit import rho_fq_ssfr_descendant
#from evolve_lineage import rho_ssfr_lineage_evol


def distance(theta, Mrcut=18): 

    rho = rho_fq_ssfr_descendant(
        nsnap_descendant = 1, 
        nsnap_ancestor = 20, 
        pq_prop = {'slope': theta[0], 'yint':theta[1]}, 
        tau_prop = {
            'name': 'line', 
            'fid_mass': 11.1, 
            'slope': theta[2], 
            'yint': theta[3]
            },
        Mrcut=18)
    #rho = rho_ssfr_lineage_evol(
    #    start_nsnap = 13, 
    #    final_nsnap = 1, 
    #    tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': theta[0], 'yint': theta[1]}, 
    #    Mrcut=Mrcut
    #    )
    #rho = rho_ssfr_cq_evol(
    #    start_nsnap = 13, 
    #    final_nsnap = 1, 
    #    sf_prop = {'name': 'average'}, 
    #    fq_prop = {'name': 'wetzelsmooth'}, 
    #    tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': theta[0], 'yint': theta[1]}, 
    #    mass_evol = 'sham', 
    #    Mrcut=Mrcut
    #    ) 
    return rho

"""covariance matrix in abc sampler"""

def covariance(theta , w , type = 'weighted'):

    if type == 'neutral':

      return np.cov(theta)

    if type == 'normalized neutral':

      return np.corrcoef(theta)

    if type == 'weighted':
      mean = np.sum(theta*w[None,:] , axis = 1)/ np.sum(w)
      tmm  = theta - mean.reshape(theta.shape[0] , 1)
      sigma2 = 1./(w.sum()) * (tmm*w[None,:]).dot(tmm.T)
      return sigma2  


from scipy.stats import uniform
from scipy.stats import norm 

class Prior(object): 
    def __init__(self, prior_dict): 
        self.prior_dict = prior_dict.copy()
        self.ordered_keys = prior_dict.keys()
        self.ordered_keys.sort()
    
    def prior(self): 
        priorz = [] 
        for key in self.ordered_keys: 
            
            prior_key = self.prior_dict[key]
            
            if prior_key['shape'] == 'uniform': 
                
                loc = prior_key['min']
                scale = prior_key['max'] - prior_key['min']
                
                priorz.append( uniform(loc, scale))
            
            elif prior_key['shape'] == 'gauss':
                
                loc = prior_key['mean']
                scale = prior_key['stddev']
                
                priorz.append( norm(loc, scale) )
                
        return priorz

prior_dict = {
        'pq_slope' : {'shape': 'uniform', 'min': 0.0, 'max': 0.25}, 
        'pq_yint' : {'shape': 'uniform', 'min': 0.0, 'max': 0.25}, 
        'tau_slope' : {'shape': 'uniform', 'min': -1.5, 'max': 0.0}, 
        'tau_yint': {'shape': 'uniform', 'min': 0.1, 'max': 1.0}
        }
n_params = len(prior_dict.keys())
prior_obj = Prior(prior_dict) 

N_threads = 10
N_particles = 500 
N_iter = 100
eps0 = 5.0

def prior_sampler(): 
    """ Sample prior distribution and return theta_star 
    """
    theta_star = np.zeros(n_params)
    
    for i in xrange(n_params): 
        np.random.seed()
        theta_star[i] = prior_obj.prior()[i].rvs(size=1)[0]

    return theta_star

def pi_priors(tmp_theta): 
    for i in xrange(n_params): 
        try:
            p_theta *= prior_obj.prior()[i].pdf(tmp_theta[i])
        except UnboundLocalError: 
            p_theta = prior_obj.prior()[i].pdf(tmp_theta[i])
            
    return p_theta 

import corner 
plot_range, plot_labels = [], [] 
for key in prior_obj.ordered_keys: 
    plot_labels.append(key)
    plot_range.append([prior_dict[key]['min'], prior_dict[key]['max']])
#print plot_range

def plot_thetas(theta , w , t): 
    fig = corner.corner(
        theta.T, weights = w.flatten() ,
        plot_datapoints=True, fill_contours=False, levels=[0.68, 0.95], 
                color='b', bins=40, smooth=1.0, 
        range=plot_range, 
        labels = plot_labels
        )
    
    plt.savefig("figure/theta_t"+str(t)+".png")
    plt.close()
    np.savetxt("dat/pmc_abc/theta_t"+str(t)+".dat" , theta.T)
    
    np.savetxt("dat/pmc_abc/w_t"+str(t)+".dat" , w.T)

def initial_pool_sampling(i_particle): 
    """ Sample theta_star from prior distribution for the initial pool
    """
    theta_star = prior_sampler()
    rho = distance(theta_star)
    
    while rho > eps0: 
        
        theta_star = prior_sampler()
        rho = distance(theta_star)
        
    pool_list = [np.int(i_particle)]
    for i_param in xrange(n_params): 
        pool_list.append(theta_star[i_param])
    pool_list.append(1./np.float(N_particles))
    pool_list.append(rho)
    
    return np.array(pool_list)


def initial_pool():

    args_list = [i for i in xrange(N_particles)]

    if N_threads > 1: 
        pool = InterruptiblePool(processes = N_threads)
        mapfn = pool.map
        results = mapfn(initial_pool_sampling, args_list)
        
        pool.close()
        pool.terminate()
        pool.join()
    else: 
        results = [] 
        for arg in args_list:  	
            results.append(initial_pool_sampling(arg))
    
    results = np.array(results).T
    theta_t = results[1:n_params+1,:]
    w_t = results[n_params+1,:]
    rhos = results[n_params+2,:]
    sig_t = covariance(theta_t , w_t)
    
    return theta_t, w_t, rhos, sig_t


def weighted_sampling(theta, w): 
    """ Given array of thetas and their corresponding weights, sample
    """
    w_cdf = w.cumsum()/w.sum() # normalized CDF
    
    np.random.seed()
    rand1 = np.random.random(1)
    cdf_closest_index = np.argmin( np.abs(w_cdf - rand1) )
    closest_theta = theta[:, cdf_closest_index]
    
    return closest_theta

def better_multinorm(theta_stst, theta_before, cov): 
    n_par, n_part = theta_before.shape
    
    sig_inv = np.linalg.inv(cov)
    x_mu = theta_before.T - theta_stst

    nrmliz = 1.0 / np.sqrt( (2.0*np.pi)**n_par * np.linalg.det(cov))

    multinorm = nrmliz * np.exp(-0.5 * np.sum( (x_mu.dot(sig_inv[None,:])[:,0,:]) * x_mu, axis=1 ) )

    return multinorm


def importance_pool_sampling(args): 
    # args = [i_particle, theta_t_1, w_t_1, sig_t_1, eps_t]
    i_particle = args[0]
    theta_t_1 = args[1]
    w_t_1 = args[2]
    sig_t_1 = args[3]
    eps_t = args[4]
    
    theta_star = weighted_sampling(theta_t_1, w_t_1)
    
    np.random.seed()
    # perturbed theta (Double check)    
    theta_starstar = np.random.multivariate_normal( theta_star, sig_t_1, 1 )[0]
    rho = distance(theta_starstar)
    
    while rho > eps_t:
        theta_star = weighted_sampling(theta_t_1, w_t_1)
        theta_starstar = np.random.multivariate_normal( theta_star, sig_t_1, 1 )[0]
        rho = distance(theta_starstar)

    p_theta = pi_priors(theta_starstar)

    w_starstar = p_theta/np.sum( w_t_1 * better_multinorm(theta_starstar, theta_t_1, sig_t_1) )    
    
    pool_list = [np.int(i_particle)]
    for i_p in xrange(n_params): 
        pool_list.append(theta_starstar[i_p])
    pool_list.append(w_starstar)
    pool_list.append(rho)
    
    return pool_list 
    
def pmc_abc(N_threads = N_threads): 
    
    # initial pool
    theta_t, w_t, rhos, sig_t = initial_pool()
    t = 0 # iternation number
    
    plot_thetas(theta_t , w_t, t)
    
    while t < N_iter: 
        
        eps_t = np.percentile(rhos, 75)
        print 'New Distance Threshold Eps_t = ', eps_t
        
        theta_t_1 = theta_t.copy()
        w_t_1 = w_t.copy()
        sig_t_1 = sig_t.copy()
        
        args_list = [[i, theta_t_1, w_t_1, sig_t_1, eps_t] for i in xrange(N_particles)]
        if N_threads > 1: 
            pool = InterruptiblePool(processes = N_threads)
            mapfn = pool.map
            results = mapfn(importance_pool_sampling, args_list)
            pool.close()
            pool.terminate()
            pool.join()
        else: 
            results = [] 
            for args in args_list: 
                pool_sample = importance_pool_sampling(args)
                results.append( pool_sample )
        
                 
        results = np.array(results).T
        theta_t = results[1:n_params+1,:]
        w_t = results[n_params+1,:]
        rhos = results[n_params+2,:]
        sig_t = covariance(theta_t , w_t)
        
        t += 1
        
        plot_thetas(theta_t, w_t , t)

pmc_abc()

