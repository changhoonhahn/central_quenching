import numpy as np
import corner 
import time
from scipy.stats import uniform
from scipy.stats import norm 
from interruptible_pool import InterruptiblePool
# --- Local --- 
from ssfr import rho_ssfr_cq_evol
from evolve_lineage import rho_ssfr_lineage_evol

def covariance(theta , w , type = 'weighted'):
    """ covariance matrix in abc sampler
    """

    if type == 'neutral':

      return np.cov(theta)

    if type == 'normalized neutral':

      return np.corrcoef(theta)

    if type == 'weighted':
      mean = np.sum(theta*w[None,:] , axis = 1)/ np.sum(w)
      tmm  = theta - mean.reshape(theta.shape[0] , 1)
      sigma2 = 1./(w.sum()) * (tmm*w[None,:]).dot(tmm.T)
      return sigma2  

class Prior(object): 
    def __init__(self, prior_dict): 
        """Class that describes the priors based on specified prior dictionary 
        """

        self.prior_dict = prior_dict.copy()
        
        self.n_params = len(self.prior_dict.keys())   # number of parameters

        self.ordered_keys = prior_dict.keys()
        self.ordered_keys.sort()

        self.priorz = self.prior()
    
    def prior(self): 
        priorz = [] 
        for key in self.prior_dict.keys(): 
            
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

    def sampler(self): 
        """ Samples the prior distribution to return a vector theta of parameters
        """
    
        theta = np.zeros(self.n_params)

        for i in xrange(self.n_params): 

            np.random.seed()

            theta[i] = self.priorz[i].rvs(size=1)[0]
                                  
        return theta
    
    def pi_priors(self, theta_in): 
        """ Get the product of the prior probability distributions for each particle 
        represented by theta_in

        $\prod p_i(theta_i)$

        """

        for i in xrange(self.n_params): 
            try:
                p_theta *= self.priorz[i].pdf(theta_in[i])
            except UnboundLocalError: 
                p_theta = self.priorz[i].pdf(theta_in[i])
        
        return p_theta 

def distance(theta, Mrcut=18): 
    rho = rho_ssfr_lineage_evol(
        start_nsnap = 13, 
        final_nsnap = 1, 
        tau_prop = {'name': 'line', 'fid_mass': 10.75, 'slope': theta[0], 'yint': theta[1]}, 
        Mrcut=Mrcut
        )

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

def initial_pool_sampling(args):
    """ Sample theta_star from prior distribution for the initial pool
    """
    i_particle, prior_obj, eps0, N_particles = args
    theta_star = prior_obj.sampler()
    rho = distance(theta_star)
    print theta_star, rho
    
    while rho > eps0: 
        
        theta_star = prior_obj.sampler()
        rho = distance(data, theta_star)
        print theta_star, rho
        
    pool_list = [np.int(i_particle)]
    for i_param in xrange(prior_obj.n_params): 
        pool_list.append(theta_star[i_param])
    pool_list.append(1./np.float(N_particles))
    pool_list.append(rho)
    
    return np.array(pool_list)

def initial_pool(prior_obj, eps0, N_particles, N_threads=1):
    """ Initial Pool
    """
    
    args_list = [[i, prior_obj, eps0, N_particles] for i in xrange(N_particles)]
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
    theta_t = results[1:prior_obj.n_params+1,:]
    w_t = results[prior_obj.n_params+1,:]
    rhos = results[prior_obj.n_params+2,:]
    sig_t = covariance(theta_t , w_t)
    
    return theta_t, w_t, rhos, sig_t

def importance_pool_sampling(args): 
    # args = [i_particle, theta_t_1, w_t_1, sig_t_1, eps_t]
    i_particle = args[0]
    prior_obj = args[1]
    theta_t_1 = args[2]
    w_t_1 = args[3]
    sig_t_1 = args[4]
    eps_t = args[5]
    
    theta_star = weighted_sampling(theta_t_1, w_t_1)
    
    # perturbed theta (Double check)    
    np.random.seed()
    theta_starstar = np.random.multivariate_normal( theta_star, sig_t_1, 1 )[0]
    rho = distance(theta_starstar)
    
    while rho > eps_t:
        theta_star = weighted_sampling(theta_t_1, w_t_1)
        theta_starstar = np.random.multivariate_normal( theta_star, sig_t_1, 1)[0]
        rho = distance(theta_starstar)
    
    p_theta = prior_obj.pi_priors(theta_starstar)

    w_starstar = p_theta/np.sum( w_t_1 * better_multinorm(theta_starstar, theta_t_1, sig_t_1) )    
    
    pool_list = [np.int(i_particle)]
    for i_p in xrange(prior_obj.n_params): 
        pool_list.append(theta_starstar[i_p])
    pool_list.append(w_starstar)
    pool_list.append(rho)
    
    return pool_list 
    
def pmc_abc(prior_dict, N_particles=100, N_iter=30, eps0=20.0, N_threads = 1): 
    """
    """
    prior_obj = Prior(prior_dict)

    # initial pool
    theta_t, w_t, rhos, sig_t = initial_pool(prior_obj, eps0, N_particles, N_threads=N_threads)
    t = 0 # iternation number
    
    #plot_thetas(theta_t , w_t, prior_dict, t)
    
    while t < N_iter: 
        
        eps_t = np.percentile(rhos, 75)
        print 'New Distance Threshold Eps_t = ', eps_t
        
        theta_t_1 = theta_t.copy()
        w_t_1 = w_t.copy()
        sig_t_1 = sig_t.copy()
    

        args_list = [[i, prior_obj, theta_t_1, w_t_1, sig_t_1, eps_t] for i in xrange(N_particles)]

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
        theta_t = results[1:prior_obj.n_params+1,:]
        w_t = results[prior_obj.n_params+1,:]
        rhos = results[prior_obj.n_params+2,:]

        sig_t = covariance(theta_t , w_t)
        
        t += 1
        
        plot_thetas(theta_t , w_t, prior_dict, t)

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

def plot_thetas(theta, w, prior_dict, iterno, basename="abc_run",
                truths=None, fig_name=None, **kwargs):
    """ Plot the parameter values of the particles
    """

    np.savetxt("{0}_{1}_thetas.dat".format(basename, str(iterno)), theta)
    np.savetxt("{0}_{1}_ws.dat".format(basename, str(iterno)), w)

    if 'range' in kwargs.keys():

        param_range = kwargs['range']

    else:
        ordered_key = prior_dict.keys()
        ordered_key.sort()

        param_range = []

        for key in ordered_key:
            if prior_dict[key]['shape'] == 'uniform':
                param_range.append(
                        (prior_dict[key]['min'], prior_dict[key]['max'])
                        )

            elif prior_dict[key]['shape'] == 'gaussian':
                param_range.append(
                        (prior_dict[key]['mean'] - 3.0 * prior_dict[key]['stddev'],
                            prior_dict[key]['mean'] + 3.0 * prior_dict[key]['stddev'])
                        )
            else:
                raise NotImplementedError()

    if 'labels' in kwargs.keys():
        plt_labels = kwargs['labels']
    else:
        plt_labels = ordered_key

    fig = corner.corner(
            theta.T,
            weights = w,
            smooth = True,
            truths = truths,
            truth_color = 'red',
            range = param_range,
            labels = plt_labels,
            levels = [0.68, 0.95],
            )

    if fig_name == None:
        fig_name = 'blah.png'

    fig.savefig(fig_name)
    plt.close()

if __name__=="__main__": 
    prior_dict = {
            'slope': {'shape': 'uniform', 'min': -1.5, 'max': 0.0}, 
            'yint': {'shape': 'uniform', 'min': 0.1, 'max': 1.0}
            }

    pmc_abc(prior_dict, N_particles=100, N_iter=30, eps0=5.0, N_threads = 5)
