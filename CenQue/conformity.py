'''
Code to explore galaxy conformity. 
'''
import numpy as np 

# Local ---
from util import util as Util
from inherit import ABCInherit


def AgeRankCentrals(tf, abcrun=None, prior_name=None, acc_precentage=80.):  
    ''' Measure age of the central halos and then rank their SFR in terms of their age
    in bins of M*.
    '''
    # load Inherited 
    if abcrun is None: 
        raise ValueError
    if prior_name is None: 
        raise ValueError
    abcinh = ABCInherit(tf, abcrun=abcrun, prior_name=prior_name) 
    inh = abcinh.PickleLoad('all') 

    des = inh['1']   # descendants at snapshot = 1
    anc = inh['ancestor']   # ancestor 
    
    anc_index, des_index = Util.intersection_index(anc.snap_index, des.ancestor15)
    if not np.array_equal(des_index, range(len(des.halo_mass))): 
        raise ValueError
    
    # Mhalo history 
    des_Mhalo_history = (anc.Mhalo_evol[0])[anc_index]
    # Mhalo history ratio
    des_Mhalo_ratio = 10**(des_Mhalo_history - des.halo_mass[:,None])
    
    acc_precentage /= 100.
    des_Mhalo_ratio[np.where(des_Mhalo_ratio == 0.)] = 999.
    # younger than ancestor (halo accretes acc_percentage of its mass after nsnap_ancestor) 
    younger = np.where(np.amin(des_Mhalo_ratio, axis=1) <= acc_percentage)
    # older than ancestor (halo accretes acc_percentage of its mass before nsnap_ancestor) 
    older = np.where(np.amin(des_Mhalo_ratio, axis=1) > acc_percentage)
    
    # get the snapshot when halos accrete the accretion_precentage of its mass  
    accretion_nsnap = np.zeros(len(des_Mhalo_ratio))
    accretion_nsnap[younger] = np.abs(
            des_Mhalo_ratio[younger] - np.tile(acc_percentage, (len(younger[0]), des_Mhalo_ratio.shape[1]))
            ).argmin(axis=1) + 1
    accretion_nsnap[older] = 999.   # older than ancestor
    
    # star-forming galaxies
    sfgal = np.where(des.sfr_class == 'star-forming')
    
    dlogm = 0.1
    m_arr = np.arange(des.mass[sfgal].min(), des.mass[sfgal].max()+dlogm, dlogm)
    m_low = m_arr[:-1]
    m_high = m_arr[1:]

    for i_m in range(len(m_low)):
        massbin = np.where(
                (des.sfr_class == 'star-forming') & 
                (des.mass >= m_low[i_m]) & 
                (des.mass < m_high[i_m])
                )
        if len(massbin[0]) == 0: 
            continue

        bin_sfr = des.sfr[massbin]
        bin_nsnap = accretion_nsnap[massbin]





if __name__=='__main__': 
    AgeRankCentrals(7, abcrun='multifq_wideprior', prior_name='updated')
