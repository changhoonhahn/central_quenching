''' 
ABC test runs for specific purposes 
Make sure everything is documented precisely 
'''
from abcee import ABC

# Test the SMF evolution's impact on the quenching parameter constraints 
obvs_list = ['ssfr', 'fqz_multi'] 
guess = [10. for i in range(len(obvs_list))]
ABC(20, guess, Npart=1000, prior_name='updated', observables=obvs_list, 
        abcrun='RHOssfrfq_TinkerFq_XtraSMF', subhalo='li-march-extreme')
