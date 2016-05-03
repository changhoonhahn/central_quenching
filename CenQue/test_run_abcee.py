''' 
ABC test runs for specific purposes 
Make sure everything is documented precisely 
'''
from abcee import ABC

# Test the SMF evolution's impact on the quenching parameter constraints 
obvs_list = ['ssfr', 'fqz_multi'] 
guess = [10. for i in range(len(obvs_list))]
ABC(3, guess, Npart=10, prior_name='updated', observables=obvs_list, 
        abcrun='multifq_wideprior_nosmfevo', subhalo='constant-li')
