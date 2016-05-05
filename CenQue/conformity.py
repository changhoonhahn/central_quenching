'''
Code to explore galaxy conformity. 
'''
import numpy as np 

def AgeRankCentrals(tf, abcrun=None, prior_name=None):  
    ''' Measure age of the central halos and then rank their SFR in terms of their age
    in bins of M*.
    '''
    if abcrun is None: 
        raise ValueError
    if prior_name is None: 
        raise ValueError
    abcinh = ABCInherit(tf, abcrun=abcrun, prior_name=prior_name) 
    inh = abcinh.PickleLoad('all') 

