'''
Simple wrapper for running ABC with commandline arguments 
'''
import sys 
from abcee import ABC

Niter = int(sys.argv[1])
print 'N_iterations = ', Niter
Npart = int(sys.argv[2])
print 'N_particle = ', Npart
ABC(Niter, 10., Npart=Npart, prior_name='try0', abcrun='firstrun')
