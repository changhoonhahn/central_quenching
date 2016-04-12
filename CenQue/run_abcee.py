'''
Simple wrapper for running ABC with commandline arguments 
'''
import sys 
from abcee import ABC

Niter = int(sys.argv[1])
print 'N_iterations = ', Niter
Npart = int(sys.argv[2])
print 'N_particle = ', Npart
abcrun = sys.argv[3]
print 'ABC run name = ', abcrun
ABC(Niter, [10., 10.], Npart=Npart, prior_name='try0', observables=['ssfr', 'fqz03'], abcrun=abcrun)
