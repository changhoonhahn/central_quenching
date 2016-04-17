'''
Simple wrapper for running ABC with commandline arguments 
'''
import sys 
from abcee import ABC

restart = int(sys.argv[1])

if restart == 0: 
    Niter = int(sys.argv[2])
    print 'N_iterations = ', Niter
    Npart = int(sys.argv[3])
    print 'N_particle = ', Npart
    abcrun = sys.argv[4]
    print 'ABC run name = ', abcrun

    ABC(Niter, [10., 10.], Npart=Npart, prior_name='try0', observables=['ssfr', 'fqz_multi'], abcrun=abcrun)
    #ABC(Niter, [10.], Npart=Npart, prior_name='try0', observables=['ssfr'], abcrun=abcrun)

elif restart == 1:  
    Niter = int(sys.argv[2])
    print 'N_iterations = ', Niter
    Npart = int(sys.argv[3])
    print 'N_particle = ', Npart
    abcrun = sys.argv[4]
    print 'ABC run name = ', abcrun
    trestart = int(sys.argv[5])
    print 'Restarting abc from t = ', str(trestart)

    ABC(Niter, [4.22460181, 0.33794247], Npart=Npart, prior_name='try0', observables=['ssfr', 'fqz03'], 
            abcrun=abcrun, t_restart=trestart, restart=True)

    #ABC(Niter, [4.28186054117], Npart=Npart, prior_name='try0', observables=['ssfr'], 
    #        abcrun=abcrun, t_restart=trestart, restart=True)
