'''
Simple wrapper for running ABC with commandline arguments 
'''
import sys 
from abcee import ABC
from abcee import FixedTauABC

restart = int(sys.argv[1])

if restart == 0: 
    Niter = int(sys.argv[2])
    print 'N_iterations = ', Niter
    Npart = int(sys.argv[3])
    print 'N_particle = ', Npart
    obvs = int(sys.argv[4])
    if obvs == 0:       # SSFR only distance metric
        print 'ssfr ONLY'
        obvs_list = ['ssfr'] 
        guess = [10.] 
    elif obvs == 1:     # SSFR + fQ at z=0.34 distance metric 
        print 'ssfr + f_Q(z=0.34)'
        obvs_list = ['ssfr', 'fqz03'] 
    elif obvs == 2:     # SSFR + fQ at redshifts
        print 'ssfr + [f_Q(z=0.05), f_Q(z=0.16), f_Q(z=0.3), f_Q(z=1.)]'
        obvs_list = ['ssfr', 'fqz_multi'] 
    elif obvs == 3:     # only fQ at redshifts
        print '[f_Q(z=0.05), f_Q(z=0.16), f_Q(z=0.3), f_Q(z=1.)]'
        obvs_list = ['fqz_multi']
    else: 
        raise ValueError 
    prior = sys.argv[5]
    print 'ABC Prior name = ', prior
    if prior not in ['try0', 'updated', 'satellite', 'longtau']:
        raise ValueError("Prior can only be 'try0', 'updated' and 'satellite'")
    abcrun = sys.argv[6]
    print 'ABC run name = ', abcrun

    guess = [10. for i in range(len(obvs_list))]
    
    if prior in ['satellite', 'longtau']: 
        FixedTauABC(Niter, guess, Npart=Npart, prior_name=prior, observables=obvs_list, abcrun=abcrun)
    else:
        ABC(Niter, guess, Npart=Npart, prior_name=prior, observables=obvs_list, abcrun=abcrun)

elif restart == 1:  
    raise NotImplementedError('Reimplement... carefully.')

    Niter = int(sys.argv[2])
    print 'N_iterations = ', Niter
    Npart = int(sys.argv[3])
    print 'N_particle = ', Npart
    abcrun = sys.argv[4]
    print 'ABC run name = ', abcrun
    trestart = int(sys.argv[5])
    print 'Restarting abc from t = ', str(trestart)

    ABC(Niter, [4.22460181, 0.33794247], Npart=Npart, prior_name='try0', 
            observables=['ssfr', 'fqz03'], abcrun=abcrun, t_restart=trestart, 
            restart=True)
