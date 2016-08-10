from abcee import FixedTauABC

obvs_list = ['fqz_multi']
guess = [10. for i in range(len(obvs_list))]
FixedTauABC(3, guess, Npart=10, prior_name='satellite', observables=obvs_list, abcrun='TestTestSatABC_TinkerFq')
