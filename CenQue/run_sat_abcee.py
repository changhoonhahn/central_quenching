from abcee import SatelliteABC

obvs_list = ['fqz_multi']
guess = [10. for i in range(len(obvs_list))]
SatelliteABC(20, guess, Npart=500, prior_name='satellite', observables=obvs_list, abcrun='SatABC_TinkerFq')
