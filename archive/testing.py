from treepm import subhalo_io 
from treepm import sham
from utilities import utility 
import numpy as np 

sub = subhalo_io.Treepm.read('subhalo', 250, zis=range(1,16))  # snapchats of galaxies from z~0.0 to z~1.0833
sham.assign(sub, 'm.star', scat=0, dis_mf=0.0, zis=range(1,16)) # assigns M* to the subhalos using SHAM

for i in [1]: 
    cen_index = utility.utility_catalog.indices_ilk(sub[i], ilk='cen') # finds the indices to the central galaxies for snapshot 1
    print 'number of central galaxies in snapshot', i, '=', len(cen_index)

    mstar = sub[i]['m.star']
    pos = sub[i]['pos']
    print len(pos)
    print pos
    print pos[0]
    print pos[0,0], pos[1,0], pos[2,0]
    print pos[0,0], pos[0,1], pos[0,2]
    print pos[:,0]
    print pos[0,:]
    print (sub[i]['m.star'])[cen_index]
    cen_mstar = mstar[cen_index] 
    print cen_mstar
    print np.min(cen_mstar), np.max(cen_mstar)

#    par_index = utility.utility_catalog.indices_tree(sub, 1, 2, cen_index) 
#    par_nonneg = par_index[ par_index >= 0 ] 
#    print len(cen_index), len(par_index), len(par_nonneg)

#mstarz1 = sub[2]['m.star']
#print len(par_nonneg), len(par_index)
#print 'negative index = ', len(par_index)-len(par_nonneg)
#print len(mstarz1)
#cen_mstarz1 = mstarz1[par_nonneg] 
#
#print cen_mstarz1
