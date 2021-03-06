#!/usr/bin/env python
import os 
import sys
import subprocess
import test_observations as Obv
from lineage import Lineage 

data_dir = raw_input('Please enter local directory to dump data : ')
u_sure = raw_input('Are you sure you want '+data_dir+' as your local data directory? [y/n]') 

if u_sure == 'y': 
    try: 
        os.symlink(data_dir, 'dat')
    except OSError: 
        os.remove('dat')
        u_sure2 = raw_input('Overwrite symlink? [y/n]')
        if u_sure2 == 'y': 
            os.symlink(data_dir, 'dat')
else: 
    raise ValueError("Don't doubt yourself next time") 

for dir in ['wetzel_tree', 'lineage', 'InheritSF', 'pmc_abc', 'galpop', 'observations']: 
    if not os.path.exists('dat/'+dir):
        print 'creating directory dat/'+dir+'/'
        os.makedirs('dat/'+dir)
if not os.path.exists('dat/pmc_abc/run'): 
        print 'creating directory dat/pmc_abc/run/'
        os.makedirs('dat/pmc_abc/run')

# observations
print 'Downloading SDSS group catalog ...'
if not os.path.exists('dat/observations/clf_groups_M18_9.4_D360.prob'): 
    subprocess.call(['wget', 'http://physics.nyu.edu/~chh327/data/groupcat.tar', '-P', 'dat/observations/']) 
    subprocess.call(['tar', '-xvf', 'dat/observations/groupcat.tar', '-C', 'dat/observations/'])
    print 'Downloading SDSS+PRIMUS iSEDfit data ...'
    subprocess.call(['wget', 'http://physics.nyu.edu/~chh327/data/primussdss.tar', '-P', 'dat/observations/']) 
    subprocess.call(['tar', '-xvf', 'dat/observations/primussdss.tar', '-C', 'dat/observations/'])
print 'Building Group Catalog hdf5 files'
[Obv.BuildGroupCat(Mrcut=Mr, position='central') for Mr in [18, 19, 20]]

# subhalos 
subhalo_file = 'subhalo_sham.central.snapshot1.ancestor15.scatter0.0.li-march.hdf5'
if not os.path.exists('dat/wetzel_tree/'+subhalo_file): 
    subprocess.call(['wget', 'http://physics.nyu.edu/~chh327/data/subhalo_sham.ancestor15.li-march.tar', '-P', 'dat/wetzel_tree/']) 
    subprocess.call(['tar', '-xvf', 'dat/wetzel_tree/subhalo_sham.ancestor15.li-march.tar', '-C', 'dat/wetzel_tree/'])

# lineage 
for scat in [0.0, 0.2]: 
    bloodline = Lineage(nsnap_ancestor=15, subhalo_prop={'scatter': scat, 'source': 'li-march'}, clobber=True)
    bloodline.Descend(clobber=True) 
    bloodline.Write()
