'''

Using Andrew Wetzel code, assigns M* to subhalos using SHAM

'''

from treepm import subhalo_io 
from treepm import sham
from utilities import utility 
import numpy as np
import pyfits

def build_centrals_snapshot(): 
    ''' Build central catalogues using SHAMed TreePM subhalos 
    '''
    snap_range = range(1,16)    # snapshot range

    # read in TreePM Subhalo snapshots from z~0.0502 to z~1.0833
    sub = subhalo_io.Treepm.read('subhalo',250,zis=snap_range) 

    # assign M* to subhalos using SHAM
    sham.assign(sub, 'm.star', scat=0, dis_mf=0.0, zis=snap_range) 

    for i in snap_range: 
        # indices of central galaxies for snapshot i 
        cen_index = utility.utility_catalog.indices_ilk(sub[i], ilk='cen') 
        print 'Number of Central Subhalos in Snapshot', i, '=', len(cen_index) 

        mstar   = (sub[i]['m.star'])[cen_index] # M* of subhalo
        pos     = (sub[i]['pos'])[cen_index] # position of subhalo
        ilk     = (sub[i]['ilk'])[cen_index] # classification flag in case we want to consider centrals and virtual centrals separately 

        # keep track of child indices for subhalo
        if i == np.min(snap_range):     
            child = np.zeros(len(cen_index))
            child.astype(int)
        else: 
            child = utility.utility_catalog.indices_tree(sub,i, i-1, cen_index) 

        if i == np.max(snap_range):     # parent index of subhalo
            parent = np.zeros(len(cen_index))
            parent.astype(int)
        else: 
            parent = utility.utility_catalog.indices_tree(sub,i, i+1, cen_index) 

        # write to fits file
        c0 = pyfits.Column(name='index', format='K', array=cen_index)
        c1 = pyfits.Column(name='pos', format='3E', array=pos) 
        c2 = pyfits.Column(name='ilk', format='I', array=ilk)
        c3 = pyfits.Column(name='mass', format='E', array=mstar)
        c4 = pyfits.Column(name='child', format='K', array=child)
        c5 = pyfits.Column(name='parent', format='K', array=parent)
        subsham_cen = pyfits.new_table([c0, c1, c2, c3, c4, c5])
        
        print np.min(mstar), np.max(mstar), np.mean(mstar)
        print i,'number of parents that have child',len(child[child >=0])
        print i,'number of child that have parent', len(parent[parent >=0]) 
       
        subsham_cen_dir = 'dat/wetzel_tree/'
        subsham_cen_file = 'subhalo_sham_centrals_snapshot'+str(i)+'.fits'
    #    subsham_cen.writeto(subsham_cen_dir+subsham_cen_file,clobber=1)
    #    output_file = open(subsham_cen_dir+subsham_cen_file, 'w') 
    #    for ii in range(len(cen_index)): 
    #        output_file.write('%f\t%f\t%f\t%d\t%f\t%d\t%d\n' % (pos[ii,0],pos[ii,1],pos[ii,2],ilk[ii],mstar[ii],child[ii],parent[ii]))
    #    output_file.close()

if __name__=="__main__": 
    build_centrals_snapshot()
