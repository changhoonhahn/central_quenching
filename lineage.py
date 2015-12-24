"""

Lineage CenQue object ancestor and descendants

"""
import os
import json 
import time
import h5py
import random 
import numpy as np
import warnings

from cenque import CenQue
from cenque import AssignSFR 
#from assign_sfr import assign_sfr
import util.cenque_utility as util

class Lineage(object): 
    def __init__(self, nsnap_ancestor = 20, **kwargs):
        '''
        Class that describes the Lineage object which contains CenQue class objects
        for snapshots from ancestor to all its decendants.        
        '''
        self.kwargs = kwargs.copy()
        self.nsnap_ancestor = nsnap_ancestor
        self.ancestor_cq = None
        self.file_name = None

        self.sfr_prop = None            # properties used for SFR assign 
        self.subhalo_prop = None        # properties used for subhalo SHAM

    def assign_sfr_ancestor(self, sfr_prop = None):
        ''' 
        Assign SFR properties to the 'ancestor_cq' object using mass_genesis, tsnap_genesis

        Notes
        -----
        ''' 
        if self.ancestor_cq is None: 
            raise ValueError('Run self.descend() method first')
        # SFR properties
        if self.sfr_prop is None: 
            if sfr_prop is not None: 
                self.sfr_prop = sfr_prop
            else: 
                raise ValueError
        else: 
            if sfr_prop and self.sfr_prop != sfr_prop: 
                raise ValueError
    
        # Star-Formation properties based on mass_genesis at redshift, zsnap_genesis
        if sfr_prop['sfr']['name'] == 'average':
            gal_type, sfr, ssfr, delta_sfr, avg_sfr = AssignSFR(
                    self.ancestor_cq.mass_genesis, self.ancestor_cq.zsnap_genesis, sfr_prop = self.sfr_prop
                    )
            self.ancestor_cq.data_columns += ['sfr', 'ssfr', 'gal_type', 'avg_sfr', 'delta_sfr']
            self.ancestor_cq.sfr = sfr
            self.ancestor_cq.ssfr = ssfr
            self.ancestor_cq.gal_type = gal_type 
            self.ancestor_cq.avg_sfr = avg_sfr
            self.ancestor_cq.delta_sfr = delta_sfr
        else: 
            raise NotImplementedError

        return None

    def file(self): 
        ''' 
        File name of Lineage object. Specifies the prescription used for 
        assigning star formation properties to the ancestor, which is the only main parameter 
        ''' 
        # properties that are specified for SF assignment
        if 'sfr_prop' not in self.__dict__.keys():
            raise ValueError
        if 'subhalo_prop' not in self.__dict__.keys():
            raise ValueError

        file_spec_str = self._file_spec(
                subhalo_prop = self.subhalo_prop, 
                sfr_prop = self.sfr_prop
                )
        
        lineage_filename = ''.join([
            'dat/lineage/', 
            'lineage_ancestor.snapshot', 
            str(self.nsnap_ancestor), 
            file_spec_str, 
            '_descendants.hdf5'
            ]) 

        return lineage_filename

    def _file_spec(self, subhalo_prop = None, sfr_prop  = None): 
        if subhalo_prop is None: 
            raise ValueError

        subhalo_str = ''.join([
            '_subhalo.scatter', str(subhalo_prop['scatter']), 
            '.source-', subhalo_prop['source']
            ])

        if sfr_prop is not None: 
            sfr_str = ''.join([
                '_sfr.fq-', sfr_prop['fq']['name'], 
                '.sfr-', sfr_prop['sfr']['name']])
        else: 
            sfr_str = ''

        return ''.join([subhalo_str, sfr_str])

    def writeout(self): 
        ''' 
        Write ancestor and descendant CenQue objects to hdf5 file 
        '''
        if self.file_name is None: 
            output_file = self.file()  
            print 'Writing ', output_file 
        else: 
            output_file = self.file_name 

        f = h5py.File(output_file, 'w')    
        
        cq_grps = ['ancestor_cq']
        for i_snap in xrange(1, self.nsnap_ancestor):
            cq_grps.append('descendant_cq_snapshot'+str(i_snap))

        for l_grp in cq_grps: 
            grp = f.create_group(l_grp)
            grp_cq = getattr(self, l_grp)
            
            # data columsn
            for column in grp_cq.data_columns:      
                column_attr = getattr(grp_cq, column)
                grp.create_dataset(column, data=column_attr)     # save to h5py data group 
            # save metadata 
            for metadatum in grp_cq.metadata: 
                if isinstance(getattr(grp_cq, metadatum), dict): 
                    grp.attrs[metadatum] = json.dumps(getattr(grp_cq, metadatum))
                else: 
                    grp.attrs[metadatum] = getattr(grp_cq, metadatum) 
            
            grp.create_dataset('data_columns', data=grp_cq.data_columns)   
            grp.create_dataset('metadata', data=grp_cq.metadata )   

        f.close()  
        return None 

    def readin(self, nsnap_descendants, subhalo_prop = None, sfr_prop = None, clobber = False): 
        ''' 
        Read in the lineage h5py object 
        '''
        if self.file_name is None: 
            if self.subhalo_prop is None: 
                self.subhalo_prop = subhalo_prop 
            else: 
                if subhalo_prop and  self.subhalo_prop != subhalo_prop: 
                    raise ValueError
            if self.sfr_prop is None: 
                self.sfr_prop = sfr_prop
            else: 
                if sfr_prop and self.sfr_prop != sfr_prop: 
                    raise ValueError

            lineage_file = self.file()   

            if not np.array([os.path.isfile(lineage_file), not clobber]).all(): 
                self.ancestor(clobber = clobber) 
                self.descend() 
                self.writeout()
        else: 
            lineage_file = self.file_name

        print 'Reading ', lineage_file
        f = h5py.File(lineage_file, 'r')
        
        cq_grps = ['ancestor_cq']
        for i_snap in nsnap_descendants:
            cq_grps.append('descendant_cq_snapshot'+str(i_snap))
        
        for l_grp in cq_grps: 
            grp = f[l_grp]
            grp_cq = CenQue()

            # read in meta data first
            for i_meta, metadatum in enumerate(grp.attrs.keys()): 
                if isinstance((grp.attrs.values())[i_meta], str):
                    try: 
                        json_dict = json.loads((grp.attrs.values())[i_meta])
                        setattr(grp_cq, metadatum, json_dict) 
                    except ValueError: 
                        setattr(grp_cq, metadatum, (grp.attrs.values())[i_meta]) 
                else:  
                    setattr(grp_cq, metadatum, (grp.attrs.values())[i_meta]) 
            grp_cq.metadata = [str(key) for key in grp.attrs.keys()]

            for i_col, column in enumerate(grp.keys()): 
                setattr(grp_cq, column, grp[column][:])
            
            setattr(self, l_grp, grp_cq) 
        
        f.close() 
        return None 

    def descend(self, subhalo_prop=None, quiet=True, clobber=False):
        '''
        Find descendants by tracking parents and children through TreePM's halos. 
        (Really only has to run once)
        '''
        if self.subhalo_prop is None:
            if subhalo_prop is not None: 
                self.subhalo_prop = subhalo_prop
            else:
                raise ValueError

        # Import subhalos for nsnap <= nsnap_ancestor 
        self.ancestor_cq = CenQue(n_snap = self.nsnap_ancestor, subhalo_prop = self.subhalo_prop)
        self.ancestor_cq.cenque_type = 'treepm_import' 
        if np.array([os.path.isfile(self.ancestor_cq.file()), not clobber]).all(): 
            print 'Reading ', self.ancestor_cq.file()
            self.ancestor_cq.readin()
        else: 
            self.ancestor_cq.import_treepm(self.nsnap_ancestor)
            self.ancestor_cq.writeout()

        child_cq_list = [] 
        for i_snap in range(1, self.nsnap_ancestor): 
            child_cq = CenQue(n_snap=i_snap, subhalo_prop=self.subhalo_prop) 
            child_cq.cenque_type = 'treepm_import'
            if np.array([os.path.isfile(child_cq.file()), not clobber]).all(): 
                print 'Reading ', child_cq.file()
                child_cq.readin()
            else: 
                child_cq.import_treepm(i_snap)
                child_cq.writeout()
            child_cq_list.append(child_cq)

        anc_nsnap_genesis = np.repeat(-999, len(self.ancestor_cq.snap_index))
        anc_tsnap_genesis = np.repeat(-999., len(self.ancestor_cq.snap_index))
        anc_zsnap_genesis = np.repeat(-999., len(self.ancestor_cq.snap_index))
        anc_mass_genesis = np.repeat(-999., len(self.ancestor_cq.snap_index))
    
        for i_snap in range(1, self.nsnap_ancestor)[::-1]:    

            child_cq = child_cq_list[i_snap-1]
            ancestor_index = getattr(child_cq, 'ancestor'+str(self.nsnap_ancestor))

            # has ancestors at nsnap_ancestor
            has_ancestor, has_descendant = util.intersection_index(ancestor_index, self.ancestor_cq.snap_index)
            print 'Snanpshot ', i_snap
            print 'Children with ancestors ', len(has_ancestor), ' All children ', len(child_cq.snap_index)
    
            nsnap_genesis = np.repeat(-999, len(child_cq.snap_index))
            tsnap_genesis = np.repeat(-999., len(child_cq.snap_index))
            zsnap_genesis = np.repeat(-999., len(child_cq.snap_index))
            mass_genesis = np.repeat(-999., len(child_cq.snap_index)) 
            ancs = ancestor_index[has_ancestor] # ancestor indices
    
            # go through higher redshift snapshots in order to determine when the
            # subhalo was first "started"
            for ii_snap in range(i_snap, self.nsnap_ancestor): 
                
                ii_child = child_cq_list[ii_snap-1]
                ii_anc, anc_ii = util.intersection_index(getattr(ii_child, 'ancestor'+str(self.nsnap_ancestor)), ancs)
                
                massive = np.where(ii_child.mass[ii_anc] > 0.0)

                nsnap_genesis[has_ancestor[anc_ii[massive]]] = ii_snap
                mass_genesis[has_ancestor[anc_ii[massive]]] = ii_child.mass[ii_anc[massive]]
                anc_nsnap_genesis[has_descendant[anc_ii[massive]]] = ii_snap
                anc_mass_genesis[has_descendant[anc_ii[massive]]] = ii_child.mass[ii_anc[massive]]
    
            massive_ancestor = np.where(self.ancestor_cq.mass[has_descendant] > 0.0)
            nsnap_genesis[has_ancestor[massive_ancestor]] = self.nsnap_ancestor 
            mass_genesis[has_ancestor[massive_ancestor]] = self.ancestor_cq.mass[has_descendant[massive_ancestor]]
            anc_nsnap_genesis[has_descendant[massive_ancestor]] = self.nsnap_ancestor 
            anc_mass_genesis[has_descendant[massive_ancestor]] = self.ancestor_cq.mass[has_descendant[massive_ancestor]]

            nonneg = np.where(nsnap_genesis[has_ancestor] > 0)
            tsnap_genesis[has_ancestor[nonneg]] = util.get_t_nsnap(nsnap_genesis[has_ancestor[nonneg]])
            zsnap_genesis[has_ancestor[nonneg]] = util.get_z_nsnap(nsnap_genesis[has_ancestor[nonneg]])

            #neg = np.where(nsnap_massive[has_ancestor] < 0)
            #print child_cq.mass[has_ancestor[neg]]

            # trim sample to only keep galaxies that have ancestors at nsnap_ancestor and 
            # and 'starts' before nsnap.
            child_cq.sample_trim(has_ancestor[nonneg])  
            setattr(child_cq, 'nsnap_genesis', nsnap_genesis[has_ancestor[nonneg]])
            setattr(child_cq, 'tsnap_genesis', tsnap_genesis[has_ancestor[nonneg]])
            setattr(child_cq, 'zsnap_genesis', zsnap_genesis[has_ancestor[nonneg]])
            setattr(child_cq, 'mass_genesis', mass_genesis[has_ancestor[nonneg]])
            child_cq.data_columns += ['nsnap_genesis', 'tsnap_genesis', 'zsnap_genesis', 'mass_genesis']

            setattr(self, 'descendant_cq_snapshot'+str(i_snap), child_cq)

        self.ancestor_cq.data_columns += ['nsnap_genesis', 'tsnap_genesis', 'zsnap_genesis', 'mass_genesis']

        positive = np.where(anc_nsnap_genesis > 0)
        anc_tsnap_genesis[positive] = util.get_t_nsnap(anc_nsnap_genesis[positive])
        anc_zsnap_genesis[positive] = util.get_z_nsnap(anc_nsnap_genesis[positive])

        setattr(self.ancestor_cq, 'nsnap_genesis', anc_nsnap_genesis)
        setattr(self.ancestor_cq, 'tsnap_genesis', anc_tsnap_genesis)
        setattr(self.ancestor_cq, 'zsnap_genesis', anc_zsnap_genesis)
        setattr(self.ancestor_cq, 'mass_genesis', anc_mass_genesis)

        return None

    def _descend_old(self, quiet=True):
        '''
        Find descendants by tracking parents and children through TreePM's halos. 
        (Really only has to run once)
        '''
        if not self.ancestor_cq: 
            raise ValueError('specify ancestor')
        parent_cq = self.ancestor_cq 
    
        for i_snap in range(1, self.nsnap_ancestor)[::-1]:    
            # import CenQue object from TreePM
            child_cq = CenQue() 
            child_cq.import_treepm(i_snap, subhalo_prop = self.subhalo_prop) 
        
            # remove galaxies below the min and max mass
            keep = np.where(
                    (child_cq.mass > np.min(child_cq.mass_bins.mass_low)) & 
                    (child_cq.mass <= np.max(child_cq.mass_bins.mass_high))
                    ) 
            child_cq.sample_trim(keep)
            n_child = len(keep[0])
            #setattr(child_cq, 'ancestor'+str(self.nsnap_ancestor)+'_index', np.repeat(-999, n_child))
            
            # parents who have children and children who have parents. The following 
            # matches the parents to children based on snap_index attribute of parent CenQue
            # and parent attribute of child CenQue 
            parents, children = util.intersection_index(parent_cq.snap_index, child_cq.parent)
            print 'Parents with children ', len(parents), ' All parents ', len(parent_cq.snap_index)
            print 'Children with parents ', len(children), ' All children ', len(child_cq.parent)

            # Deal with parent to children inheritance (children inherit the parents' attributes)
            ancestor_index_attr = ''.join(['ancestor', str(self.nsnap_ancestor), '_index'])
            try: 
                ancestor_index = getattr(parent_cq, ancestor_index_attr)[parents]
                getattr(child_cq, ancestor_index_attr)[children] = ancestor_index

                if not quiet: 
                    print ''
                    print child_cq.nsnap, '---------' 
                    print 'Parent with ancestors ', len(np.where(parent_cq.ancestor20_index >= 0)[0]) 
                    print 'Children with ancestors ', len(np.where(child_cq.ancestor20_index >= 0)[0])
                    print len(np.where(parent_cq.ancestor13_index >= 0)[0]) -\
                            len(np.where(child_cq.ancestor13_index >= 0)[0]), ' ancestors lost'
            except AttributeError:  # highest snapshot 
                ancestor_index = parent_cq.snap_index[parents]
                getattr(child_cq, ancestor_index_attr)[children] = ancestor_index

            setattr(self, 'descendant_cq_snapshot'+str(i_snap), child_cq)
            parent_cq = child_cq
        
        child_cq.data_columns.append(ancestor_index_attr)
        for i_snap in xrange(2, self.nsnap_ancestor):    
            child.data_columns.append(ancestor_index_attr) 
        '''
        # For every snapshot CenQue object, only keep galaxies with ancestors at ancestor_nsnap
        # n snapshot = 0 
        child_ancestor_index = getattr(child_cq, ancestor_index_attr) 
        has_ancestor = np.where(child_ancestor_index >= 0)
        final = child_ancestor_index[has_ancestor]
        
        ancestor, descend = intersection_index(self.ancestor_cq.snap_index, final)
        print 'Number of last descendants who have ancestors ', len(descend)
        print 'Number of last descendants who have ancestors ', len(child_ancestor_index)
        
        descendant = has_ancestor[0][descend]
        child_cq.sample_trim(descendant)
        setattr(self, 'descendant_cq_snapshot1', child_cq)
        self.ancestor_cq.sample_trim(ancestor)

        for i_snap in xrange(2, self.nsnap_ancestor):    
            child = getattr(self, 'descendant_cq_snapshot'+str(i_snap))
            child.data_columns.append(ancestor_index_attr) 
            child_ancestor_index = getattr(child, ancestor_index_attr) 

            ancestor_prime, descend_prime = intersection_index(child_ancestor_index, final)
            if not np.array_equal(descend, descend_prime):
                raise ValueError

            child.sample_trim(ancestor_prime)
            setattr(self, 'descendant_cq_snapshot'+str(i_snap), child)
        
        # check to make sure that all the CenQue object data columns are ordered appropriately 
        for i_snap in xrange(1, self.nsnap_ancestor):    
            child = getattr(self, 'descendant_cq_snapshot'+str(i_snap))
            #print getattr(child, 'ancestor'+str(self.nsnap_ancestor)+'_index') 
        #print self.ancestor_cq.snap_index
        return None
        '''



if __name__=="__main__": 
    for scat in [0.0, 0.2]:
        start_time = time.time()
        bloodline = Lineage(nsnap_ancestor = 20)
        bloodline.descend(subhalo_prop = {'scatter': scat, 'source': 'li-march'})#, clobber=True) 
        bloodline.assign_sfr_ancestor(sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}})
        bloodline.writeout()
        print 'lineage construction and write out takes ', (time.time() - start_time)/60.0
