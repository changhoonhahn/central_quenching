"""

Lineage CenQue object ancestor and descendants

"""
import os
import json 
import h5py
import random 
import numpy as np
import warnings

from cenque import CenQue
from assign_sfr import assign_sfr
from util.cenque_utility import intersection_index

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

    def ancestor(self, subhalo_prop = None, sfr_prop = None, clobber = False):
        ''' 
        Specify ancestor of lineage with a SFR assigned CenQue obj at snapshot nsnap_ancestor. 
        ''' 
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

        self.ancestor_cq = CenQue(
                n_snap = self.nsnap_ancestor, 
                subhalo_prop = self.subhalo_prop, 
                sfr_prop = self.sfr_prop)
        self.ancestor_cq.cenque_type = 'sf_assigned' 
        
        if np.array([os.path.isfile(self.ancestor_cq.file()), not clobber]).all(): 
            print 'Reading ', self.ancestor_cq.file()
            self.ancestor_cq.readin()
        else: 
            self.ancestor_cq.import_treepm(self.nsnap_ancestor)
            print 'Assigning SFR'
            self.ancestor_cq.assign_sfr(quiet=True)
            self.ancestor_cq.writeout()

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

    def descend(self, quiet=True):
        '''
        Find descendants by tracking parents and children through TreePM's halos. 
        (Really only has to run once)
        '''
        if not self.ancestor_cq: 
            raise ValueError('specify ancestor')
        parent_cq = self.ancestor_cq 
        
        child_cq_list = [] 
        for i_snap in range(1, self.nsnap_ancestor):    
            # import CenQue object from TreePM
            child_cq = CenQue() 
            child_cq.import_treepm(i_snap, subhalo_prop = self.subhalo_prop) 
            child_cq_list.append(child_cq)
    
        for i_snap in range(1, self.nsnap_ancestor)[::-1]:    
            # remove galaxies below the min and max mass
            child_cq = CenQue() 
            child_cq.import_treepm(i_snap, subhalo_prop = self.subhalo_prop) 
            #keep = np.where(
            #        (child_cq.mass > np.min(child_cq.mass_bins.mass_low)) & 
            #        (child_cq.mass <= np.max(child_cq.mass_bins.mass_high))
            #        ) 
            #child_cq.sample_trim(keep)
            
            #nozero = np.where(child_cq.mass > 0)
            ancestor_index = getattr(child_cq, 'ancestor'+str(self.nsnap_ancestor))
            has_ancestor = np.where(ancestor_index > 0)
            print 'Snanpshot ', i_snap
            print 'Children with ancestors ', len(has_ancestor[0]), ' All children ', len(child_cq.parent)
        
            for ii in range(i_snap+1, self.nsnap_ancestor)[::-1]:
                ii_ancestor = getattr(child_cq, 'ancestor'+str(ii))[has_ancestor]
                ancestors = range(len(child_cq_list[ii-1].snap_index))[np.in1d(child_cq_list[ii-1].snap_index, ii_ancestor)]
                #intersection_index(
                #        child_cq_list[ii-1].snap_index,
                #        ii_ancestor)
                print ii, child_cq_list[ii-1].nsnap
                print 'these should be the same', len(children), len(ii_ancestor)
                print child_cq_list[ii-1].mass[ancestors]
                #print len(np.where(ii_ancestor < 0.)[0])

            #child_cq.sample_trim(has_ancestor)  # trim sample
            ancestors, children = intersection_index(self.ancestor_cq.snap_index, child_cq.ancestor20[has_ancestor])
            print len(self.ancestor_cq.snap_index[ancestors]), len(self.ancestor_cq.snap_index)
            print len(child_cq.ancestor20[has_ancestor[0][children]]), len(child_cq.ancestor20[has_ancestor])
            child_cq.sample_trim(has_ancestor[0][children])  # trim sample

            setattr(self, 'descendant_cq_snapshot'+str(i_snap), child_cq)
            #parent_cq = child_cq

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
            parents, children = intersection_index(parent_cq.snap_index, child_cq.parent)
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
        bloodline = Lineage(nsnap_ancestor = 20)
        bloodline.ancestor(
                subhalo_prop = {'scatter': scat, 'source': 'li-march'}, 
                sfr_prop = {
                    'fq': {'name': 'wetzelsmooth'}, 
                    'sfr': {'name': 'average'}
                    }, 
                clobber=True) 
        bloodline.descend() 
        bloodline.writeout()
