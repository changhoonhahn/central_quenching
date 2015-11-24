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

    def __init__(self, nsnap_ancestor = 13, quiet = True, **kwargs):
        """
        Class that describes the galaxy ancestors and descendants of 
        CenQue snapshots. 
        """

        self.kwargs = kwargs.copy()
        self.quiet = quiet
        self.nsnap_ancestor = nsnap_ancestor
        #self.nsnap_descendant = nsnap_descendant

        self.ancestor_cq = None
        self.descendant_cq = None
        self.file_name = None

    def ancestor(self, cenque_type='sf_assigned', scatter=0.0, sf_prop={'name': 'average'}, clobber=False):
        """
        Specify ancestor of lineage with a CenQue obj at snapshot nsnap_ancestor. 
        """

        self.ancestor_cq = CenQue(n_snap = self.nsnap_ancestor, scatter = scatter)
        self.ancestor_cenque_type = cenque_type 
        self.ancestor_cq.cenque_type = cenque_type 
        self.mass_scatter = scatter
        
        if np.array([os.path.isfile(self.ancestor_cq.file()), not clobber]).all(): 
            print 'Reading ', self.ancestor_cq.file()
            self.ancestor_cq.readin()
        else: 
            self.ancestor_cq.import_treepm(self.nsnap_ancestor, scatter=scatter)

            if cenque_type == 'sf_assigned': 
                print 'Assigning SFR'
                self.ancestor_cq = assign_sfr(self.ancestor_cq, sf_prop=sf_prop)

            self.ancestor_cq.writeout()
        
        if self.ancestor_cenque_type == 'sf_assigned': 
            self.ancestor_sf_prop = self.ancestor_cq.sf_prop
            self.ancestor_fq_prop = self.ancestor_cq.fq_prop

        return None

    def file(self): 
        """ 
        File name of Lineage object. Specifies the prescription used for 
        assigning star formation properties to the ancestor, which is the only main parameter 
        """

        if self.ancestor_cenque_type == 'treepm_import':
            
            ancestor_str = ''

        elif self.ancestor_cenque_type == 'sf_assigned':
    
            # properties that are specified for SF assignment
            if 'ancestor_sf_prop' not in self.__dict__.keys():
                raise ValueError
                #self.sf_prop = {'name': 'average'}
            if 'ancestor_fq_prop' not in self.__dict__.keys():
                raise ValueError
                #self.fq_prop = {'name': 'wetzelsmooth'}
            
            ancestor_str = '_'
            if self.ancestor_fq_prop['name'] == 'wetzelsmooth': 
                ancestor_str += 'wetzelsmooth' 
            ancestor_str +='_fq'

            ancestor_str += '_'
            if self.ancestor_sf_prop['name'] == 'average': 
                ancestor_str += 'average'
            elif self.ancestor_sf_prop['name'] == 'average_noscatter': 
                ancestor_str += 'average_noscatter'
            else: 
                raise NotImplementedError
            ancestor_str += '_sfassign'
                
        else: 
            raise NameError

        lineage_filename = ''.join([
            'dat/lineage/', 
            'lineage_ancestor_', 
            str(self.nsnap_ancestor), 
            ancestor_str, 
            '_descendants',
            '_mass_scatter', 
            str(round(self.mass_scatter, 1)), 
            '.hdf5'
            ]) 

        return lineage_filename

    def writeout(self): 
        """ 
        Write ancestor and descendant CenQue objects to hdf5 file 
        """
        
        if self.file_name is None: 
            output_file = self.file()  
            print 'Writing ', output_file 
        else: 
            output_file = self.file_name 

        f = h5py.File(output_file, 'w')    
        
        cq_grps = ['ancestor']
        for i_snap in xrange(1, self.nsnap_ancestor):
            cq_grps.append('descendant_cq_snapshot'+str(i_snap))

        for l_grp in cq_grps: 
            grp = f.create_group(l_grp)

            if l_grp == 'ancestor': 
                grp_cq = self.ancestor_cq
            else: 
                grp_cq = getattr(self, l_grp)
            
            # data columsn
            n_cols = len(grp_cq.data_columns)      
    
            for column in grp_cq.data_columns:      
                column_attr = getattr(grp_cq, column)
                grp.create_dataset( column, data=column_attr )     # save to h5py data group 

            grp.create_dataset( 'data_columns', data=grp_cq.data_columns )   

            # save metadata 
            for metadatum in grp_cq.metadata: 
                if isinstance(getattr(grp_cq, metadatum), dict): 
                    grp.attrs[metadatum] = json.dumps(getattr(grp_cq, metadatum))
                else: 
                    grp.attrs[metadatum] = getattr(grp_cq, metadatum) 

            grp.create_dataset( 'metadata', data=grp_cq.metadata )   

        f.close()  

        return None 

    def readin(self, 
            nsnap_descendants, 
            ancestor_cenque_type = 'sf_assigned', 
            ancestor_sf_prop = {'name': 'average'},
            ancestor_fq_prop = {'name': 'wetzelsmooth'},
            scatter=0.0,
            clobber = False
            ): 
        """
        Read in the lineage h5py object 
        """

        if self.file_name is None: 

            self.ancestor_cenque_type = ancestor_cenque_type 
            if self.ancestor_cenque_type == 'sf_assigned':
                self.ancestor_sf_prop = ancestor_sf_prop
                self.ancestor_fq_prop = ancestor_fq_prop
            self.mass_scatter = scatter

            lineage_file = self.file()   

            if not np.array([os.path.isfile(lineage_file), not clobber]).all(): 
                self.ancestor(clobber = clobber) 
                self.descend() 
                self.writeout()
        else: 
            lineage_file = self.file_name
        print 'Reading ', lineage_file
        f = h5py.File(lineage_file, 'r')
        
        cq_grps = ['ancestor']
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
            
            if l_grp == 'ancestor': 
                self.ancestor_cq = grp_cq
            else: 
                setattr(self, l_grp, grp_cq) 
        
        f.close() 

        return None 

    def descend(self):
        """
        Find descendants by tracking parents and children through TreePM's halos. 
        (Really only has to run once)
        """
        
        if not self.ancestor_cq: 
            raise ValueError('specify ancestor')

        parent_cq = self.ancestor_cq 
    
        for i_snap in range(1, self.nsnap_ancestor)[::-1]:    
    
            # import CenQue object from TreePM
            child_cq = CenQue() 
            child_cq.import_treepm(i_snap, scatter=self.mass_scatter) 
    
            # remove galaxies below the min and max mass
            mass_bins = child_cq.mass_bins
            within_massbin = np.where((child_cq.mass > min(mass_bins.mass_low))) 
            child_cq.sample_trim(within_massbin)                   
    
            n_child = len(within_massbin[0])     
            
            if not self.quiet: 
                print n_parent, ' parents ', n_child, ' children'
        
            setattr(
                    child_cq, 
                    'ancestor'+str(self.nsnap_ancestor)+'_index', 
                    np.array([-999 for i in xrange(n_child)])
                    )

            # parent children match which assigns indices into dictionaries and then get dictionary values
            parents, children = intersection_index(parent_cq.snap_index, child_cq.parent)

            # Deal with parent to children inheritance (children inherit the parents' attributes)
            try: 
                ancestor_index = getattr(
                        parent_cq, 
                        'ancestor'+str(self.nsnap_ancestor)+'_index'
                        )[parents]
                getattr(
                        child_cq, 
                        'ancestor'+str(self.nsnap_ancestor)+'_index'
                        )[children] = ancestor_index

                if not self.quiet: 
                    print ''
                    print child_cq.nsnap, '---------' 
                    print 'Parent with ancestors ', len(np.where(parent_cq.ancestor13_index >= 0)[0]) 
                    print 'Children with ancestors ', len(np.where(child_cq.ancestor13_index >= 0)[0])
                    print len(np.where(parent_cq.ancestor13_index >= 0)[0]) -\
                            len(np.where(child_cq.ancestor13_index >= 0)[0]), ' ancestors lost'

            except AttributeError: 

                ancestor_index = parent_cq.snap_index[parents]
                getattr(
                        child_cq, 
                        'ancestor'+str(self.nsnap_ancestor)+'_index'
                        )[children] = ancestor_index

            setattr(self, 'descendant_cq_snapshot'+str(i_snap), child_cq)

            parent_cq = child_cq
        
        child_cq.data_columns.append('ancestor'+str(self.nsnap_ancestor)+'_index')
        child_ancestor_index = getattr(child_cq, 'ancestor'+str(self.nsnap_ancestor)+'_index') 

        has_ancestor = np.where(child_ancestor_index >= 0)
        final_descendant = child_ancestor_index[has_ancestor]
        
        ancestor, descend = intersection_index(
                self.ancestor_cq.snap_index, 
                final_descendant
                )

        for i_snap in xrange(2, self.nsnap_ancestor):    
            child = getattr(self, 'descendant_cq_snapshot'+str(i_snap))
            child.data_columns.append('ancestor'+str(self.nsnap_ancestor)+'_index') 
            child_ancestor_index = getattr(child, 'ancestor'+str(self.nsnap_ancestor)+'_index') 

            ancestor_prime, descend_prime = intersection_index(
                    child_ancestor_index,
                    final_descendant
                    )
            if not np.array_equal(descend, descend_prime):
                raise ValueError

            child.sample_trim(ancestor_prime)
            setattr(self, 'descendant_cq_snapshot'+str(i_snap), child)
        
        descendant = has_ancestor[0][descend]

        # trim both ancestor and descendants to only keep ones that have lineage
        child_cq.sample_trim(descendant)
        setattr(self, 'descendant_cq_snapshot1', child_cq)
        self.ancestor_cq.sample_trim(ancestor)
    
        # check to make sure that all the CenQue object data columns are ordered appropriately 
        #for i_snap in xrange(1, self.nsnap_ancestor):    
        #    child = getattr(self, 'descendant_cq_snapshot'+str(i_snap))
        #    print getattr(child, 'ancestor'+str(self.nsnap_ancestor)+'_index') 
        #
        #print self.ancestor_cq.snap_index
        
        return None

if __name__=="__main__": 
    for scat in [0.0, 0.2]:
        bloodline = Lineage(nsnap_ancestor = 20)
        bloodline.ancestor(scatter = scat, sf_prop={'name': 'average'}, clobber=True) 
        bloodline.descend() 
        bloodline.writeout()
