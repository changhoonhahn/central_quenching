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

from galpop import CGPop 
from galpop import AssignSFR 
#from assign_sfr import assign_sfr
import util.cenque_utility as util

class Lineage(object): 
    def __init__(self, nsnap_ancestor=20, subhalo_prop=None, clobber=False, **kwargs):
        ''' Class that describes the Lineage object which contains CGPop class objects
        for each snapshots from the ancestor to snapshot 1 (z = 0). 

        Parameters
        ----------
        nsnap_ancestor : int 
            Snapshot number of ancestors. Use snapshot to redshift conversion table to 
            determine what redshift the snapshot corresponds to. 

        clobber : bool
            If True, recompute the ancestor CGPop snapshot along with all the child CGPop 
            snapshots. 
        '''
        self.kwargs = kwargs.copy()
        self.nsnap_ancestor = nsnap_ancestor    # snapshot of ancestor
        self.ancestor = None
        self.clobber=clobber

        self.sfr_prop = None            # properties used for SFR assign 
        self.subhalo_prop = subhalo_prop.copy()     # properties used for subhalo SHAM
        if 'nsnap_ancestor' not in self.subhalo_prop.keys():
            self.subhalo_prop['nsnap_ancestor'] = nsnap_ancestor
    
        # Read/Import central galaxy population at snapshot = nsnap_ancestor
        self.ancestor = CGPop(n_snap=self.nsnap_ancestor, subhalo_prop=self.subhalo_prop)
        if os.path.isfile(self.ancestor.File()): 
            print 'Reading ', self.ancestor.File()
            self.ancestor.Read()
        else: 
            self.ancestor.ImportSubhalo(self.nsnap_ancestor)
            self.ancestor.Write()

    def AssignSFR_ancestor(self, sfr_prop=None):
        ''' Assign SFR properties to the 'ancestor_cq' object using 
        mass_genesis, tsnap_genesis

        Notes
        -----
        ''' 
        # SFR properties
        if self.sfr_prop is None: 
            if sfr_prop is not None: 
                self.sfr_prop = sfr_prop
            else: 
                raise ValueError
        else: 
            if sfr_prop: 
                if self.sfr_prop == 'default': 
                    self.sfr_prop = sfr_prop
                elif self.sfr_prop != sfr_prop: 
                    raise ValueError

        if 'mass_genesis' not in self.ancestor.__dict__.keys(): 
            raise ValueError
        
        if 'subhalogrowth' not in sfr_prop.keys(): 
            # Star-Formation properties based on mass_genesis at redshift, 
            # zsnap_genesis
            sfr_class, sfr, ssfr, delta_sfr, avg_sfr = AssignSFR(
                    self.ancestor.mass_genesis, 
                    self.ancestor.zsnap_genesis, 
                    sfr_prop=self.sfr_prop
                    )
        else: 
            desc_obj = getattr(self, 
                        'descendant_snapshot'+str(sfr_prop['subhalogrowth']['nsnap_descendant']))
            sfr_class, sfr, ssfr, delta_sfr, avg_sfr = AssignSFR(
                    self.ancestor.mass_genesis, 
                    self.ancestor.zsnap_genesis, 
                    sfr_prop=self.sfr_prop, 
                    ancestor=self.ancestor,
                    descendant=desc_obj
                    )

        for prop in ['sfr', 'ssfr', 'sfr_class', 'avg_sfr', 'delta_sfr']: 
            if prop not in self.ancestor.data_columns: 
                if isinstance(self.ancestor.data_columns, list): 
                    self.ancestor.data_columns.append(prop) 
                elif isinstance(self.ancestor.data_columns, np.ndarray): 
                    self.ancestor.data_columns = np.append(self.ancestor.data_columns, np.array([prop]))
                else: 
                    ValueError

        self.ancestor.sfr = sfr
        self.ancestor.ssfr = ssfr
        self.ancestor.sfr_class = sfr_class 
        self.ancestor.avg_sfr = avg_sfr
        self.ancestor.delta_sfr = delta_sfr
        self.ancestor.sfr_prop = self.sfr_prop.copy()

        return None

    def File(self): 
        ''' File name of Lineage object. The name specifies the properties
        of the subhalo catalog. 
        ''' 
        # properties that are specified for SF assignment
        if 'subhalo_prop' not in self.__dict__.keys():
            raise ValueError

        file_spec_str = self._file_spec(subhalo_prop = self.subhalo_prop)
        
        lineage_filename = ''.join([
            'dat/lineage/', 
            'lineage', 
            '.ancestor', str(self.nsnap_ancestor), 
            file_spec_str, 
            '.hdf5']) 
        return lineage_filename

    def _file_spec(self, subhalo_prop=None): 
        ''' file specifier string. Specifies the subhalo properties of the lienage object
        '''
        if subhalo_prop is None: 
            raise ValueError
        subhalo_str = ''.join([
            '.subhalo',
            '_scatter', str(subhalo_prop['scatter']), 
            '_', subhalo_prop['source']
            ])
        return ''.join([subhalo_str])#, fq_str, sfms_str])

    def Write(self, filename=None): 
        ''' Write all the Lineage data recording the  ancestor and descendants objects 
        to an hdf5 file 
        '''
        if filename is None: 
            filename = self.File()

        f = h5py.File(filename, 'w')    
        # each data group represents an object
        cq_grps = ['ancestor']
        for i_snap in xrange(1, self.nsnap_ancestor):
            cq_grps.append('descendant_snapshot'+str(i_snap))

        for l_grp in cq_grps: 
            try: 
                grp_cq = getattr(self, l_grp)
            except AttributeError: 
                continue
            grp = f.create_group(l_grp)
            
            # data columsn
            for column in grp_cq.data_columns:      
                column_attr = getattr(grp_cq, column)
                #print column
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

    def Read(self, nsnap_descendants, filename=None, clobber=False): 
        ''' Read in lineage object from stored hdf5 file. More specifically, 
        read in specified snapshots of descendants. 
        
        Parameters
        ----------
        nsnap_descendants : list 
            List of descendant N_snaps to read. 

        subhalo_prop : dict
            Dictionary that specifies the subhalo properties 

        sfr_prop : dict
            Dictionary that specifies the star forming properties of the ancestor object. 
        '''
        if filename is not None: 
            f = h5py.File(filename, 'r')
            print 'Reading ', filename 
        else: 
            f = h5py.File(self.File(), 'r')
            print 'Reading ', self.File() 
        
        cq_grps = ['ancestor']
        for i_snap in nsnap_descendants:
            cq_grps.append('descendant_snapshot'+str(i_snap))
        
        for l_grp in cq_grps: 
            grp = f[l_grp]
            grp_obj = CGPop()

            # read in meta data first
            for i_meta, metadatum in enumerate(grp.attrs.keys()): 
                if isinstance((grp.attrs.values())[i_meta], str):
                    try: 
                        json_dict = json.loads((grp.attrs.values())[i_meta])
                        setattr(grp_obj, metadatum, json_dict) 
                    except ValueError: 
                        setattr(grp_obj, metadatum, (grp.attrs.values())[i_meta]) 
                else:  
                    setattr(grp_obj, metadatum, (grp.attrs.values())[i_meta]) 
            grp_obj.metadata = [str(key) for key in grp.attrs.keys()]

            for i_col, column in enumerate(grp.keys()): 
                setattr(grp_obj, column, grp[column][:])
            
            setattr(self, l_grp, grp_obj) 
        
        f.close() 
        return None 

    def Descend(self, subhalo_prop=None, quiet=True, clobber=False):
        ''' Find descendants by tracking parents and children through TreePM's halos. 
        '''
        if self.subhalo_prop is None:
            if subhalo_prop is not None: 
                self.subhalo_prop = subhalo_prop
            else:
                raise ValueError
    
        # import Central galaxy populations at all snapshots after nsnap_ancestor
        child_list = [] 
        for i_snap in range(1, self.nsnap_ancestor): 
            child = CGPop(n_snap=i_snap, subhalo_prop=self.subhalo_prop) 
            if np.array([os.path.isfile(child.File()), not self.clobber]).all(): 
                print 'Reading ', child.File()
                child.Read()
            else: 
                child.ImportSubhalo(i_snap)
                child.Write()
            child_list.append(child)

        anc_nsnap_genesis = np.repeat(-999, len(self.ancestor.snap_index))
        anc_tsnap_genesis = np.repeat(-999., len(self.ancestor.snap_index))
        anc_zsnap_genesis = np.repeat(-999., len(self.ancestor.snap_index))
        anc_mass_genesis = np.repeat(-999., len(self.ancestor.snap_index))
        anc_halomass_genesis = np.repeat(-999., len(self.ancestor.snap_index))

        anc_Msham = np.repeat(-999., len(self.ancestor.snap_index)*len(child_list)).reshape(len(self.ancestor.snap_index), len(child_list))
    
        # go down the snapshots from nsnap_ancestor and track subhalos
        for i_snap in range(1, self.nsnap_ancestor)[::-1]:    

            child = child_list[i_snap-1]
            ancestor_index = getattr(child, 'ancestor'+str(self.nsnap_ancestor))

            # has ancestors at nsnap_ancestor
            has_ancestor, has_descendant = util.intersection_index(
                    ancestor_index, self.ancestor.snap_index)
            print 'Snanpshot ', i_snap
            print 'Children with ancestors ', len(has_ancestor), \
                    ' All children ', len(child.snap_index)

            # save SHAM masses
            anc_Msham[has_descendant, i_snap-1] = child.mass[has_ancestor]
    
            # snapshot, t_cosmic, and redshift where the subhalo starts 
            # hosting a galaxy. Aslo the mass of the new galaxy
            nsnap_genesis = np.repeat(-999, len(child.snap_index))
            mass_genesis = np.repeat(-999., len(child.snap_index)) 
            halomass_genesis = np.repeat(-999., len(child.snap_index)) 
            ancs = ancestor_index[has_ancestor] # ancestor indices

    
            # go through higher redshift snapshots in order to determine when the
            # subhalo was first "started"
            for ii_snap in range(i_snap, self.nsnap_ancestor): 
                
                child_ii = child_list[ii_snap-1]
                ii_anc, anc_ii = util.intersection_index(getattr(child_ii, 'ancestor'+str(self.nsnap_ancestor)), ancs)
                
                massive = np.where(child_ii.mass[ii_anc] > 0.0)

                nsnap_genesis[has_ancestor[anc_ii[massive]]] = ii_snap
                mass_genesis[has_ancestor[anc_ii[massive]]] = child_ii.mass[ii_anc[massive]]
                halomass_genesis[has_ancestor[anc_ii[massive]]] = child_ii.halo_mass[ii_anc[massive]]
                anc_nsnap_genesis[has_descendant[anc_ii[massive]]] = ii_snap
                anc_mass_genesis[has_descendant[anc_ii[massive]]] = child_ii.mass[ii_anc[massive]]
                anc_halomass_genesis[has_descendant[anc_ii[massive]]] = child_ii.halo_mass[ii_anc[massive]]
    
            massive_ancestor = np.where(self.ancestor.mass[has_descendant] > 0.0)
            nsnap_genesis[has_ancestor[massive_ancestor]] = self.nsnap_ancestor 
            mass_genesis[has_ancestor[massive_ancestor]] = self.ancestor.mass[has_descendant[massive_ancestor]]
            halomass_genesis[has_ancestor[massive_ancestor]] = self.ancestor.halo_mass[has_descendant[massive_ancestor]]
            anc_nsnap_genesis[has_descendant[massive_ancestor]] = self.nsnap_ancestor 
            anc_mass_genesis[has_descendant[massive_ancestor]] = self.ancestor.mass[has_descendant[massive_ancestor]]
            anc_halomass_genesis[has_descendant[massive_ancestor]] = self.ancestor.halo_mass[has_descendant[massive_ancestor]]

            nonneg = np.where(nsnap_genesis[has_ancestor] > 0)
            tsnap_genesis = np.repeat(-999., len(nsnap_genesis))
            zsnap_genesis = np.repeat(-999., len(nsnap_genesis))
            tsnap_genesis[has_ancestor[nonneg]] = util.get_t_nsnap(nsnap_genesis[has_ancestor[nonneg]])
            zsnap_genesis[has_ancestor[nonneg]] = util.get_z_nsnap(nsnap_genesis[has_ancestor[nonneg]])

            #neg = np.where(nsnap_massive[has_ancestor] < 0)
            #print child.mass[has_ancestor[neg]]

            # trim sample to only keep galaxies that have ancestors at nsnap_ancestor and 
            # and 'starts' before nsnap.
            child.sample_trim(has_ancestor[nonneg])  
            setattr(child, 'nsnap_genesis', nsnap_genesis[has_ancestor[nonneg]])
            setattr(child, 'tsnap_genesis', tsnap_genesis[has_ancestor[nonneg]])
            setattr(child, 'zsnap_genesis', zsnap_genesis[has_ancestor[nonneg]])
            setattr(child, 'mass_genesis', mass_genesis[has_ancestor[nonneg]])
            setattr(child, 'halomass_genesis', halomass_genesis[has_ancestor[nonneg]])
            child.data_columns += ['nsnap_genesis', 'tsnap_genesis', 'zsnap_genesis', 'mass_genesis', 'halomass_genesis']

            setattr(self, 'descendant_snapshot'+str(i_snap), child)

        positive = np.where(anc_nsnap_genesis > 0)
        self.ancestor.sample_trim(positive[0])
        self.ancestor.data_columns += ['nsnap_genesis', 'tsnap_genesis', 'zsnap_genesis', 'mass_genesis', 'halomass_genesis', 'Msham_evol']

        anc_tsnap_genesis[positive] = util.get_t_nsnap(anc_nsnap_genesis[positive])
        anc_zsnap_genesis[positive] = util.get_z_nsnap(anc_nsnap_genesis[positive])

        setattr(self.ancestor, 'nsnap_genesis', anc_nsnap_genesis[positive])
        setattr(self.ancestor, 'tsnap_genesis', anc_tsnap_genesis[positive])
        setattr(self.ancestor, 'zsnap_genesis', anc_zsnap_genesis[positive])
        setattr(self.ancestor, 'mass_genesis', anc_mass_genesis[positive])
        setattr(self.ancestor, 'halomass_genesis', anc_halomass_genesis[positive])
        setattr(self.ancestor, 'Msham_evol', anc_Msham[positive, :])

        return None

"""
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
"""


if __name__=="__main__": 
    for nsnap in [10]:#, 15, 20]: 
        for scat in [0.0, 0.2]:
            start_time = time.time()
            bloodline = Lineage(nsnap_ancestor=nsnap, 
                    subhalo_prop={'scatter': scat, 'source': 'li-march'}, clobber=True)
            bloodline.Descend(clobber=True) 
            bloodline.Write()
            print 'lineage construction and write out takes ', (time.time() - start_time)/60.0
