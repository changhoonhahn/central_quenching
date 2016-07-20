'''


Create catalog of central galaxies at the snapshots of just
the SFMS. 


'''
import h5py 
import numpy as np 

# Local 
from inherit import Inherit
from abcee import PlotABC
from abcee import ReadABCrun 

from util import util as UT


def Build_SFMScentrals(tf=7, abcrun='RHOssfrfq_TinkerFq_Std', prior_name='updated'): 
    ''' Build the SFMS central galaxy catalog using the best-fit quenching 
    probability parameters, and quenching timescale parameters from CenQue. 

    Creates two separate catalogs. One for quiescent and quenched galaxies 
    that use SHAM masses. The other with star forming and quenching galaxies
    that will ultimatley have integrated stellar masses. 

    Parameters
    ----------
    tf : int 
        timestep of the ABC run

    abcrun : string
        String that specifies the ABC run. 

    prior_name : string 
        String that specifies the name of the prior 

    '''
    # Get the best-fit parameters from CenQue 
    ppp = PlotABC(tf, abcrun=abcrun, prior_name=prior_name)
    gv_slope, gv_offset, fudge_slope, fudge_offset, tau_slope, tau_offset = ppp.med_theta
    
    # Using the best-fit parameters run Inherit to generate a
    # CenQue catalog
    sfinherit_kwargs, abcrun_flag = ReadABCrun(abcrun)
    sim_kwargs = sfinherit_kwargs.copy()
    sim_kwargs['sfr_prop']['gv'] = \
            {'slope': gv_slope, 'fidmass': 10.5, 'offset': gv_offset}
    sim_kwargs['evol_prop']['fudge'] = \
            {'slope': fudge_slope, 'fidmass': 10.5, 'offset': fudge_offset}
    sim_kwargs['evol_prop']['tau'] = \
            {'name': 'line', 'slope': tau_slope, 'fid_mass': 11.1, 'yint': tau_offset}
    print sim_kwargs['subhalo_prop']
    inh = Inherit([1], 
            nsnap_ancestor=sim_kwargs['nsnap_ancestor'],
            subhalo_prop=sim_kwargs['subhalo_prop'], 
            sfr_prop=sim_kwargs['sfr_prop'], 
            evol_prop=sim_kwargs['evol_prop'])
    des_dict = inh() 
    des = des_dict['1']     # descendant 

    # ancestor
    anc = inh.ancestor 
    anc_snap_index = anc.snap_index
    succession, will = UT.intersection_index(
            getattr(des, 'ancestor'+str(sim_kwargs['nsnap_ancestor'])), 
            anc_snap_index
            )
    M_sham_hist = anc.Msham_evol[0][will]
    M_halo_hist = anc.Mhalo_evol[0][will]

    sfr_class = des.sfr_class   # star-forming or quiescent
    t_quench = des.t_quench     # quenching time (if not quenching t_quench = 999) 
    quenched = des.quenched

    tsnap_genesis = des.tsnap_genesis   # When the central subhalo enters the simulation

    M_sham = des.mass 

    # quiescent galaxies and *quenched* galaxies
    q_only = np.where(
            ((sfr_class == 'quiescent') & (t_quench == 999.)) | 
            ((t_quench != 999.) & (quenched == 1))
            )
    N_q = len(q_only[0]) 
    
    # star forming and quenching galaxies (galaxies still quenching)
    others = np.where(
            ((t_quench != 999.) & (quenched == 0)) |
            (sfr_class == 'star-forming')
            ) 
    N_others = len(others[0]) 

    # Make sure the accounting checks out! 
    if N_q + N_others != len(des.mass): 
        raise ValueError 

    # Write quenched catalog to hdf5 file 
    quenched_file = ''.join([
        '/data1/hahn/centralMS/cenque/', 
        'quenched.centrals.', 
        'tf', str(tf), 
        '.abc_', abcrun, 
        '.prior_', prior_name, 
        '.hdf5'])
    f = h5py.File(quenched_file, 'w')    
    grp = f.create_group('data')
    # first save meta data 
    # ABC metadata
    grp.attrs['tf'] = tf
    grp.attrs['abcrun'] = abcrun
    grp.attrs['prior'] = prior_name 
    # cenque metadata
    grp.attrs['nsnap_ancestor'] = str(sim_kwargs['nsnap_ancestor']) 
    subhalo_prop_str = ','.join([   # subhalo prop
        ':'.join([key, str(sim_kwargs['subhalo_prop'][key])]) 
        for key in sim_kwargs['subhalo_prop'].keys()])
    grp.attrs['subhalo'] = subhalo_prop_str    
    fq_prop_str = ','.join([        # quiescent fraction prop 
        ':'.join([key, str(sim_kwargs['sfr_prop']['fq'][key])]) 
        for key in sim_kwargs['sfr_prop']['fq'].keys()])
    grp.attrs['fq'] = fq_prop_str
    gv_prop_str = ','.join([        # green valley prop 
        ':'.join([key, str(sim_kwargs['sfr_prop']['gv'][key])]) 
        for key in sim_kwargs['sfr_prop']['gv'].keys()])
    grp.attrs['gv'] = gv_prop_str
    sfms_prop_str = ','.join([          # SFMS prop 
        ':'.join([key, str(sim_kwargs['sfr_prop']['sfms'][key])]) 
        for key in sim_kwargs['sfr_prop']['sfms'].keys()])
    grp.attrs['sfms'] = sfms_prop_str
    pq_prop_str = ','.join([            # quenching probably correction factor 
        ':'.join([key, str(sim_kwargs['evol_prop']['fudge'][key])]) 
        for key in sim_kwargs['evol_prop']['fudge'].keys()])
    grp.attrs['pq'] = pq_prop_str 
    tau_prop_str = ','.join([            # quenching probably correction factor 
        ':'.join([key, str(sim_kwargs['evol_prop']['tau'][key])]) 
        for key in sim_kwargs['evol_prop']['tau'].keys()])
    grp.attrs['tau'] = pq_prop_str 
    # hardcoded columns for the catalogs
    for col in ['mass', 'sfr', 'ssfr', 'sfr_class', 'halo_mass', 
            'tsnap_genesis', 'nsnap_genesis', 'mass_genesis', 'halomass_genesis', 
            'Msham_hist', 'Mhalo_hist']: 
        if col == 'mass': 
            grp.create_dataset(col, data = M_sham[q_only]) 
        elif col == 'sfr_class': 
            grp.create_dataset(col, data = sfr_class[q_only]) 
        elif col == 'Msham_hist': 
            grp.create_dataset(col, data = M_sham_hist[q_only]) 
        elif col == 'Mhalo_hist': 
            grp.create_dataset(col, data = M_halo_hist[q_only]) 
        else: 
            grp.create_dataset(col, data=getattr(des, col)[q_only])     # save to h5py data group 
    f.close()

    # Write Star forming and Quenching catalog 
    sfms_file = ''.join([
        '/data1/hahn/centralMS/cenque/', 
        'sfms.centrals.', 
        'tf', str(tf), 
        '.abc_', abcrun, 
        '.prior_', prior_name, 
        '.hdf5'])
    f = h5py.File(sfms_file, 'w')    
    grp = f.create_group('data')
    # first save meta data 
    grp.attrs['tf'] = tf
    grp.attrs['abcrun'] = abcrun
    grp.attrs['prior'] = prior_name 
    # cenque metadata
    grp.attrs['nsnap_ancestor'] = str(sim_kwargs['nsnap_ancestor']) 
    grp.attrs['subhalo_prop'] = subhalo_prop_str    
    grp.attrs['subhalo'] = subhalo_prop_str    
    grp.attrs['fq'] = fq_prop_str
    grp.attrs['gv'] = gv_prop_str
    grp.attrs['sfms'] = sfms_prop_str
    grp.attrs['pq'] = pq_prop_str 
    grp.attrs['tau'] = tau_prop_str 

    # hardcoded columns for the catalogs
    for col in ['mass', 'sfr', 'ssfr', 'sfr_class', 'halo_mass',
            'tsnap_genesis', 'nsnap_genesis', 'mass_genesis', 't_quench', 'halomass_genesis', 
            'Msham_hist', 'Mhalo_hist']: 
        if col == 'mass': 
            grp.create_dataset(col, data = M_sham[others]) 
        elif col == 'sfr_class': 
            grp.create_dataset(col, data = sfr_class[others]) 
        elif col == 'Msham_hist': 
            grp.create_dataset(col, data = M_sham_hist[others]) 
        elif col == 'Mhalo_hist': 
            grp.create_dataset(col, data = M_halo_hist[others]) 
        else: 
            grp.create_dataset(col, data=getattr(des, col)[others])     # save to h5py data group 
    f.close()

    return None



if __name__=='__main__': 
    Build_SFMScentrals(tf=8, abcrun='RHOssfrfq_TinkerFq_NOSMFevol', prior_name='updated')
