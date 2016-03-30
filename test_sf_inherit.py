'''

Module to test the sf_inherit module. 

'''

import numpy as np
import time

from gal_prop import Fq 
from gal_prop import SMF
from lineage import Lineage
from sf_inherit import InheritSF 
from util.util import get_t_nsnap

from plotting.plots import PlotSMF

# --- plotting --- 
import matplotlib.pyplot as plt
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors

pretty_colors = prettycolors()

def TestInheritSF(nsnap_descendant, nsnap_ancestor=20, 
        subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, 
        sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
        evol_prop = {
            'pq': {'slope': 0.05, 'yint': 0.0}, 
            'tau': {'name': 'line', 'fid_mass': 10.75, 'slope': -0.6, 'yint': 0.6}, 
            'sfr': {'dutycycle': {'name': 'notperiodic'}}, 
            'mass': {'name': 'sham'}}, 
        ssfr=True, fq=False, tau=False, smhm=False, sfms=False, smf=False):

    bloodline = InheritSF(nsnap_descendant, nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop, 
            sfr_prop=sfr_prop, evol_prop=evol_prop)
    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant))
    
    if evol_prop['mass']['name'] == 'sham': 
        massevol_str = 'sham'
    else: 
        massevol_str = ''.join([evol_prop['mass']['type'], evol_prop['mass']['name']])

    lineage_str = ''.join([ 
        '.ancestor', str(nsnap_ancestor), 
        '.descendant', str(nsnap_descendant), 
        '.subhalo_Mscatter', str(round(subhalo_prop['scatter'], 1)), 
        '.sfr_dutycycle_', evol_prop['sfr']['dutycycle']['name'], 
        '.massevol_', massevol_str])

    attr_list = []
    if ssfr: 
        attr_list.append('ssfr')
    if fq: 
        attr_list.append('fq')
    if tau: 
        attr_list.append('tau')
    if sfms: 
        attr_list.append('sfms')
    if smhm: 
        attr_list.append('smhm')
    if smf: 
        attr_list.append('smf')

    for attr in attr_list: 
        fig_name = ''.join([
            'figure/test/', 
            'TestInheritSF_', attr, lineage_str, '.png'])

        if attr == 'ssfr':                          # SSFR
            descendant.plotSsfr(line_color='red', line_width=4, 
                    groupcat=True, savefig=fig_name)
        elif attr == 'fq':                          # Quiescent Fraction 
            descendant.plotFq(model=True, line_color='r', label=None,
                    savefig=fig_name)
        elif attr == 'tau':                         # Quenching Timescale
            descendant.plotTau(evol_prop['tau'], savefig=fig_name)
        elif attr == 'smf':                         # Stellar Mass Function
            descendant.plotSMF(savefig=fig_name)
        elif attr == 'sfms':                        # Star Forming Main Sequence 
            descendant.plotSFMS(bovyplot=True, model=True, savefig=fig_name)
        elif attr == 'smhm':                # Stellar Mass - Halo Mass 
            descendant.plotSMHM(color=descendant.nsnap, savefig=fig_name)
    return None

def TestInheritSF_ABC(nsnap_descendant, nsnap_ancestor=20, n_step=29, 
        subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, 
        sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
        evol_prop = {
            'sfr': {'dutycycle': {'name': 'notperiodic'}}, 
            'mass': {'name': 'sham'}}, 
        ssfr=True, fq=False, tau=False, smhm=False, sfms=False, smf=False):
    ''' Plot InheritSF output for median posterior from AB
    '''
    # quenching probabilyt and tau properties from ABC posterior
    med_theta = abc_posterior_median(n_step)
    pq_dict = {'slope': med_theta[0], 'yint': med_theta[1]}
    tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}
    evol_prop['pq'] = pq_dict
    evol_prop['tau'] = tau_dict

    TestInheritSF(nsnap_descendant, nsnap_ancestor=nsnap_ancestor, 
            subhalo_prop = subhalo_prop, 
            sfr_prop = sfr_prop,
            evol_prop = evol_prop, 
            ssfr=ssfr, fq=fq, tau=tau, smhm=smhm, sfms=sfms, smf=smf)
    return None

def AnecstorPlots(abc_step=29, n_descendant=1, **sfinh_kwargs): 
    ''' FQ for the SF Inherited Ancestor galaxy population. This test is mainly 
    to see what the discrepancy between the empirical SF/Q classifcation and the
    parameterized FQ model. Ideally, the model and the 
    '''
    sfinherit_file = InheritSF_file(n_descendant, abc_step=abc_step, **sfinh_kwargs)
    bloodline = Lineage(
            nsnap_ancestor=sfinh_kwargs['nsnap_ancestor'], 
            subhalo_prop=sfinh_kwargs['subhalo_prop']
            )
    bloodline.Read([n_descendant], filename=sfinherit_file)

    ancestor = getattr(bloodline, 'ancestor') 
    ancestor.sfr_prop = sfinh_kwargs['sfr_prop']
    
    # FQ
    fig_file = lambda prop: ''.join(['figure/test/', 
        'Ancestor', prop.upper(), '.',
        '.'.join((sfinherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.png'])
    ancestor.plotFq(model=True, savefig=fig_file('fq'))      
    ancestor.plotSFMS(sfqcut=True, gal_class='all', bovyplot=False, sigSFR=False, model=False, savefig=fig_file('sfms'))

    return None

def DescendantSSFR(nsnap_descendant, nsnap_ancestor=20, n_step=29, 
    subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, 
    sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfms': {'name': 'linear'}},
    evol_prop = {'sfr': {'dutycycle': {'name': 'notperiodic'}}, 'mass': {'name': 'sham'}}): 
    ''' Halo Mass Function composition at z = z_final by initial halo mass at 
    nsnap_ancestor. 
    '''
    # import SF inherited Lineage
    sf_inherit_file = InheritSF_file(
            nsnap_descendant, 
            nsnap_ancestor=nsnap_ancestor, 
            abc_step=n_step, 
            subhalo_prop=subhalo_prop, 
            sfr_prop=sfr_prop, 
            evol_prop=evol_prop
            )
    bloodline = Lineage(
            nsnap_ancestor=nsnap_ancestor, 
            subhalo_prop=subhalo_prop
            )
    bloodline.Read([nsnap_descendant], filename=sf_inherit_file)
    # descendant
    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 

    fig_name = ''.join([
        'figure/test/', 
        'SSFR.',
        '.'.join((sf_inherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.png'])

    descendant.plotSsfr(line_color='red', line_width=4, 
            groupcat=True, savefig=fig_name)
    return None

def DescendantHMF_composition(nsnap_descendant, nsnap_ancestor=20, n_step=29, 
    subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, 
    sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfms': {'name': 'linear'}},
    evol_prop = {'sfr': {'dutycycle': {'name': 'notperiodic'}}, 'mass': {'name': 'sham'}}): 
    ''' Halo Mass Function composition at z = z_final by initial halo mass at 
    nsnap_ancestor. 
    '''
    # import SF inherited Lineage
    sf_inherit_file = InheritSF_file(
            nsnap_descendant, 
            nsnap_ancestor=nsnap_ancestor, 
            abc_step=n_step, 
            subhalo_prop=subhalo_prop, 
            sfr_prop=sfr_prop, 
            evol_prop=evol_prop
            )
    bloodline = Lineage(
            nsnap_ancestor=nsnap_ancestor, 
            subhalo_prop=subhalo_prop
            )
    bloodline.Read([nsnap_descendant], filename=sf_inherit_file)
    # descendant
    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    smf_plot = PlotSMF()
    mf = SMF()

    # SMF composition based on their initial mass at nsnap_ancestor 
    started_here = np.where(descendant.nsnap_genesis == nsnap_ancestor)
    start_mass = descendant.halomass_genesis[started_here]
    
    mass_bins = np.arange(10., 17.5, 0.5)
    mass_bin_low = mass_bins[:-1]
    mass_bin_high = mass_bins[1:]

    for i_m in range(len(mass_bin_low)): 

        mbin = np.where( 
                (start_mass > mass_bin_low[i_m]) & 
                (start_mass <= mass_bin_high[i_m])
                )

        mass, phi = mf._smf(descendant.halo_mass[started_here[0][mbin]], 
                m_arr= np.arange(10., 17.1, 0.1))

        try: 
            phi_tot += phi 
            smf_plot.sub.fill_between(mass, phi_tot - phi, phi_tot, 
                    color=pretty_colors[i_m], 
                    label=str(round(mass_bin_low[i_m], 2))+'-'+str(round(mass_bin_high[i_m], 2)))
        except UnboundLocalError: 
            phi_tot = phi 
            smf_plot.sub.fill_between(mass, np.zeros(len(phi)), phi, 
                    color=pretty_colors[i_m], 
                    label=str(round(mass_bin_low[i_m], 2))+'-'+str(round(mass_bin_high[i_m], 2)))

    #smf_plot.set_axes()
    smf_plot.sub.set_yscale('log') 
    smf_plot.sub.set_xlim([10.0, 17.0])
    smf_plot.sub.legend(loc='upper right') 
    mass_evol_str = '.mass_evol'

    fig_file = ''.join([
        'figure/test/', 
        ''.join((sf_inherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.HMF_composition', 
        mass_evol_str, '.png'])
    smf_plot.save_fig(fig_file)
    return None

def DescendantSMF_composition(nsnap_descendant, mass_evol=True, nsnap_ancestor=20, abc_step=29, 
    subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, 
    sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
    evol_prop = {'sfr': {'dutycycle': {'name': 'notperiodic'}}, 'mass': {'name': 'sham'}}): 
    
    sf_inherit_file = InheritSF_file(nsnap_descendant, nsnap_ancestor=nsnap_ancestor, 
        abc_step=abc_step, subhalo_prop=subhalo_prop, sfr_prop=sfr_prop, evol_prop=evol_prop)

    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read([nsnap_descendant], filename=sf_inherit_file)
    
    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    smf_plot = descendant.plotSMF()
    
    mf = SMF()

    if mass_evol: 
        # SMF composition based on their initial mass at nsnap_ancestor 
        started_here = np.where(descendant.nsnap_genesis == nsnap_ancestor)
        start_mass = descendant.mass_genesis[started_here]
        
        mass_bins = np.arange(7., 12.5, 0.5)
        mass_bin_low = mass_bins[:-1]
        mass_bin_high = mass_bins[1:]

        for i_m in range(len(mass_bin_low)): 

            mbin = np.where( 
                    (start_mass > mass_bin_low[i_m]) & 
                    (start_mass <= mass_bin_high[i_m])
                    )

            mass, phi = mf._smf(descendant.mass[started_here[0][mbin]])

            try: 
                phi_tot += phi 
                smf_plot.sub.fill_between(mass, phi_tot - phi, phi_tot, 
                        color=pretty_colors[i_m], 
                        label=str(round(mass_bin_low[i_m], 2))+'-'+str(round(mass_bin_high[i_m], 2)))
            except UnboundLocalError: 
                phi_tot = phi 
                smf_plot.sub.fill_between(mass, np.zeros(len(phi)), phi, 
                        color=pretty_colors[i_m], 
                        label=str(round(mass_bin_low[i_m], 2))+'-'+str(round(mass_bin_high[i_m], 2)))
    else:
        for isnap in range(2, nsnap_ancestor+1)[::-1]:
            started_here = np.where(descendant.nsnap_genesis == isnap)
            mass, phi = mf._smf(descendant.mass[started_here])

            try: 
                phi_tot += phi 
                smf_plot.sub.fill_between(mass, phi_tot - phi, phi_tot, color=pretty_colors[isnap], label='Nsnap='+str(isnap))
            except UnboundLocalError: 
                phi_tot = phi 
                smf_plot.sub.fill_between(mass, np.zeros(len(phi)), phi, color=pretty_colors[isnap], label='Nsnap='+str(isnap))
    smf_plot.set_axes()
    mass_evol_str = ''
    if mass_evol: 
        mass_evol_str = '.mass_evol'

    fig_file = ''.join([
        'figure/test/', 
        ''.join((sf_inherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.SMF_composition', 
        mass_evol_str, '.png'])
    smf_plot.save_fig(fig_file)
    return None

def Compare_DescendantSMF_composition(nsnap_descendant, kwarg_list, n_step=29):
    ''' Compare the M* distribution at z_final of galaxies from same initial 
    M* bin for different prescriptions of SF evolution. More specifically
    this is to see how well the integrated SF stellar masses match the SHAM 
    stellar masses. 
    '''
    descendants = [] 
    for kwarg in kwarg_list:    # import SF inherited Lineages for different set of arguments
        sf_inherit_file = InheritSF_file(nsnap_descendant, abc_step=n_step, 
                **kwarg)
        bloodline = Lineage(nsnap_ancestor=kwarg['nsnap_ancestor'], subhalo_prop=kwarg['subhalo_prop'])
        bloodline.Read([nsnap_descendant], filename=sf_inherit_file)

        descendants.append(getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)))
    
    mf = SMF()
    
    mass_bins = np.arange(7., 12.5, 0.5)
    mass_bin_low = mass_bins[:-1]
    mass_bin_high = mass_bins[1:]

    prettyplot()
    pretty_colors = prettycolors()
    lstyles = ['-', '--', '-.']
    for i_m in range(len(mass_bin_low)): 
        fig = plt.figure(1)
        sub = fig.add_subplot(111)
        mass_evol_str = ''
        for i_desc, descendant in enumerate(descendants): 
            # SMF composition based on their initial mass at nsnap_ancestor 
            started_here = np.where(descendant.nsnap_genesis == kwarg_list[i_desc]['nsnap_ancestor'])
            start_mass = descendant.mass_genesis[started_here]

            mbin = np.where( 
                    (start_mass > mass_bin_low[i_m]) & 
                    (start_mass <= mass_bin_high[i_m])
                    )

            mass, phi = mf._smf(descendant.mass[started_here[0][mbin]])

            sub.plot(
                    mass, 
                    phi, 
                    c=pretty_colors[i_desc+1], 
                    lw=4, 
                    ls=lstyles[i_desc],
                    label=(kwarg_list[i_desc])['evol_prop']['mass']['name']+' Masses'
                    )
            mass_evol_str += (kwarg_list[i_desc])['evol_prop']['mass']['name']
        sub.set_xlim([7.5, 12.0])
        sub.set_ylim([0.0, 0.012])
        sub.set_ylabel(r'$\Phi$', fontsize=25)
        sub.set_xlabel(r'$\mathtt{log\;M_*}$', fontsize=25)
        sub.legend(loc='upper right')
        sub.set_title(''.join([
            'log M* = ', str(round(mass_bin_low[i_m],2)), '-', str(round(mass_bin_high[i_m],2)),
            ' at Snapshot ', str(kwarg_list[i_desc]['nsnap_ancestor'])]), fontsize=25)

        fig_file = ''.join([
            'figure/test/', 
            'DescendantSMF_composition', 
            '.massevol_', 
            mass_evol_str, 
            '.Mbin', str(round(mass_bin_low[i_m],2)), '_', str(round(mass_bin_high[i_m],2)),
            '.png'])
        fig.savefig(fig_file, bbox_inches='tight')
        plt.close()
    return None

def DescendantSFMS_composition(nsnap_descendant, abc_step=29, bovyplot=False, 
        **sfinh_kwargs): 
    ''' SFMS for the SF Inherited Descendant galaxy population  
    '''
    sfinherit_file = InheritSF_file(nsnap_descendant, abc_step=abc_step, **sfinh_kwargs)
    bloodline = Lineage(
            nsnap_ancestor=sfinh_kwargs['nsnap_ancestor'], 
            subhalo_prop=sfinh_kwargs['subhalo_prop']
            )
    bloodline.Read([nsnap_descendant], filename=sfinherit_file)

    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    if not bovyplot: 
        sfms_plot = descendant.plotSFMS(bovyplot=False, scatter=False) 

        # SMF composition based on their initial mass at nsnap_ancestor 
        started_here = np.where(descendant.nsnap_genesis == sfinh_kwargs['nsnap_ancestor'])
        start_mass = descendant.mass_genesis[started_here]
            
        mass_bins = np.arange(7., 12.5, 0.5)
        mass_bin_low = mass_bins[:-1]
        mass_bin_high = mass_bins[1:]

        for i_m in range(len(mass_bin_low)): 

            mbin = np.where( 
                    (start_mass > mass_bin_low[i_m]) & 
                    (start_mass <= mass_bin_high[i_m])
                    )

            sfms_plot.plot(
                    mass=descendant.mass[started_here[0][mbin]], 
                    sfr=descendant.sfr[started_here[0][mbin]], 
                    sfr_class=descendant.sfr_class[started_here[0][mbin]], 
                    gal_class='quiescent', 
                    bovyplot=False, 
                    sigSFR=False, 
                    color=i_m, 
                    label=r"$\mathtt{M_{*,i}=}$"+str(round(mass_bin_low[i_m],2))+"-"+str(round(mass_bin_high[i_m],2)))
        
        qfrac = Fq()
        m_arr = np.arange(9.0, 12.5, 0.5)
        sfms_plot.sub.plot(
                m_arr, qfrac.SFRcut(m_arr, descendant.zsnap, sfms_prop=(sfinh_kwargs['sfr_prop'])['sfms']), 
                c='k', ls='--', lw=4)
        bovyplot_str = ''
    else: 
        sfms_plot = descendant.plotSFMS(bovyplot=True) 
        bovyplot_str = '.bovy'

    sfms_plot.sub.legend(loc='lower right', scatterpoints=1)
    fig_file = ''.join([
        'figure/test/', 
        'SFMS_composition.', 
        ''.join((sfinherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        bovyplot_str, 
        '.png'])
    sfms_plot.save_fig(fig_file)
    return None

def DescendantSMHM_composition(nsnap_descendant, abc_step=29, **sfinh_kwargs):
    ''' Plot the Stellar Mass to Halo Mass relation. 
    '''
    sfinherit_file = InheritSF_file(nsnap_descendant, abc_step=abc_step, **sfinh_kwargs)
    bloodline = Lineage(
            nsnap_ancestor=sfinh_kwargs['nsnap_ancestor'], 
            subhalo_prop=sfinh_kwargs['subhalo_prop']
            )
    bloodline.Read([nsnap_descendant], filename=sfinherit_file)

    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    smhm_plot = descendant.plotSMHM(bovyplot=False, scatter=False)

    # SMF composition based on their initial mass at nsnap_ancestor 
    started_here = np.where(descendant.nsnap_genesis == sfinh_kwargs['nsnap_ancestor'])
    start_mass = descendant.mass_genesis[started_here]
        
    mass_bins = np.arange(7., 12.5, 0.5)
    mass_bin_low = mass_bins[:-1]
    mass_bin_high = mass_bins[1:]

    for i_m in range(len(mass_bin_low)): 

        mbin = np.where( 
                (start_mass > mass_bin_low[i_m]) & 
                (start_mass <= mass_bin_high[i_m])
                )

        smhm_plot.plot(
                stellarmass=descendant.mass[started_here[0][mbin]], 
                halomass=descendant.halo_mass[started_here[0][mbin]], 
                bovyplot=False, 
                color=i_m, 
                label=r"$\mathtt{M_{*,i}=}$"+str(round(mass_bin_low[i_m],2))+"-"+str(round(mass_bin_high[i_m],2)))

    smhm_plot.plotSummary(stellarmass=descendant.mass, halomass=descendant.halo_mass)

    smhm_plot.sub.legend(loc='upper left', scatterpoints=1)
    smhm_plot.set_axes()
    fig_file = ''.join([
        'figure/test/', 
        'SMHM_composition.', 
        ''.join((sfinherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.png'])
    smhm_plot.save_fig(fig_file)
    return None

def DescendantFQ(nsnap_descendant, abc_step=29, **sfinh_kwargs): 
    ''' FQ for the SF Inherited Descendant galaxy population  
    '''
    sfinherit_file = InheritSF_file(nsnap_descendant, abc_step=abc_step, **sfinh_kwargs)
    bloodline = Lineage(
            nsnap_ancestor=sfinh_kwargs['nsnap_ancestor'], 
            subhalo_prop=sfinh_kwargs['subhalo_prop']
            )
    bloodline.Read([nsnap_descendant], filename=sfinherit_file)

    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    descendant.sfr_prop = sfinh_kwargs['sfr_prop']
    
    fig_file = ''.join([
        'figure/test/', 
        'Fq.', 
        '.'.join((sfinherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.png'])
    descendant.plotFq(model=True, savefig=fig_file)
    return None

def Read_InheritSF(nsnap_descendant, nsnap_ancestor=20, n_step=29, 
    subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, 
    sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
    evol_prop = {'sfr': {'dutycycle': {'name': 'notperiodic'}}, 'mass': {'name': 'sham'}}): 
    
    sf_inherit_file = InheritSF_file(nsnap_descendant, nsnap_ancestor=nsnap_ancestor, 
        abc_step=n_step, subhalo_prop=subhalo_prop, sfr_prop=sfr_prop, evol_prop=evol_prop)

    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read([nsnap_descendant], filename=sf_inherit_file)

    return bloodline

def Save_InheritSF(nsnap_descendant, nsnap_ancestor=20, n_step=29, 
    subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, sfr_prop=None, 
    evol_prop=None): 
    '''
    '''
    med_theta = abc_posterior_median(n_step)
    pq_dict = {'slope': med_theta[0], 'yint': med_theta[1]}
    tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}
    if evol_prop['type'] == 'pq_based': 
        evol_prop['pq'] = pq_dict
    evol_prop['tau'] = tau_dict

    bloodline = InheritSF(
            nsnap_descendant, 
            nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop, 
            sfr_prop=sfr_prop, evol_prop=evol_prop)
    
    sfinherit_file = InheritSF_file(
                nsnap_descendant, 
                nsnap_ancestor=nsnap_ancestor, 
                abc_step=n_step, 
                subhalo_prop=subhalo_prop, 
                sfr_prop=sfr_prop, 
                evol_prop=evol_prop)
    print 'Writing ', sfinherit_file 
    bloodline.Write(sfinherit_file)
    
    return None

def InheritSF_file(nsnap_descendant, nsnap_ancestor=20, abc_step=None, 
    subhalo_prop={'scatter': 0.2, 'source': 'li-march'}, 
    sfr_prop={'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
    evol_prop={'sfr': {'dutycycle': {'name': 'notperiodic'}}, 'mass': {'name': 'sham'}}):
    '''
    '''
    abc_str = ''
    if abc_step is not None:
        abc_str = '.ABC'+str(abc_step)
    # SFMS string
    sfms_str = '.sfms_'+sfr_prop['sfms']['name']
    if sfr_prop['sfms']['name'] == 'linear': 
        sfms_str += ''.join([
            '_M', str(round(sfr_prop['sfms']['mslope'],2)), 
            'z', str(round(sfr_prop['sfms']['zslope'],2))])
    elif sfr_prop['sfms']['name'] == 'kinked': 
        sfms_str += ''.join([
            '_Mhigh', str(round(sfr_prop['sfms']['mslope_highmass'],2)), 
            '_Mlow', str(round(sfr_prop['sfms']['mslope_lowmass'],2)), 
            'z', str(round(sfr_prop['sfms']['zslope'],2))])
    # SFR assigned based on subhalo growth 
    if 'subhalogrowth' in sfr_prop.keys(): 
        sfms_str += '.halogrowth_SFRassign'
    # mass evolution string 
    if evol_prop['mass']['name'] == 'sham': 
        mass_str = '.Mevo_SHAM'
    elif evol_prop['mass']['name'] == 'integrated':
        if evol_prop['mass']['type'] == 'euler': 
            mass_str = '.Mevo_EulerInt'
        elif evol_prop['mass']['type'] == 'rk4': 
            mass_str = '.Mevo_RK4Int'

    hdf5_file = ''.join([
        'dat/InheritSF/'
        'InheritSF', 
        '.Snap', str(nsnap_descendant), '_', str(nsnap_ancestor), 
        '.subhalo_sig', str(round(subhalo_prop['scatter'], 2)), '_', subhalo_prop['source'], 
        '.fq_', sfr_prop['fq']['name'], 
        sfms_str, 
        '.', evol_prop['type'], '_evol',
        '.SF_Duty_', evol_prop['sfr']['dutycycle']['name'],  
        mass_str, 
        abc_str, '.hdf5'])
    return hdf5_file 

def abc_posterior_median(n_step): 
    '''
    Median theat from ABC posterior 
    '''
    # read in ABC posterior file 
    posterior_file = ''.join([
        'dat/pmc_abc/', 
        'theta_t', 
        str(n_step), 
        '.dat'])

    theta = np.loadtxt(
            posterior_file, 
            unpack = True
            ) 

    med_theta = [] 
    for i_param in xrange(len(theta)):
        med_theta.append( 
                np.median(theta[i_param])
                )

    return med_theta 

def kwargs_InheritSF(nsnap_ancestor=10, subhalo_prop=None, fq_prop=None, 
        sfms_evol=None, evol_type=None, dutycycle_evol=None, mass_evol=None, 
        subhalogrowth=None): 
    ''' Generate kwargs for Inherit SF
    '''
    # subhalo
    if subhalo_prop is None: 
        subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}
    # quiescent fraction
    if fq_prop is None: 
        fq_prop = {'name': 'wetzelsmooth'}
    # SFMS  
    if sfms_evol is None: 
        sfms_prop = {'name': 'linear', 'mslope': 0.55, 'zslope': 1.1}
    else: 
        if not isinstance(sfms_evol, dict):
            if sfms_evol == 'linear': 
                sfms_prop = {'name': 'linear', 'mslope': 0.55, 'zslope': 1.1}
            elif sfms_evol == 'kinked': 
                sfms_prop = {'name': 'kinked', 'mslope_highmass': 0.55, 'mslope_lowmass': 0.7, 'zslope': 1.1}
            else: 
                raise ValueError
        else: 
            sfms_prop = sfms_evol
    # dutycycle
    if evol_type is None: 
        evol_prop = 'pq_based'
    else: 
        evol_prop = 'simult'
    if dutycycle_evol is None: 
        dutycycle_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3}
    else: 
        if not isinstance(dutycycle_evol, dict):
            if dutycycle_evol == 'newamp_squarewave': 
                dutycycle_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'sigma': 0.3}
            elif dutycycle_evol == 'notperiodic': 
                dutycycle_prop = {'name': 'notperiodic'}
        else: 
            dutycycle_prop = dutycycle_evol
    if mass_evol is None: 
        raise ValueError
    else: 
        if mass_evol == 'euler':
            mass_prop = {'name': 'integrated', 'type':'euler', 'f_retain': 0.6, 't_step': 0.05}
        if mass_evol == 'rk4':
            mass_prop = {'name': 'integrated', 'type':'rk4', 'f_retain': 0.6, 't_step': 0.05}
        elif mass_evol == 'sham': 
            mass_prop = {'name': 'sham'} 
    kwargs = {
            'nsnap_ancestor': nsnap_ancestor, 
            'subhalo_prop': subhalo_prop, 
            'sfr_prop': {
                'fq': fq_prop, 
                'sfms': sfms_prop
                },
            'evol_prop': {
                'type': evol_type,
                'sfr': {'dutycycle': dutycycle_prop},
                'mass': mass_prop                 
                }
            }
    # Assign SFR based on subhalo growth 
    if subhalogrowth is not None: 
        if subhalogrowth:  
            kwargs['sfr_prop']['subhalogrowth'] = {} 
    return kwargs


"""
    def qaplot_sf_inherit_average_scatter(
            nsnap_descendants, nsnap_ancestor = 20, scatter = 0.0, n_step = 29,
            sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
            massevol_prop = {'name': 'sham'}
            ):
        '''
        Test Lineage SF Inherit function
        '''
        med_theta = abc_posterior_median(n_step)         # Get media theta from ABC posterior

        # quenching probabilyt properties
        pq_dict = {'slope': med_theta[0], 'yint': med_theta[1]}
        # tau properties
        tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}

        start_time = time.time()
        bloodline = sf_inherit(
                nsnap_descendants, 
                nsnap_ancestor = nsnap_ancestor, 
                ancestor_sf_prop = {'name': 'average'}, 
                pq_prop = pq_dict, 
                tau_prop = tau_dict,            
                sfrevol_prop = sfrevol_prop, 
                massevol_prop = massevol_prop, 
                quiet = True, 
                scatter = scatter
                )
        print 'SF inherit took ', (time.time() - start_time)/60.0, ' minutes' 
        
        for nsnap_descendant in nsnap_descendants: 

            descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))
            
            if massevol_prop['name'] == 'sham': 
                massevol_str = massevol_prop['name']
            else: 
                massevol_str = massevol_prop['type'] + massevol_prop['name']

            lineage_str = ''.join([ 
                '_', 
                str(nsnap_ancestor), 'ancestor_', 
                str(nsnap_descendant), 'descendant_', 
                str(round(scatter)), 'Mscatter_', 
                sfrevol_prop['name'], '_sfrevol_SFRMt_',
                massevol_str, '_massevol'
                ])

            for justsf in [True, False]: 
                justsf_str = ''
                if justsf: 
                    justsf_str = '_justsf'

                fig_name = ''.join([
                    'figure/', 
                    'qaplot_sf_inherit_sfms', lineage_str, justsf_str, '_twoslope_sfms.png'
                    ])

                sfms_plot = descendant.plotSFMS(justsf=justsf, bovyplot=True)
                sfms_plot.param_sfms()
                
                print fig_name
                sfms_plot.save_fig(fig_name)

            plt.close()

    def sf_inherited_lineage(
            n_step, 
            nsnap_ancestor = 20, 
            scatter = 0.0, 
            sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
            massevol_prop = {'name': 'sham'}
            ): 
        '''
        Write out SF Inherited Lineage to file for record keeping purposes and faster access for plotting 
        '''

        med_theta = abc_posterior_median(n_step)

        nsnap_range = np.arange(1, nsnap_ancestor)
        
        bloodline = sf_inherit(
                list(nsnap_range),
                nsnap_ancestor = nsnap_ancestor, 
                pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
                tau_prop = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}, 
                sfrevol_prop = sfrevol_prop, 
                massevol_prop = massevol_prop, 
                quiet = True, 
                scatter = scatter
                )

        if massevol_prop['name'] == 'sham': 
            massevol_str = ''.join([massevol_prop['name'], '_Mevol_'])
        else: 
            massevol_str = ''.join([
                massevol_prop['type'], massevol_prop['name'], '_Mevol_', 
                str(round(massevol_prop['t_step'], 3)), '_tstep_'
                ])

        lineage_str = ''.join([ 
            '_', 
            str(round(scatter)), 'Mscatter_', 
            sfrevol_prop['name'], '_sfrevol_SFRMt_', 
            massevol_str, str(n_step), '_abcpost'
            ])

        bloodline.file_name = ''.join([
                'dat/lineage/', 
                'lineage_ancestor_', 
                str(bloodline.nsnap_ancestor), 
                '_descendants',
                lineage_str, 
                '.hdf5'
                ]) 
        bloodline.writeout()

        return None

    def track_sf_pop_sfms_evol(
            n_gal, nsnap_ancestor = 20, n_step=29, scatter = 0.0, justsf=True, 
            sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
            massevol_prop = {'name': 'sham'}
            ):
        '''
        Directly track the evolution of n_gal galaxies in SFR versus M* paradigm. 
        Makes cool plots, 10 out of 10 would do again.
        '''
        if massevol_prop['name'] == 'sham': 
            massevol_str = ''.join([massevol_prop['name'], '_Mevol_'])
        else: 
            massevol_str = ''.join([
                massevol_prop['type'], massevol_prop['name'], '_Mevol_', 
                str(round(massevol_prop['t_step'], 3)), '_tstep_'
                ])

        lineage_str = ''.join([ 
            '_', 
            str(round(scatter)), 'Mscatter_', 
            sfrevol_prop['name'], '_sfrevol_SFRMt_', 
            massevol_str, str(n_step), '_abcpost'
            ])

        nsnap_range = np.arange(1, nsnap_ancestor)
        bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
        bloodline.file_name = ''.join([
                'dat/lineage/', 
                'lineage_ancestor_', 
                str(bloodline.nsnap_ancestor), 
                '_descendants',
                lineage_str, 
                '.hdf5'
                ]) 
        bloodline.readin(nsnap_range)

        for i_snap in nsnap_range: 
            descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(i_snap))
            
            if i_snap == nsnap_range[0]: 
                if justsf: 
                    gal_assigned = np.where(descendant.gal_type != 'quiescent')
                else: 
                    gal_assigned = np.where(descendant.gal_type != '')

            if i_snap == nsnap_range[0]: 
                sf_masses = descendant.mass[gal_assigned]
                sf_sfrs = descendant.sfr[gal_assigned]
            else: 
                sf_masses = np.vstack((sf_masses, descendant.mass[gal_assigned]))
                sf_sfrs = np.vstack((sf_sfrs, descendant.sfr[gal_assigned]))
        
        prettyplot()
        pretty_colors = prettycolors()
        fig = plt.figure(figsize=(10,10))
        sub = fig.add_subplot(111)
        
        for i in xrange(n_gal): #sf_masses.shape[1]):
            sub.plot(
                    sf_masses[:,i], 
                    sf_sfrs[:,i],
                    color=pretty_colors[i % 20],
                    lw=2
                    )
            sub.scatter(
                    sf_masses[:,i], 
                    sf_sfrs[:,i],
                    color=pretty_colors[i % 20],
                    s=16
                    )
        sub.set_xlim([9.0, 12.0])
        sub.set_ylim([-5.0, 2.0])
        sub.set_xlabel(r'$\mathtt{M_*}$')
        sub.set_ylabel(r'$\mathtt{log\;SFR}$')
        
        justsf_str = ''
        if justsf: 
            justsf_str = '_justsf'
        sfms_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_sfms', lineage_str, justsf_str, '_galaxytrack.png'
            ])
        fig.savefig(sfms_fig_file, bbox_inches="tight")
        plt.close()
        
        fig = plt.figure(figsize=(15,15))
        fig.subplots_adjust(hspace=0., wspace=0.)
        
        for i_attr, attr in enumerate(['sfr', 'mass', 'ssfr']): 
            sub = fig.add_subplot(3, 1, i_attr+1) 

            for i in xrange(n_gal): #sf_masses.shape[1]):
                
                if attr == 'sfr': 
                    attr_val = sf_sfrs[:,i]
                elif attr == 'mass': 
                    attr_val = sf_masses[:,i]
                elif attr == 'ssfr': 
                    attr_val = sf_sfrs[:,i] - sf_masses[:,i]

                sub.plot( get_t_nsnap(nsnap_range), 
                        attr_val, 
                        color=pretty_colors[i % 20],
                        lw=2
                        )

            sub.set_xlim([3, 14.0])
            sub.set_ylabel(r'$\mathtt{log\;'+attr.upper()+'}$')

            if i_attr == 2:
                sub.set_ylim([-11.0, -8.0])
                sub.set_xlabel(r'$\mathtt{t_{cosmic}}$', fontsize=30)
            elif i_attr == 1: 
                sub.set_ylim([7.5, 11.0])
                sub.set_xticklabels([])
            elif i_attr == 0: 
                sub.set_ylim([-2.0, 1.5])
                sub.set_xticklabels([])
        
        attr_fig_file = ''.join([
            'figure/', 
            'qaplot_sf_inherit_ssfr', lineage_str, justsf_str, '_galaxytrack.png'
            ])
        fig.savefig(attr_fig_file, bbox_inches="tight")
        #plt.show()
        plt.close()

    def qaplot_sf_inherited_lineage(
            nsnap_range=None, nsnap_ancestor = 20, scatter = 0.0, n_step = 29, 
            sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
            massevol_prop = {'name': 'sham'}, 
            ssfr=True, fq=False, tau=False, sfms=False, mass_scatter=False, smf=False
            ): 
        '''
        '''
        if massevol_prop['name'] == 'sham': 
            massevol_str = ''.join([massevol_prop['name'], '_Mevol_'])
        else: 
            massevol_str = ''.join([
                massevol_prop['type'], massevol_prop['name'], '_Mevol_', 
                str(round(massevol_prop['t_step'], 3)), '_tstep_'
                ])

        lineage_str = ''.join([ 
            '_', 
            str(round(scatter)), 'Mscatter_', 
            sfrevol_prop['name'], '_sfrevol_SFRMt_', 
            massevol_str, str(n_step), '_abcpost'
            ])
        
        if nsnap_range is None: 
            nsnap_range = np.arange(1, nsnap_ancestor)

        bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
        bloodline.file_name = ''.join([
                'dat/lineage/', 
                'lineage_ancestor_', 
                str(bloodline.nsnap_ancestor), 
                '_descendants',
                lineage_str, 
                '.hdf5'
                ]) 
        bloodline.readin(nsnap_range)
        
        attr_list = []
        if ssfr: 
            attr_list.append('ssfr')
        if fq: 
            attr_list.append('fq')
        if tau: 
            attr_list.append('tau')
        if sfms: 
            attr_list.append('sfms')
        if smf: 
            attr_list.append('smf')
        if mass_scatter: 
            attr_list.append('mass_scatter')
        
        for i_snap in nsnap_range: 
            descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(i_snap))

            for attr in attr_list: 
                fig_name = ''.join([
                    'figure/', 
                    'qaplot_sf_inherit_lineage_', attr, lineage_str, '_', str(i_snap), 'of', str(nsnap_ancestor),'.png'
                    ])

                if attr == 'ssfr':                          # SSFR
                    descendant.plotSsfr(
                            line_color='red', 
                            line_width=4, 
                            label=r'Median $\mathtt{\tau_{Q}}$ of Posterior', 
                            groupcat = True, 
                            savefig = fig_name                    
                            )
                elif attr == 'fq':                          # Quiescent Fraction 
                    descendant.plotFq(
                            param=True, 
                            line_color='r', 
                            label = None,
                            savefig= fig_name                    
                            )
                elif attr == 'tau':                         # Quenching Timescale
                    descendant.plotTau(savefig = fig_name)
                elif attr == 'sfms':                        # Star Forming Main Sequence 
                    descendant.plotSFMS(justsf=True, bovyplot=True, savefig = fig_name)
                elif attr == 'smf': 
                    descendant.plotSMF(savefig = fig_name)
                elif attr == 'mass_scatter':                # Stellar Mass - Halo Mass 
                    descendant.plotMstarMhalo(savefig = fig_name)
        return None


        def qaplot_sf_inherit_nosfr_scatter(
                n_step, 
                nsnap_descendants,
                nsnap_ancestor = 20, 
                scatter = 0.0, 
                sfrevol_prop = {'name': 'squarewave', 'freq_range': [0.0, 2*np.pi], 'phase_range': [0, 1]},
                sfrevol_massdep = False,
                massevol_prop = {'name': 'sham'}
                ):
            '''
            Test Lineage SF Inherit function for the case where there's no scatter in 
            the SFMS
            '''
            
            if sfrevol_massdep: 
                sfrevol_massdep_str = 'SFRMt_'
            else: 
                sfrevol_massdep_str = 'SFRM0t_'

            posterior_file = ''.join([
                'dat/pmc_abc/', 
                'theta_t', 
                str(n_step), 
                '.dat'])

            theta = np.loadtxt(
                    posterior_file, 
                    unpack = True
                    ) 

            med_theta = [] 
            for i_param in xrange(len(theta)):
                med_theta.append( 
                        np.median(theta[i_param])
                        )

            bloodline = sf_inherit(
                    nsnap_descendants, 
                    nsnap_ancestor = nsnap_ancestor, 
                    ancestor_sf_prop = {'name': 'average_noscatter'}, 
                    pq_prop = {'slope': med_theta[0], 'yint': med_theta[1]}, 
                    tau_prop = {
                        'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]
                        }, 
                    sfrevol_prop = sfrevol_prop, 
                    sfrevol_massdep = sfrevol_massdep,
                    massevol_prop = massevol_prop, 
                    quiet = True, 
                    scatter = scatter
                    )
            
            prettyplot()
            pretty_colors = prettycolors()
            
            sfr_mstar_z, sig_sfr_mstar_z = get_param_sfr_mstar_z()

            for nsnap_descendant in nsnap_descendants: 
                descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(nsnap_descendant))

                fig = plt.figure()
                sub = fig.add_subplot(111)

                sf_gal = np.where(descendant.gal_type == 'star-forming')

                sub.scatter(
                        descendant.mass[sf_gal], 
                        descendant.sfr[sf_gal], 
                        color=pretty_colors[nsnap_descendant],
                        s=10
                        )
                sub.plot(np.arange(8.5, 12.0, 0.1), 
                        sfr_mstar_z(np.arange(8.5, 12.0, 0.1), get_z_nsnap(nsnap_descendant)), 
                        c = 'k', 
                        ls = '--', 
                        lw = 3 
                        )
                sub.set_xlim([9.0, 12.0])
                sub.set_ylim([-5.0, 2.0])
                sub.set_xlabel(r'$\mathtt{M_*}$')
                sub.set_ylabel(r'$\mathtt{log\;SFR}$')
            
                file_str = ''.join([ 
                    '_', 
                    str(nsnap_ancestor), 'ancestor_', 
                    str(nsnap_descendant), 'descendant_', 
                    str(round(scatter)), 'Mscatter_', 
                    sfrevol_prop['name'], '_sfrevol_',
                    sfrevol_massdep_str, 
                    massevol_prop['name'], '_massevol'
                    ])

                sfms_fig_file = ''.join([
                    'figure/', 
                    'qaplot_sf_inherit_sfms', file_str, '_no_sfr_scatter_twoslope_sfms.png'
                    ])
                fig.savefig(sfms_fig_file, bbox_inches='tight')
                plt.close()

            for scat in [0.2]: 
                start_time = time.time()
                bloodline = Lineage(nsnap_ancestor = 20)
                bloodline.descend(subhalo_prop = {'scatter': scat, 'source': 'li-march'}, clobber=True) 
                bloodline.assign_sfr_ancestor(sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}})
                bloodline.writeout()
                print 'lineage construction and write out takes ', (time.time() - start_time)/60.0

                for id in [1, 5, 9, 15, 19]:#, 3, 5, 7, 9, 11, 13, 15, 17, 19]:
                    #qaplot_sf_inherit(
                    #    nsnap_ancestor = 20, nsnap_descendant = id, 
                    #    subhalo_prop = {'scatter': scat, 'source': 'li-march'}, 
                    #    sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
                    #    evol_prop = {
                    #        'sfr': {'name': 'newamp_squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
                    #        'mass': {'name': 'sham'} 
                    #        },
                    #    ssfr=True, fq=True, tau=False, mass_scatter=True, sfms=True, smf=True
                    #    )
                    qaplot_sf_inherit(
                        nsnap_ancestor = 20, nsnap_descendant = id, 
                        subhalo_prop = {'scatter': scat, 'source': 'li-march'}, 
                        sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
                        evol_prop = {
                            'sfr': {'name': 'newamp_squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
                            'mass': {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.05} 
                            },
                        ssfr=True, fq=True, tau=False, mass_scatter=True, sfms=True, smf=True
                        )

            # {'name': 'newamp_squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},

            # {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.05}
            #'name': 'squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1]
            #qaplot_sf_inherit_average_scatter(
            #        [1],
            #        sfrevol_prop = 
            #        massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.01}
            #        )   # {'name': 'notperiodic'}

            #for integ in ['euler', 'rk4']:
            #start_time = time.time()
            #sf_inherited_lineage(
            #        29, 
            #        nsnap_ancestor = 20, 
            #        scatter = 0.0, 
            #        sfrevol_prop = {
            #            'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3
            #            },
            #        massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.0025}
            #        )
            #print (time.time() - start_time)/60.0, ' minutes'
            #for tstep in [0.01, 0.0025]: 
            #    for integ in ['euler', 'rk4']: 
            #        qaplot_sf_inherited_lineage(
            #                nsnap_range=[1],
            #                sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3}, 
            #                massevol_prop = {'name': 'integrated', 'type': integ, 'f_retain': 0.6, 't_step': tstep}, 
            #                ssfr=False, fq=False, tau=False, sfms=False, mass_scatter=False, smf=True
            #                )
            #track_sf_pop_sfms_evol(
            #        10,
            #        sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
            #        massevol_prop = {'name': 'integrated', 'type': 'rk4', 'f_retain': 0.6, 't_step': 0.0025}, 
            #        justsf=False
            #        )
"""


if __name__=="__main__":
    #print abc_posterior_median(29)
    mevo = 'euler'
    for duty in ['notperiodic']:#, 'newamp_squarewave']: 
        if duty == 'notperiodic':  
            shgrow = True
        else: 
            shgrow = False
        for esemeffes in ['kinked']: #['linear']: # , 'kinked']:
            kwargs = kwargs_InheritSF(
                    nsnap_ancestor=20,
                    sfms_evol=esemeffes,
                    dutycycle_evol=duty, 
                    mass_evol=mevo, 
                    evol_type='simult',
                    subhalogrowth=shgrow)
            for nsnap in [1]: #range(1,20)[::-1]: 
                Save_InheritSF(nsnap, **kwargs)
                #DescendantSMHM_composition(nsnap, **kwargs)
                DescendantSFMS_composition(nsnap, **kwargs)
                DescendantFQ(nsnap, **kwargs)
                #DescendantSMF_composition(nsnap, **kwargs)
            AnecstorPlots(n_descendant=9, **kwargs)
