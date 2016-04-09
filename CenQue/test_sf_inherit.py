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
from observations import GroupCat
from central_subhalo import CentralSubhalos

from plotting.plots import PlotSMF

# --- plotting --- 
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors

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

def AncestorPlots(**sfinh_kwargs): 
    ''' FQ for the SF Inherited Ancestor galaxy population. This test is mainly 
    to see what the discrepancy between the empirical SF/Q classifcation and the
    parameterized FQ model. Ideally, the model and the 
    '''
    bloodline = Lineage(
            nsnap_ancestor=sfinh_kwargs['nsnap_ancestor'], 
            subhalo_prop=sfinh_kwargs['subhalo_prop']
            )
    bloodline.Read([1])
    bloodline.AssignSFR_ancestor(sfr_prop=sfinh_kwargs['sfr_prop'])

    ancestor = getattr(bloodline, 'ancestor') 
    ancestor.sfr_prop = sfinh_kwargs['sfr_prop']
    fig_file = lambda prop: ''.join(['figure/test/', 
        'Ancestor', str(sfinh_kwargs['nsnap_ancestor']), 
        '.', prop.upper(), '.png'])
    ancestor.plotSsfr(savefig=fig_file('ssfr')) 
    plt.close()
    # FQ
    ancestor.plotFq(model=True, savefig=fig_file('fq'))      
    ancestor.plotSFMS(sfqcut=True, gal_class='all', bovyplot=False, sigSFR=False, model=False, savefig=fig_file('sfms'))
    plt.close()
    return None

def DescendantQAplot(nsnap_descendant, **sfinh_kwargs): 
    ''' The ultimate QAplot t rule them all. 4 panels showing all the properties.
    '''
    sfinherit_file = InheritSF_file(nsnap_descendant, **sfinh_kwargs)
    bloodline = Lineage(
            nsnap_ancestor=sfinh_kwargs['nsnap_ancestor'], 
            subhalo_prop=sfinh_kwargs['subhalo_prop']
            )
    bloodline.Read([nsnap_descendant], filename=sfinherit_file)

    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    descendant.sfr_prop = sfinh_kwargs['sfr_prop']
    started_here = np.where(descendant.nsnap_genesis == sfinh_kwargs['nsnap_ancestor'])
    start_mass = descendant.mass_genesis[started_here]

    # Mass bins
    mass_bins = np.arange(7., 12.5, 0.5)
    mass_bin_low = mass_bins[:-1]
    mass_bin_high = mass_bins[1:]

    plt.close()
    prettyplot()
    fig = plt.figure(1, figsize=[25,6])
    for i_sub in range(1,5): 
        sub_i = fig.add_subplot(1,4,i_sub) 
        
        if i_sub == 1:  # SMF
            mf = SMF()
            mass, phi = mf.Obj(descendant)

            sub_i.plot(mass, phi, lw=4, c=pretty_colors[descendant.nsnap], label=r'Simulated')

            censub = CentralSubhalos()
            censub.Read(descendant.nsnap, 
                    scatter=sfinh_kwargs['subhalo_prop']['scatter'], 
                    source=sfinh_kwargs['subhalo_prop']['source'], 
                    nsnap_ancestor=sfinh_kwargs['nsnap_ancestor'])

            mass, phi = mf._smf(censub.mass)
            sub_i.plot(mass, phi, c='k', lw=4, ls='--', label='Central Subhalos') 

            sub_i.set_ylim([10**-5, 10**-1])
            sub_i.set_xlim([7.5, 12.0])
            plt.xticks([8., 9., 10., 11., 12.])
            sub_i.set_yscale('log')

            # x,y labels
            sub_i.set_xlabel(r'Mass $\mathtt{M_*}$', fontsize=25) 
            sub_i.set_ylabel(r'Stellar Mass Function $\mathtt{\Phi}$', fontsize=25) 
            sub_i.legend(loc='upper right', frameon=False)

        elif i_sub == 2: # SFMS
            # SMF composition based on their initial mass at nsnap_ancestor 
            for i_m in range(len(mass_bin_low)): 
                mbin = np.where((start_mass > mass_bin_low[i_m]) & (start_mass <= mass_bin_high[i_m]))

                sub_i.scatter(
                        descendant.mass[started_here[0][mbin]], 
                        descendant.sfr[started_here[0][mbin]], 
                        color=pretty_colors[i_m], 
                        label=r"$\mathtt{M_{*,i}=}$"+str(round(mass_bin_low[i_m],2))+"-"+str(round(mass_bin_high[i_m],2))
                        )
            
            qfrac = Fq()
            m_arr = np.arange(9.0, 12.5, 0.5)
            sub_i.plot(
                    m_arr, 
                    qfrac.SFRcut(m_arr, descendant.zsnap, sfms_prop=(sfinh_kwargs['sfr_prop'])['sfms']), 
                    c='k', ls='--', lw=4)

            sub_i.set_xlim([9.0, 12.0])
            sub_i.set_ylim([-5.0, 2.0])
            sub_i.set_xlabel(r'$\mathtt{log\;M_*}$', fontsize=25)
            sub_i.set_ylabel(r'$\mathtt{log\;SFR}$', fontsize=25)

        elif i_sub == 3: #SMHM
            for i_m in range(len(mass_bin_low)): 
                mbin = np.where((start_mass > mass_bin_low[i_m]) & (start_mass <= mass_bin_high[i_m]))

                sub_i.scatter(
                        descendant.halo_mass[started_here[0][mbin]], 
                        descendant.mass[started_here[0][mbin]], 
                        color=pretty_colors[i_m], 
                        label=r"$\mathtt{M_{*,i}=}$"+str(round(mass_bin_low[i_m],2))+"-"+str(round(mass_bin_high[i_m],2))
                        )

            stellarmass = descendant.mass[started_here]
            halomass = descendant.halo_mass[started_here]
            mbin = np.arange(halomass.min(), halomass.max(), 0.25) 
            mlow = mbin[:-1]
            mhigh = mbin[1:]
            
            muMstar = np.zeros(len(mlow))
            sigMstar = np.zeros(len(mlow))

            for im in range(len(mlow)): 
                mbin = np.where((halomass > mlow[im]) & (halomass <= mhigh[im]))
                muMstar[im] = np.mean(stellarmass[mbin])
                sigMstar[im] = np.std(stellarmass[mbin])
            sub_i.errorbar(0.5*(mlow+mhigh), muMstar, yerr=sigMstar, color='k', lw=3, fmt='o', capthick=2)

            sub_i.set_ylim([9.0, 12.0])
            sub_i.set_xlim([10.0, 15.0])
            sub_i.set_ylabel(r'Stellar Mass $\mathtt{M_*}$', fontsize=25) 
            sub_i.set_xlabel(r'Halo Mass $\mathtt{M_{Halo}}$', fontsize=25) 

            #sub_i.legend(loc='upper left', frameon=False, scatterpoints=1)
        elif i_sub == 4: # Fq
            #mass, fq = descendant.Fq()
            sfq = qfrac.Classify(descendant.mass, descendant.sfr, descendant.zsnap, 
                    sfms_prop=descendant.sfr_prop['sfms'])
            gc = GroupCat(Mrcut=18, position='central')
            gc.Read()
            gc_sfq = qfrac.Classify(gc.mass, gc.sfr, np.mean(gc.z), 
                    sfms_prop=descendant.sfr_prop['sfms'])
            #sub_i.plot(mass, fq, color=pretty_colors[descendant.nsnap], lw=3, ls='--', 
            #        label=r'$\mathtt{z =} '+str(descendant.zsnap)+'$')
            M_bin = np.array([9.7, 10.1, 10.5, 10.9, 11.3])
            M_low = M_bin[:-1]
            M_high = M_bin[1:]
            M_mid = 0.5 * (M_low + M_high)

            fq = np.zeros(len(M_low))
            gc_fq = np.zeros(len(M_low))
            for i_m in xrange(len(M_low)):
                mlim = np.where((descendant.mass > M_low[i_m]) & (descendant.mass <= M_high[i_m]))
                gc_mlim = np.where((gc.mass > M_low[i_m]) & (gc.mass <= M_high[i_m]))
                ngal = np.float(len(mlim[0]))
                gc_ngal = np.float(len(gc_mlim[0]))

                if ngal != 0:  # no galaxy in mass bin 
                    ngal_q = np.float(len(np.where(sfq[mlim] == 'quiescent')[0]))
                    fq[i_m] = ngal_q/ngal
                if gc_ngal != 0:
                    gc_ngal_q = np.float(len(np.where(gc_sfq[gc_mlim] == 'quiescent')[0]))
                    gc_fq[i_m] = gc_ngal_q/gc_ngal

            sub_i.plot(M_mid, fq, color=pretty_colors[descendant.nsnap], lw=3, label=r'$\mathtt{z =} '+str(descendant.zsnap)+'$')
            
            fq_model = qfrac.model(M_bin, descendant.zsnap, lit=sfinh_kwargs['sfr_prop']['fq']['name'])
            sub_i.plot(M_bin, fq_model, color='k', lw=4, ls='--', label=sfinh_kwargs['sfr_prop']['fq']['name'])
            sub_i.scatter(M_mid, gc_fq, color='k', s=100, lw=0, label='Group Catalog')

            sub_i.set_xlim([9.0, 12.0])
            sub_i.set_ylim([0.0, 1.0])

            sub_i.set_xlabel(r'Mass $\mathtt{M_*}$') 
            sub_i.set_ylabel(r'Quiescent Fraction $\mathtt{f_Q}$', fontsize=20) 

            sub_i.legend(loc='upper left', frameon=False, scatterpoints=1, markerscale=0.75)

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0) 
    fig_name = ''.join([
        'figure/test/', 
        'QAplot.',
        '.'.join((sfinherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.png'])
    fig.savefig(fig_name, bbox_inches='tight') 
    plt.close()
    return None 

def DescendantSSFR(nsnap_descendant, **kwargs):
    ''' Halo Mass Function composition at z = z_final by initial halo mass at 
    nsnap_ancestor. 
    '''
    # import SF inherited Lineage
    sf_inherit_file = InheritSF_file(nsnap_descendant, **kwargs)
    bloodline = Lineage(nsnap_ancestor=kwargs['nsnap_ancestor'], subhalo_prop=kwargs['subhalo_prop'])
    bloodline.Read([nsnap_descendant], filename=sf_inherit_file)
    # descendant
    descendant = getattr(bloodline, 'descendant_snapshot'+str(nsnap_descendant)) 
    fig_name = ''.join([
        'figure/test/', 
        'SSFR.',
        '.'.join((sf_inherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
        '.png'])

    descendant.plotSsfr(line_color='red', line_width=4, 
            sfms_prop=kwargs['sfr_prop']['sfms'], z=descendant.zsnap, 
            groupcat=True, savefig=fig_name)
    plt.close()
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
        '.'.join((sfinherit_file.rsplit('/')[-1]).rsplit('.')[:-1]), 
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
    descendant.plotFq(model=sfinh_kwargs['sfr_prop']['fq']['name'], savefig=fig_file)
    return None

def Read_InheritSF(nsnap_descendant, nsnap_ancestor=20, n_step=29, 
    subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}, 
    sfr_prop = {'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
    evol_prop = {'sfr': {'dutycycle': {'name': 'notperiodic'}}, 'mass': {'name': 'sham'}}, 
    flag=None): 
    
    sf_inherit_file = InheritSF_file(nsnap_descendant, nsnap_ancestor=nsnap_ancestor, 
        abc_step=n_step, subhalo_prop=subhalo_prop, sfr_prop=sfr_prop, evol_prop=evol_prop, 
        flag=flag)

    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read([nsnap_descendant], filename=sf_inherit_file)

    return bloodline

def Save_InheritSF(nsnap_descendant, **kwargs): 
    '''
    '''
    bloodline = InheritSF(
            nsnap_descendant, 
            nsnap_ancestor=kwargs['nsnap_ancestor'], 
            subhalo_prop=kwargs['subhalo_prop'], 
            sfr_prop=kwargs['sfr_prop'], 
            evol_prop=kwargs['evol_prop'])
    
    sfinherit_file = InheritSF_file(nsnap_descendant, **kwargs)
    print 'Writing ', sfinherit_file 
    bloodline.Write(sfinherit_file)
    
    return None

def InheritSF_file(nsnap_descendant, nsnap_ancestor=20, tau='abc', abc_step=None, 
    subhalo_prop={'scatter': 0.2, 'source': 'li-march'}, 
    sfr_prop={'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}},
    evol_prop={'sfr': {'dutycycle': {'name': 'notperiodic'}}, 'mass': {'name': 'sham'}}, 
    flag=None):
    '''
    '''
    tau_str = ''
    if abc_step is not None:
        tau_str = '.tau_ABC'+str(abc_step)
    else: 
        tau_str = evol_prop['tau']['name']

    # SFMS string
    sfms_str = '.SFMS'+sfr_prop['sfms']['name'].lower()
    if sfr_prop['sfms']['name'] == 'linear': 
        sfms_str += ''.join([
            '_z', str(round(sfr_prop['sfms']['zslope'],2))])
    elif sfr_prop['sfms']['name'] == 'kinked': 
        sfms_str += ''.join([
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
    flag_str = ''
    if flag is not None: 
        flag_str = '.'+flag

    hdf5_file = ''.join([
        code_dir(), 'dat/InheritSF/'
        'InheritSF', 
        '.Snap', str(nsnap_descendant), '_', str(nsnap_ancestor), 
        '.subhalo_sig', str(round(subhalo_prop['scatter'], 2)), '_', subhalo_prop['source'], 
        '.fq_', sfr_prop['fq']['name'], 
        sfms_str, 
        '.', evol_prop['type'], '_evol',
        '.SF_Duty_', evol_prop['sfr']['dutycycle']['name'],  
        mass_str, 
        tau_str, 
        flag_str, '.hdf5'])
    return hdf5_file 

def abc_posterior_median(n_step): 
    '''
    Median theat from ABC posterior 
    '''
    # read in ABC posterior file 
    posterior_file = ''.join([
        code_dir, 'dat/pmc_abc/', 
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

def kwargs_InheritSF(nsnap_ancestor=10, tau='abc', subhalo_prop=None, fq_prop=None, gv_prop=None, fudge=None, 
        sfms_evol=None, evol_type=None, dutycycle_evol=None, mass_evol=None, 
        subhalogrowth=None): 
    ''' Generate kwargs for Inherit SF
    '''
    # subhalo
    if subhalo_prop is None: 
        subhalo_prop = {'scatter': 0.2, 'source': 'li-march'}
    # quiescent fraction
    if fq_prop is None: 
        fq_prop = {'name': 'wetzel'}
    if gv_prop is None: 
        gv_prop = {'slope': 0.0, 'fidmass':10.5, 'offset': 0.0} 
    # SFMS  
    if sfms_evol is None: 
        sfms_prop = {'name': 'linear', 'zslope': 1.5}
    else: 
        if not isinstance(sfms_evol, dict):
            if sfms_evol == 'linear': 
                sfms_prop = {'name': 'linear', 'zslope': 1.5}
            elif sfms_evol == 'kinked': 
                sfms_prop = {'name': 'kinked', 'mslope_lowmass': 0.7, 'zslope': 1.5}
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
            mass_prop = {'name': 'integrated', 'type':'euler', 'f_retain': 0.6, 't_step': 0.1}
        if mass_evol == 'rk4':
            mass_prop = {'name': 'integrated', 'type':'rk4', 'f_retain': 0.6, 't_step': 0.1}
        elif mass_evol == 'sham': 
            mass_prop = {'name': 'sham'} 

    evol_dict = { 
            'type': evol_type,
            'sfr': {'dutycycle': dutycycle_prop},
            'mass': mass_prop                 
            }
    if tau == 'abc': 
        med_theta = abc_posterior_median(29)
        tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': med_theta[2], 'yint': med_theta[3]}
    elif tau == 'satellite':
        tau_dict={'name': 'satellite'} 
    elif isinstance(tau, list): 
        tau_dict = {'name': 'line', 'fid_mass': 11.1, 'slope': tau[0], 'yint': tau[1]}
    evol_dict['tau'] = tau_dict 
    if fudge is None:
        evol_dict['fudge'] = {'slope': 0.0, 'fidmass': 10.5, 'offset': 1.0}
    else: 
        evol_dict['fudge'] = fudge 


    kwargs = {
            'nsnap_ancestor': nsnap_ancestor, 
            'subhalo_prop': subhalo_prop, 
            'sfr_prop': {
                'fq': fq_prop, 
                'sfms': sfms_prop, 
                'gv': gv_prop
                },
            'evol_prop': evol_dict
            }
    if tau == 'abc': 
        kwargs['abc_step'] = 29
    # Assign SFR based on subhalo growth 
    if subhalogrowth is not None: 
        if subhalogrowth:  
            kwargs['sfr_prop']['subhalogrowth'] = {} 
    return kwargs



if __name__=="__main__":
    print abc_posterior_median(29)
    mevo = 'sham'
    for nsnap_a in [15]: 
        for duty in ['newamp_squarewave']: 
            if duty == 'notperiodic':  
                shgrow = True
            else: 
                shgrow = False
            for esemeffes in ['kinked']: #'linear', 'kinked']:
                #        tau='satellite', 
                #        fudge={'slope': -1.5, 'fidmass': 10.5, 'offset': 1.} ,
                #        tau='abc', 
                #        fudge={'slope': -1.25, 'fidmass': 10.5, 'offset': 1.75} ,
                kwargs = kwargs_InheritSF(
                        nsnap_ancestor=nsnap_a,
                        sfms_evol=esemeffes,
                        tau='abc', 
                        fudge={'slope': -1.25, 'fidmass': 10.5, 'offset': 1.75} ,
                        gv_prop={'slope': 0.3, 'fidmass': 10.5, 'offset': 0.3},
                        dutycycle_evol=duty, 
                        mass_evol=mevo, 
                        evol_type='simult',
                        subhalogrowth=shgrow)
                flag = 'testtest'
                #AncestorPlots(**kwargs)
                for nsnap in [1]: #range(1,nsnap_a)[::-1]: 
                    Save_InheritSF(nsnap, flag=flag, **kwargs)
                    DescendantQAplot(nsnap, flag=flag, **kwargs)
                    DescendantSSFR(nsnap, flag=flag, **kwargs)
