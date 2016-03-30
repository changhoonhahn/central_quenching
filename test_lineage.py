'''

Module to test Lineage class object. More specifically, 
to test the ancestor/descendant system 

'''
import numpy as np 
import matplotlib.pyplot as plt

from gal_prop import SMF
from plotting import plots
from lineage import Lineage
from central_subhalo import Subhalos
from central_subhalo import CentralSubhalos

from defutility.plotting import quick_hist
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 


def LineageSMF(nsnap_ancestor, descendants=None, 
        subhalo_prop={'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop={'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}): 
    ''' Plot the SMF of lineage galaxy population (both ancestor and descendants).
    Compare the lineage SMF to the central subhalo and analytic SMFs for all 
    subhalos. The main agreement, is between the lineage SMF and the central subhalo
    SMF. If nsnap_ancestor is high, then there should be more discrepancy 
    between LIneage SMF and CentralSubhalos SMF because more galaxies are lost. 
    '''
    # read in desendants from the lineage object
    if descendants is None: 
        descendants = range(1, nsnap_ancestor)
    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read(descendants, sfr_prop=sfr_prop) 
    
    smf_plot = plots.PlotSMF()
    # plot ancestor against the analytic SMF 
    ancestor = getattr(bloodline, 'ancestor')
    smf_plot.cenque(ancestor)   # SMF of lineage ancestor
    smf_plot.analytic(
            ancestor.zsnap, 
            source=ancestor.subhalo_prop['source'], 
            line_style='-', lw=1,
            label='All Subhalos')

    smf = SMF()
    subh = CentralSubhalos()
    subh.Read(
            nsnap_ancestor, 
            scatter=ancestor.subhalo_prop['scatter'], 
            source=ancestor.subhalo_prop['source'])
    subh_mass, subh_phi = smf.centralsubhalos(subh)
    smf_plot.sub.plot(subh_mass, subh_phi, lw=4, ls='--', c='gray', label='Central Subhalos')

    for isnap in descendants: 
        descendant = getattr(bloodline, 'descendant_snapshot'+str(isnap))
        smf = SMF()
        subh = CentralSubhalos()
        subh.Read(
                isnap, 
                scatter=descendant.subhalo_prop['scatter'], 
                source=descendant.subhalo_prop['source'])
        d_mass, d_phi = smf_plot.cenque(descendant)
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        smf_plot.sub.plot(subh_mass, subh_phi, lw=4, ls='--', c='gray')
        #print (subh_phi - d_phi)/d_phi

        smf_plot.analytic(
                descendant.zsnap, 
                source=descendant.subhalo_prop['source'], 
                line_style='-', lw=1,
                label=None)

    smf_plot.set_axes()
    smf_plot_file = ''.join([
        'figure/test/', 
        'LineageSMF', 
        '.ancestor', str(nsnap_ancestor), 
        bloodline._file_spec(subhalo_prop=bloodline.subhalo_prop, sfr_prop=bloodline.sfr_prop), 
        '.png'
        ])
    smf_plot.save_fig(smf_plot_file)
    return None

def LineageAncestorSFMS(nsnap_ancestor, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
        ): 
    '''Plot the SF-MS of the lineage ancestor object. This is mainly to make sure 
    that the AssignSFR routine is working properly 
    '''
    # read in lineage
    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read([1], sfr_prop=sfr_prop)
    ancestor = getattr(bloodline, 'ancestor')

    sfms_plot = plots.PlotSFMS()
    sfms_plot.cenque(ancestor, justsf=True)   # lineage ancestor
    sfms_plot.param_sfms(nsnap=nsnap_ancestor) 

    sfms_plot_file = ''.join([
        'figure/test/', 
        'LineageAncestorSFMS', str(nsnap_ancestor), 
        bloodline._file_spec(subhalo_prop=subhalo_prop, sfr_prop=sfr_prop),
        '.png'])
    sfms_plot.save_fig(sfms_plot_file)
    return None

def LineageFinalDescendantSMF(nsnap_ancestor, 
        subhalo_prop={'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop={'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}): 
    ''' Plot SMF of final descendant. Also mark the composition of the final SMF from Subhalos 
    that gain stellar mass at different snapshots in the middle of the simulation. 
    '''
    prettyplot()
    pretty_colors=prettycolors()

    bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
    bloodline.Read([1], sfr_prop=sfr_prop)

    final_desc = bloodline.descendant_snapshot1
    
    mf = SMF()
    
    smf_plot = plots.PlotSMF()
    smf_plot.cgpop(final_desc, line_color='k')
    
    for isnap in range(2, nsnap_ancestor+1)[::-1]:
        started_here = np.where(final_desc.nsnap_genesis == isnap)
        mass, phi = mf._smf(final_desc.mass[started_here])

        try: 
            phi_tot += phi 
            smf_plot.sub.fill_between(mass, phi_tot - phi, phi_tot, color=pretty_colors[isnap], label='Nsnap='+str(isnap))
        except UnboundLocalError: 
            phi_tot = phi 
            smf_plot.sub.fill_between(mass, np.zeros(len(phi)), phi, color=pretty_colors[isnap], label='Nsnap='+str(isnap))

    smf_plot.set_axes() 
    fig_file = ''.join([
        'figure/test/',
        'LineageFinalDescendantSMF', 
        '.ancestor', str(nsnap_ancestor),
        bloodline._file_spec(subhalo_prop=subhalo_prop, sfr_prop=sfr_prop), 
        '.png'])
    smf_plot.save_fig(fig_file)

    return None

"""
    def ancestoral_past(nsnap_ancestor, 
            subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
            sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
            ): 
        ''' Check the ancestoral past of descendants that supposedly have ancestors
        '''
        # read in lineage class 
        bloodline = Lineage(nsnap_ancestor=nsnap_ancestor, subhalo_prop=subhalo_prop)
        bloodline.Read(range(1, nsnap_ancestor), sfr_prop=sfr_prop) 

        ancestor = getattr(bloodline, 'ancestor')
        descendant = getattr(bloodline, 'descendant_snapshot1')
        ancestor_index = getattr(descendant, 'ancestor'+str(nsnap_ancestor))
        
        noancestor = ancestor_index[np.invert(np.in1d(ancestor_index, ancestor.snap_index))]
        print 'No Ancestor Ngal=', len(noancestor)
        
        centsub = CentralSubhalos()
        centsub.Read(nsnap_ancestor, scatter=subhalo_prop['scatter'], source=subhalo_prop['source'])
        
        prettyplot()
        pretty_colors = prettycolors()
        fig = plt.figure()
        sub = fig.add_subplot(111)
        x, y = quick_hist(centsub.mass[noancestor], normed=True) 
        sub.plot(x, y, c='k', lw=2, label = str(nsnap_ancestor))
        for isnap in range(1, nsnap_ancestor)[::-1]: 
            centsub = Subhalos()
            centsub.Read(isnap, scatter=subhalo_prop['scatter'], source=subhalo_prop['source'])

            ancestor_index_i = getattr(centsub, 'ancestor'+str(nsnap_ancestor))
            no_i = np.in1d(ancestor_index_i, noancestor)
            print len(centsub.mass[no_i])
            x, y = quick_hist(centsub.mass[no_i], normed=True)
            sub.scatter(x, y, c=pretty_colors[isnap], lw=2, label = str(isnap))
        sub.set_xlim([0,12])
        sub.set_yscale("log")
        sub.legend()
        plt.show()
"""



if __name__=='__main__': 
    #LineageSMF(10, descendants = [1, 5],
    #        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
    #        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
    #        )
    #LineageFinalDescendantSMF(20, 
    #        subhalo_prop={'scatter': 0.0, 'source': 'li-march'}, 
    #        sfr_prop={'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}})
    #LineageAncestorSFMS(20, 
    #        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
    #        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}})
