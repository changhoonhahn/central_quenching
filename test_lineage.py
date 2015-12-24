
from smf import SMF
from plotting import plots
from lineage import Lineage
from central_subhalo import Subhalos
from central_subhalo import CentralSubhalos

from defutility.plotting import quick_hist
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 


def lineage_smf(nsnap_ancestor, descendants = None, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
        ): 
    '''
    Test SMF of lineage 
    '''
    # read in lineage
    if descendants is None: 
        descendants = range(1, nsnap_ancestor)
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.readin(range(1, nsnap_ancestor), subhalo_prop = subhalo_prop, sfr_prop = sfr_prop) 
    
    smf_plot = plots.PlotSMF()
    # plot ancestor 
    ancestor = getattr(bloodline, 'ancestor_cq')
    smf_plot.cenque(ancestor)   # lineage ancestor
    smf_plot.analytic(
            ancestor.zsnap, 
            source=ancestor.subhalo_prop['source'], 
            label='All Subhalos')
    smf = SMF()
    subh = CentralSubhalos()
    subh.read(
            nsnap_ancestor, 
            scatter=ancestor.subhalo_prop['scatter'], 
            source=ancestor.subhalo_prop['source'])
    subh_mass, subh_phi = smf.centralsubhalos(subh)
    smf_plot.sub.plot(subh_mass, subh_phi, lw=4, ls='--', c='gray', label='Central Subhalos')

    for isnap in descendants: 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(isnap))
        smf = SMF()
        subh = CentralSubhalos()
        subh.read(
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
                label=None)

    smf_plot.set_axes()
    smf_plot_file = ''.join(['figure/', 
                'qaplot_lineage_ancestor', str(nsnap_ancestor), 
                descendant._file_spec(subhalo_prop = descendant.subhalo_prop, sfr_prop = descendant.sfr_prop), 
                'smf.png'
                ])
    smf_plot.save_fig(smf_plot_file)

def lineage_ancestor_sfms(nsnap_ancestor, descendants = None, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
        ): 
    '''
    Test SFMS of lineage ancestor_cq
    '''
    # read in lineage
    if descendants is None: 
        descendants = range(1, nsnap_ancestor)
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.readin(range(1, nsnap_ancestor), subhalo_prop = subhalo_prop, sfr_prop = sfr_prop) 
    
    ancestor = getattr(bloodline, 'ancestor_cq')
    sfms_plot = plots.PlotSFMS()
    sfms_plot.cenque(ancestor, justsf=True)   # lineage ancestor
    descendant = getattr(bloodline, 'descendant_cq_snapshot1')
    sfms_plot_file = ''.join(['figure/', 
                'qaplot_lineage_ancestor', str(nsnap_ancestor), 
                descendant._file_spec(subhalo_prop = descendant.subhalo_prop, sfr_prop = descendant.sfr_prop), 
                'sfms.png'
                ])
    sfms_plot.save_fig(sfms_plot_file)

def ancestoral_past(nsnap_ancestor, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
        ): 
    '''
    Check the ancestoral past of descendants that supposedly have ancestors
    '''
    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.ancestor(subhalo_prop = subhalo_prop, sfr_prop = sfr_prop)
    bloodline.readin(range(1, nsnap_ancestor)) 
    ancestor = getattr(bloodline, 'ancestor_cq')
    descendant = getattr(bloodline, 'descendant_cq_snapshot1')
    ancestor_index = getattr(descendant, 'ancestor'+str(nsnap_ancestor))
    
    noancestor = ancestor_index[np.invert(np.in1d(ancestor_index, ancestor.snap_index))]
    print 'No Ancestor Ngal=', len(noancestor)
    
    centsub = CentralSubhalos()
    centsub.read(nsnap_ancestor, scatter=subhalo_prop['scatter'], source=subhalo_prop['source'])
    
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    x, y = quick_hist(centsub.mass[noancestor], normed=True) 
    sub.plot(x, y, c='k', lw=2, label = str(nsnap_ancestor))
    for isnap in range(1, nsnap_ancestor)[::-1]: 
        centsub = Subhalos()
        centsub.read(isnap, scatter=subhalo_prop['scatter'], source=subhalo_prop['source'])

        ancestor_index_i = getattr(centsub, 'ancestor'+str(nsnap_ancestor))
        no_i = np.in1d(ancestor_index_i, noancestor)
        print len(centsub.mass[no_i])
        x, y = quick_hist(centsub.mass[no_i], normed=True)
        sub.scatter(x, y, c=pretty_colors[isnap], lw=2, label = str(isnap))
    sub.set_xlim([0,12])
    sub.set_yscale("log")
    sub.legend()
    plt.show()



if __name__=='__main__': 
    #lineage_smf(20, 
    #        descendants = [1, 10],
    #        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
    #        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
    #        )
    lineage_ancestor_sfms(20, 
            descendants = [1, 10],
            subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
            sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
            )

    #ancestoral_past(
    #        20, 
    #        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
    #        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
    #        )
