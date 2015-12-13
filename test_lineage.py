
from smf import SMF
from plotting import plots
from lineage import Lineage
from central_subhalo import CentralSubhalos


def lineage_smf(nsnap_ancestor, descendants = None, 
        subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
        sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
        ): 
    '''
    Test SMF of lineage 
    '''
    if descendants is None: 
        descendants = range(1, nsnap_ancestor)

    bloodline = Lineage(nsnap_ancestor = nsnap_ancestor)
    bloodline.ancestor(subhalo_prop = subhalo_prop, sfr_prop = sfr_prop)
    bloodline.readin(range(1, nsnap_ancestor)) 

    smf_plot = plots.PlotSMF()
    ancestor = getattr(bloodline, 'ancestor_cq')
    smf_plot.cenque(ancestor)
    smf_plot.analytic(ancestor.zsnap, source=ancestor.subhalo_prop['source'], label=None)

    smf = SMF()
    subh = CentralSubhalos()
    subh.read(nsnap_ancestor, scatter=ancestor.subhalo_prop['scatter'], source=ancestor.subhalo_prop['source'])
    subh_mass, subh_phi = smf.centralsubhalos(subh)
    smf_plot.sub.plot(subh_mass, subh_phi, lw=4, c='gray')

    for isnap in descendants: 
        descendant = getattr(bloodline, 'descendant_cq_snapshot'+str(isnap))
        #descendant.plotSMF(show = True)
        smf = SMF()
        subh = CentralSubhalos()
        subh.read(isnap, scatter=descendant.subhalo_prop['scatter'], source=descendant.subhalo_prop['source'])
        subh_mass, subh_phi = smf.centralsubhalos(subh)
        smf_plot.sub.plot(subh_mass, subh_phi, lw=4, c='gray')

        smf_plot.cenque(descendant)
        smf_plot.analytic(descendant.zsnap, source=descendant.subhalo_prop['source'], label=None)
            

    smf_plot.set_axes()
    smf_plot.save_fig(
            ''.join([
                'figure/', 
                'qaplot_lineage_ancestor', str(nsnap_ancestor), 
                descendant._file_spec(subhalo_prop = descendant.subhalo_prop, sfr_prop = descendant.sfr_prop), 
                'smf.png'
                ]))




if __name__=='__main__': 
    lineage_smf(20, 
            descendants = [19],
            subhalo_prop = {'scatter': 0.0, 'source': 'li-march'}, 
            sfr_prop = { 'fq': {'name': 'wetzelsmooth'}, 'sfr': {'name': 'average'}}
            )
