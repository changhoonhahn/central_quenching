'''

Plots to test SF-MS module 


'''
import bovy_plot as bovy
from util.gal_classify import sfr_cut
from sfms.sf_mainseq import get_sfr_mstar_z_groupcat
from sfms.sf_mainseq import get_sfr_mstar_z_envcount
from sfms.fitting import get_param_sfr_mstar_z
from sfms.fitting import get_bestfit_sfms_groupcat
from sfms.fitting import get_bestfit_sfms_envcount
from group_catalog.group_catalog import sf_centrals
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors
from defutility.fitstables import mrdfits

def qaplot_sfms_groupcat_fitting(Mrcut=18, sfq_test=True):
    """ Function tests the SF-MS fitting functions and examines the 
    SF-MS of group catalog altogether as a sanity check. 

    If sfq_test is specified, then SFR(M*,z) cutoff for SF/Q classification
    is plotted on top of the SF-MS plots
    """

    # star-forming central galaxies from the SDSS group catalog
    sf_cen = sf_centrals(Mrcut=Mrcut) 

    #prettyplot()        
    pretty_colors = prettycolors() 

    bovy.scatterplot(
            sf_cen.mass,
            sf_cen.sfr,
            scatter = True,
            levels = [0.68, 0.95, 0.997],
            color = pretty_colors[1],
            s = 3,
            xrange = [9.5, 12.0],
            yrange = [-1.5, 1.5],
            xlabel = 'log \; M_{*}',
            ylabel = 'log \; SFR'
            )

    # SDSS group catalog best fit 
    gc_zmid, gc_slope, gc_yint = get_bestfit_sfms_groupcat(Mrcut=Mrcut, clobber=True)

    mass_bin = np.arange(9.0, 12.0, 0.25)       # mass bin 

    plt.plot(
            mass_bin, 
            gc_slope * (mass_bin-10.5) + gc_yint, 
            c='k', 
            lw=6, 
            ls='--'
            ) 

    # Envcount catalog best fit 
    ec_zmid, ec_slope, ec_yint = get_bestfit_sfms_envcount()

    for i_z in xrange(len(ec_zmid)):
        plt.plot(
                mass_bin, 
                ec_slope[i_z] * (mass_bin-10.5) + ec_yint[i_z], 
                c = pretty_colors[i_z+3], 
                lw = 4, 
                label = str(ec_zmid[i_z])
                ) 

    avg_sfrs, sig_sfrs, ngals = get_sfr_mstar_z_groupcat(
            mass_bin,
            Mrcut=Mrcut
            ) 
    enough_gal = np.where(np.array(avg_sfrs) > -998)

    plt.errorbar(
            mass_bin[enough_gal], 
            np.array(avg_sfrs)[enough_gal], 
            yerr = np.array(sig_sfrs)[enough_gal], 
            lw = 4, 
            c = pretty_colors[1], 
            label='Average SFR'
            )

    if sfq_test: 

        plt.plot(mass_bin, sfr_cut(mass_bin, 0.05), c='k', lw=4, ls='--')
    
    plt.legend(loc='lower right')
    
    fig_name = ''.join(['figure/', 'qaplot_sfms_fitting_groupcat_', str(Mrcut), '.png'])
    plt.savefig(fig_name, bbox_inches='tight')
    plt.close()
    
def qaplot_sfms_envcount_fitting(Mrcut=18, sfq_test=True):
    """ Test functions of the SF-MS module that fits the group catalog SFMS
    
    If sfq_test is specified, then SFR(M*,z) cutoff for SF/Q classification
    is plotted on top of the SF-MS plots
    """
    #prettyplot()        
    pretty_colors = prettycolors() 

    # read SF galaxies from the envcount catalog  
    file_dir = 'dat/wetzel_tree/envcount/'
    sdss_envcount_file = ''.join([file_dir, 
        'envcount_cylr2.5h35_thresh75_sdss_active_z0.05_0.12_primuszerr.fits']) 
    sdss_sf_data = mrdfits(sdss_envcount_file) 
    primus_envcount_file = ''.join([file_dir, 
        'envcount_cylr2.5h35_thresh75_active_z0.2_1.0_lit.fits']) 
    primus_sf_data = mrdfits(primus_envcount_file) 
    
    for i_z, z_mid in enumerate([0.1, 0.3, 0.5, 0.7, 0.9]): 

        if z_mid < 0.2: 
            sf_data = sdss_sf_data

            # impose isolation criteria, mass completeness limit (should I ?) and edge cuts
            centrals = np.where(
                    (sf_data.envcount == 0.0) & 
                    (sf_data.mass > sf_data.masslimit) & 
                    (sf_data.edgecut == 1) 
                    )
        else: 
            sf_data = primus_sf_data
            # impose isolation criteria, mass completeness limit (should I ?) and edge cuts
            centrals = np.where(
                    (sf_data.redshift >= z_mid - 0.1) &
                    (sf_data.redshift < z_mid + 0.1) &
                    (sf_data.envcount == 0.0) & 
                    (sf_data.mass > sf_data.masslimit) & 
                    (sf_data.edgecut == 1) 
                    )

        print sf_data.weight[centrals]

        bovy.scatterplot(
                sf_data.mass[centrals],
                sf_data.sfr[centrals],
                scatter = True,
                levels = [0.68, 0.95, 0.997],
                weights = sf_data.weight[centrals], 
                s = 3,
                xrange = [9.5, 12.0],
                yrange = [-1.5, 1.5],
                xlabel = 'log \; M_{*}',
                ylabel = 'log \; SFR'
                )
                #color = pretty_colors[1],

        # SDSS group catalog best fit 
        gc_zmid, gc_slope, gc_yint = get_bestfit_sfms_groupcat(Mrcut=Mrcut, clobber=True)

        mass_bin = np.arange(9.0, 12.0, 0.25)       # mass bin 

        plt.plot(
                mass_bin, 
                gc_slope * (mass_bin-10.5) + gc_yint, 
                c='k', 
                lw=6, 
                ls='--'
                ) 

        # Envcount catalog best fit 
        ec_zmid, ec_slope, ec_yint = get_bestfit_sfms_envcount(clobber=True)

        plt.plot(
                mass_bin, 
                ec_slope[i_z] * (mass_bin-10.5) + ec_yint[i_z], 
                c = pretty_colors[i_z+3], 
                lw = 4, 
                label = str(ec_zmid[i_z])
                ) 

        avg_sfrs, sig_sfrs, ngals = get_sfr_mstar_z_envcount(
            mass_bin,
            [z_mid for i in xrange(len(mass_bin))]
            ) 
        enough_gal = np.where(np.array(avg_sfrs) > -998)

        plt.errorbar(
                mass_bin[enough_gal], 
                np.array(avg_sfrs)[enough_gal], 
                yerr = np.array(sig_sfrs)[enough_gal], 
                lw = 4, 
                c = pretty_colors[1], 
                label='Average SFR'
                )

        plt.legend(loc='lower right')

        if sfq_test: 

            plt.plot(mass_bin, sfr_cut(mass_bin, z_mid), c='k', lw=4, ls='--')
        
        fig_name = ''.join(['figure/', 'qaplot_sfms_fitting_envcount_z', str(z_mid), '.png'])
        plt.savefig(fig_name, bbox_inches='tight')
        plt.close()

    return None

def qaplot_parameterized_sfms(): 
    """ QAplot comparison for the parameterized SFMS redshift 
    evolution in comparison to observed SF-MS redshift evolution.

    """
    gc_zmid, gc_slope, gc_yint = get_bestfit_sfms_groupcat()
    
    ec_zmid, ec_slope, ec_yint = get_bestfit_sfms_envcount(
            fid_mass = 10.5
            )

    avg_sfr_sfms, sig_sfr_sfms = get_param_sfr_mstar_z()

    mass_bin = np.arange(9.0, 12.0, 0.25)       # mass bin 

    fig = plt.figure()
    sub = fig.add_subplot(111)
    
    sub.scatter(
            0.03, 
            gc_yint,
            c='k', 
            s=50, 
            label='Observed GroupCatalog'
            )

    sub.scatter(
            ec_zmid, 
            ec_yint, 
            c=prettycolors()[3], 
            label='Observed EnvCount'
            )
        
    sfr_fidmass = np.array([
        avg_sfr_sfms(10.5, zmid) 
        for zmid in np.arange(0.0, 1.0, 0.05)
        ])
    sub.plot(
            np.arange(0.0, 1.0, 0.05),
            sfr_fidmass, 
            c='k', 
            lw=4, 
            ls='--',
            label='Parameterized'
            )
    
    sub.plot(
            np.arange(0.0, 1.0, 0.05),
            sfr_fidmass+0.15, 
            c=prettycolors()[3], 
            lw=4, 
            ls='--'
            )

    sfr_cutoff = np.array([
        sfr_cut(10.5, ec_zmid[i_z]) 
        for i_z in xrange(len(ec_zmid))
        ])
    sub.plot(ec_zmid, sfr_cutoff, c='k', lw=3, ls='--', label='SF/Q Classification')

    sub.set_xlim([0.0, 1.0])
    sub.set_ylim([-1.0, 2.0])
    sub.set_ylabel('SFR(M=10.5,z)', fontsize=20)
    sub.set_xlabel('Redshift (z)', fontsize=20)
    sub.legend(scatterpoints=1, loc='upper left')

    plt.show()

if __name__=="__main__":
    #qaplot_sfms_groupcat_fitting(Mrcut=18)
    #qaplot_sfms_envcount_fitting(Mrcut=18)
    qaplot_parameterized_sfms()
