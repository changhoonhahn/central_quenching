'''

Plots to test SF-MS module 


'''
import bovy_plot as bovy
from sfms.sf_mainseq import get_sfr_mstar_z_groupcat
from sfms.sf_mainseq import get_sfr_mstar_z_envcount
from sfms.fitting import get_bestfit_sfms_groupcat
from sfms.fitting import get_bestfit_sfms_envcount
from group_catalog.group_catalog import sf_centrals
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors
from defutility.fitstables import mrdfits

def qaplot_sfms_groupcat_fitting(Mrcut=18):
    """ Function tests the SF-MS fitting functions and examines the 
    SF-MS of group catalog altogether as a sanity check. 
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
    
    plt.legend(loc='lower right')
    
    fig_name = ''.join(['figure/', 'qaplot_sfms_fitting_groupcat_', str(Mrcut), '.png'])
    plt.savefig(fig_name, bbox_inches='tight')
    plt.close()
    
def qaplot_sfms_envcount_fitting(Mrcut=18):
    """ Test functions of the SF-MS module that fits the group catalog SFMS
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
        
        fig_name = ''.join(['figure/', 'qaplot_sfms_fitting_envcount_z', str(z_mid), '.png'])
        plt.savefig(fig_name, bbox_inches='tight')
        plt.close()

    return None

if __name__=="__main__":
    qaplot_sfms_groupcat_fitting(Mrcut=18)
    qaplot_sfms_envcount_fitting(Mrcut=18)
