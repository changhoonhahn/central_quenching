'''

Plots to test SF-MS module 


'''
import bovy_plot as bovy
from sfms.fitting import get_bestfit_sfms_groupcat
from sfms.fitting import get_bestfit_sfms_envcount
from group_catalog.group_catalog import sf_centrals
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors

def qaplot_sfms_groupcat_fitting(Mrcut=18):
    """ Test functions of the SF-MS module that fits the group catalog SFMS
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
    gc_zmid, gc_slope, gc_yint = get_bestfit_sfms_groupcat(Mrcut=Mrcut)

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
    plt.legend()
    
    fig_name = ''.join(['figure/', 'qaplot_sfms_fitting_groupcat_', str(Mrcut), '.png'])
    plt.savefig(fig_name, bbox_inches='tight')

    #mass, avg_sfrs, var_sfrs = [], [], [] 
    #for i_mass in range(len(mass_bin.mass_low)): 
    #    avg_sfr, var_sfr, ngal = sfms.get_sfr_mstar_z_groupcat(mass_bin.mass_mid[i_mass], 
    #            Mrcut=Mrcut) 

    #    if ngal < 10: 
    #        continue    # skip mass bin if there aren't many galaxies
    #    if mass_bin.mass_mid[i_mass] < 9.5: 
    #        continue    # skip low mass

    #    mass.append(mass_bin.mass_mid[i_mass]) 
    #    avg_sfrs.append(avg_sfr)
    #    var_sfrs.append(var_sfr)
    #    
    #subs.errorbar(mass, avg_sfrs, yerr=var_sfrs, 
    #        lw=4, c=pretty_colors[1], label='Average SFR')

    #bestfit_params = sfms.get_bestfit_groupcat_sfms(Mrcut=Mrcut, clobber=True) 
    #subs.plot(mass, util.line(np.array(mass)-10.5, bestfit_params), 
    #        lw=4, ls='--', c=pretty_colors[2], label='MPfit line') 
    #subs.text(11.25, 0.0, 'Slope = '+str("%.2f" % bestfit_params[0]))
    #subs.text(11.25, -0.5, 'Y-int = '+str("%.2f" % bestfit_params[1]))

    #subs.legend(loc='upper right') 
    

if __name__=="__main__":
    qaplot_sfms_groupcat_fitting(Mrcut=18)
