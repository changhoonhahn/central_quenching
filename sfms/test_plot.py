'''
'''
import os 
import h5py
import numpy as np
import scipy as sp 

# --- Local ---
from defutility.plotting import prettyplot
from defutility.plotting import prettycolors
from sfms.fitting import get_param_sfr_mstar_z

def plot_sfms_wetzel_vs_lee():
    '''
    Plot Wetzel parameterized SFMS and Lee et al. (2015) parameterized SFMS. 
    '''
    
    # wetzel 
    wetzel_sfr_mstar_z, wetzel_sig_sfr_mstar_z = get_param_sfr_mstar_z()

    # Lee et al. (2015)

    z_mid = np.array([0.36, 0.55, 0.70, 0.85, 0.99, 1.19])
    logM_lim = np.array([8.50, 9.00, 9.00, 9.30, 9.30, 9.30])

    masses = np.arange(7.5, 12.0, 0.1) 

    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for i_z, z in enumerate(z_mid): 

        sub.scatter(
                masses, 
                lee_et_al_2015(masses, z),
                c = pretty_colors[i_z],
                label = r"$\mathtt{z = "+str(round(z, 2))+"}$"
                )

        sub.plot(
                masses, 
                wetzel_sfr_mstar_z(masses, z),
                c = pretty_colors[i_z], 
                ls = '--', 
                lw = 3 
                )
        #sub.vlines(logM_lim[i_z], -5.0, 2.0, color='k', linestyle='--', lw=2)

    sub.set_xlim([7.0, 12.0])
    sub.set_ylim([-3.0, 2.0])
    sub.set_xlabel(r'$\mathtt{log\;M_*}$', fontsize=30)
    sub.set_ylabel(r'$\mathtt{log\;SFR}$', fontsize=30)
    sub.legend(loc='lower right', scatterpoints=1)
    fig_file = ''.join([
        'figure/', 
        'qaplot_sfms_wetzel_vs_lee.png'
        ])
    fig.savefig(fig_file, bbox_inches='tight')

    plt.show()

def lee_et_al_2015(M, z): 
    '''
    Lee et al. (2015) SFMS parameterization
    '''

    z_min = np.array([0.25, 0.46, 0.63, 0.78, 0.93, 1.11])
    z_max = np.array([0.46, 0.63, 0.78, 0.93, 1.11, 1.30])

    z_mid = np.array([0.36, 0.55, 0.70, 0.85, 0.99, 1.19])

    logM_lim = np.array([8.50, 9.00, 9.00, 9.30, 9.30, 9.30])

    S0 = np.array([0.80, 0.99, 1.23, 1.35, 1.53, 1.72])
    S0_err = np.array([0.019, 0.015, 0.016, 0.014, 0.017, 0.024])

    logM0 = np.array([10.03, 9.82, 9.93, 9.96, 10.10, 10.31])
    logM0_err = np.array([0.042, 0.031, 0.031, 0.025, 0.029, 0.043])

    gamma = np.array([0.92, 1.13, 1.11, 1.28, 1.26, 1.07])
    gamma_err = np.array([0.017, 0.033, 0.025, 0.034, 0.032, 0.028])

    i_z = np.where( (z_min < z) & (z_max >= z) ) 

    if len(i_z[0]) != 1: 
        raise ValueError
 
    S = S0[i_z] - np.log10( 1. + (10.**M / 10.**logM0[i_z])**(-1.0 * gamma[i_z]))
    
    return S

if __name__=="__main__":
    plot_sfms_wetzel_vs_lee()
