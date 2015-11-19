'''

Test integrated mass evolution 

'''
import numpy as np
from scipy import interpolate

import sfr_evol
import mass_evol

from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 
from sfms.fitting import get_param_sfr_mstar_z

def test_mass_evol(mass0, t0, tf, 
        delmass0 = 0.1, 
        t_initial = 4.0079, 
        sfrevol_param = [1.0, 1.0],
        sfrevol_prop = {'name': 'notperiodic'}, 
        massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.025}, 
        sfr_kwargs = 'm0'
        ): 
    '''

    Testing integrated mass evolution

    logSFR(z) function is implemented within this module

    '''

    # spline between z and t_cosmic
    z, t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
    z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 
    z0 = z_of_t(t0)
    zf = z_of_t(tf)
    
    # average SFR of SF MS
    logsfr_mstar_z, sig_logsfr_mstar_z = get_param_sfr_mstar_z()
            
    def logsfr_notquenched(logmass, t_input, **kwargs): 

        zi = kwargs['z0']
        mass_i = kwargs['mass0']
        delmass_i = kwargs['deltamass0']

        if kwargs['sfr_kwargs'] == 'm0': 
            avglogsfr = logsfr_mstar_z(mass_i, zi)
        elif kwargs['sfr_kwargs'] == 'm': 
            avglogsfr = logsfr_mstar_z(logmass, zi)
        #print 'SFMS(', mass0[0], ',', str(zi), ') ', avglogsfr[0]
    
        # Contribution from the SF evolution of the SFMS 
        logsfr_sfms = sfr_evol.logsfr_sfms_evol(zi, z_of_t(t_input))
        
        # Contribution from the SF Duty cycle 
        logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                t_initial, 
                t_input, 
                delta_sfr=delmass_i,
                sfrevol_param=sfrevol_param, 
                **sfrevol_prop
                )
        #print logsfr_sfduty

        logsfr_f = avglogsfr + logsfr_sfms + logsfr_sfduty
        #print np.max(logsfr_f - logsfr_mstar_z(logmass, z_of_t(t_input))), np.min(logsfr_f - logsfr_mstar_z(logmass, z_of_t(t_input)))

        return logsfr_f
    
    sfrkwargs = {'z0': z0, 'mass0': mass0, 'deltamass0': delmass0, 'sfr_kwargs': sfr_kwargs}
    mass_f, sfr_f = mass_evol.integrated_rk4(
            logsfr_notquenched, 
            mass0, 
            t0, 
            tf, 
            f_retain = massevol_prop['f_retain'], 
            delt = massevol_prop['t_step'], 
            **sfrkwargs
            )
    #print sfr_f[:3] 
    #print logsfr_notquenched(mass_f, tf, **sfrkwargs)[:3]
    #print logsfr_mstar_z(mass0, zf)[:3]

    return mass_f, sfr_f

def plot_integrated_mass_evol(mass0, t0, tf, 
        sfrevol_prop = {'name': 'notperiodic'},
        massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.025}, 
        sfr_kwargs = 'm0',
        title=None, fig_file=None):
    '''
    
    Plot of the integrated mass and SFR evolution as a function of cosmic time 

    '''
    if title is None: 
        raise ValueError 
    delt = (tf-t0)/10.0
    tfs = np.arange(t0+delt, tf+delt, delt)

    z, t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
    z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 
    
    # SF MS 
    logsfr_mstar_z, sig_logsfr_mstar_z = get_param_sfr_mstar_z()

    # figure
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(10,12))
    # mass evolution 
    m_sub = fig.add_subplot(211)
    # SFR evolution 
    sub = fig.add_subplot(212)

    delta_masses = np.arange(-0.3, 0.3, 0.05) 
    #delta_masses = [0.0]    # test case 
    for delmass in delta_masses: 
        # SFR evolution parameters (for SF duty cycle)
        sfrevol_param = sfr_evol.get_sfrevol_param(1, np.array([0]), **sfrevol_prop)

        masses = [mass0]
        sfmssfrs = [logsfr_mstar_z(mass0, z_of_t(t0))]
        sfrs = [logsfr_mstar_z(mass0, z_of_t(t0)) + delmass]

        mass_i_1 = mass0
        t_i_1 = t0

        for tf_i in tfs: 

            mass_i, sfr_i = test_mass_evol(
                    mass0, 
                    t0, 
                    tf_i, 
                    delmass0 = delmass,
                    t_initial = t0,
                    sfrevol_param=sfrevol_param, 
                    sfrevol_prop = sfrevol_prop,
                    sfr_kwargs = sfr_kwargs,
                    massevol_prop=massevol_prop
                    )

            masses.append(mass_i)
            sfrs.append(sfr_i)
            sfmssfrs.append(logsfr_mstar_z(mass0, z_of_t(tf_i)))
            
            t_i_1 = tf_i
            mass_i_1 = mass_i
        
        tfs = np.hstack([t0, tfs])
        for i in xrange(len(masses[0])):
            if delmass == delta_masses[0]: 
                m_label = str(mass0[i])
            else: 
                m_label = None
            m_sub.plot(tfs, np.array(masses)[:,i], c=pretty_colors[i], lw=3, alpha=0.5, label=m_label)
        
        for i in xrange(len(masses[0])):
            if delmass == delta_masses[0]: 
                sfr_label = str(mass0[i])
            else: 
                sfr_label = None
            sub.plot(tfs, np.array(sfrs)[:,i], c=pretty_colors[i], lw=2, alpha=0.5, label=sfr_label)
            sub.plot(tfs, np.array(sfmssfrs)[:,i], c='k', lw=2, ls='--')
        
    m_sub.set_xlim([t0-2.5, tf+0.5])
    m_sub.set_ylim([8.5, 12.5])
    m_sub.set_ylabel(r"$\mathtt{log\;M_*}$", fontsize=20)
    m_sub.set_title(title)

    sub.set_xlim([t0-2.5, tf+0.5])
    sub.set_ylim([-1.0, 2.5])
    sub.set_xlabel(r"$\mathtt{t_{cosmic}}$", fontsize=20)
    sub.set_ylabel(r"$\mathtt{log\;SFR}$", fontsize=20)
    sub.legend(loc='upper left')

    if fig_file is not None: 
        figfile = ''.join([
            'figure/', 
            fig_file
            ])
        fig.savefig(figfile, bbox_inches='tight')
    plt.show()
    plt.close()

if __name__=="__main__": 
    mass0 = np.arange(7.5, 12.5, 0.5)
    '''
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; SFR(M,t); Constant SF Duty Cycle",  
            sfrevol_prop = {'name': 'notperiodic'},
            massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.05},
            sfr_kwargs = 'm',
            fig_file = 'test_massevol_sfr_m_t_notperiodic_dutycycle.png'
            )
    
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; SFR(M0,t); Constant SF Duty Cycle",  
            sfrevol_prop = {'name': 'notperiodic'},
            massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.05},
            sfr_kwargs = 'm0',
            fig_file = 'test_massevol_sfr_m0_t_notperiodic_dutycycle.png'
            )
    '''
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; SFR(M,t); Periodic SF Duty Cycle",  
            sfrevol_prop = {'name': 'squarewave', 'freq_range': [100.0, 1000.0], 'phase_range': [0, 1]},
            massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.01},
            sfr_kwargs = 'm',
            fig_file = 'test_massevol_sfr_m_t_periodic_dutycycle.png'
            )
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; SFR(M,t); Periodic SF Duty Cycle",  
            sfrevol_prop = {'name': 'squarewave', 'freq_range': [100.0, 1000.0], 'phase_range': [0, 1]},
            massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.01},
            sfr_kwargs = 'm',
            fig_file = 'test_massevol_sfr_m_t_periodic_dutycycle.png'
            )
