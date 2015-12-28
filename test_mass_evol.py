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

def test_analytical(): 
    '''
    Quick analystical test for the integration scheme
    '''

    def logsfr(logmass, t0, **kwargs): 
        return np.log10(2.0 * np.sqrt(10.0**logmass))
    
    mass0 = 0.0 
    t_arr = np.arange(0.0, 14., 1.0)
    
    rk4_masses, euler_masses = [0.0], [0.0]
    rk4_sfrs, euler_sfrs = [0.0], [0.0]

    for tt in t_arr[1:]:

        rk4_mass_f, rk4_sfr_f = mass_evol.integrated_rk4(
                logsfr, 
                mass0, 
                np.min(t_arr), 
                tt, 
                f_retain = 1.0, 
                delt = 0.01
                )
        euler_mass_f, euler_sfr_f = mass_evol.integrated_euler(
                logsfr, 
                mass0, 
                np.min(t_arr), 
                tt, 
                f_retain = 1.0, 
                delt = 0.005
                )

        rk4_masses.append(rk4_mass_f)
        euler_masses.append(euler_mass_f)
        
        rk4_sfrs.append(rk4_sfr_f)
        euler_sfrs.append(euler_sfr_f)

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(1, figsize=(10,12))
    sub = fig.add_subplot(111)
    
    analytic = np.log10((t_arr * 10**9.)**2)
    sub.plot(t_arr, (rk4_masses-analytic)/analytic, lw=4, c=pretty_colors[1], label='RK4')
    sub.plot(t_arr, (euler_masses-analytic)/analytic, lw=4, ls='--', c=pretty_colors[3], label='Euler')
    sub.set_xlim([0.0, 14.0])
    sub.set_xlabel('fake t', fontsize=25)
    sub.set_ylabel('(fake integrated M - fake analytic M)/ fake analytic M', fontsize=25)
    #sub.set_ylim([17.0, 20.0])
    sub.legend(loc='lower right')
    plt.show()
    return None

def test_mass_evol(mass0, t0, tf, 
        delmass0 = 0.1, 
        t_initial = 4.0079, 
        sfrevol_param = [1.0, 1.0],
        sfrevol_prop = {'name': 'notperiodic'}, 
        massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.025}, 
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
    
    sfrkwargs = {'z0': z0, 'mass0': mass0, 'deltamass0': delmass0}
    mass_f, sfr_f = mass_evol.integrated(
            massevol_prop['type'], 
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
        massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.05},
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
        
        for tf_i in tfs: 
            sub.vlines(tf_i, -1.0, 2.5, color='k', linewidth=1, linestyle=':')

    m_sub.set_xlim([t0-2.5, tf+0.5])
    m_sub.set_ylim([8.5, 12.5])
    m_sub.set_ylabel(r"$\mathtt{log\;M_*}$", fontsize=20)
    m_sub.set_title(title)

    sub.set_xlim([t0-2.5, tf+0.5])
    sub.set_ylim([-2.0, 2.0])
    sub.set_xlabel(r"$\mathtt{t_{cosmic}}$", fontsize=20)
    sub.set_ylabel(r"$\mathtt{log\;SFR}$", fontsize=20)
    sub.legend(loc='upper left')

    if fig_file is not None: 
        figfile = ''.join([
            'figure/', 
            fig_file
            ])
        fig.savefig(figfile, bbox_inches='tight')
    #plt.show()
    plt.close()

if __name__=="__main__": 
    #test_analytical()
    mass0 = np.arange(7.5, 12.5, 0.5)
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; SFR(M,t); Constant SF Duty Cycle",  
            sfrevol_prop = {'name': 'notperiodic'},
            massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.05},
            fig_file = 'test_euler_massevol_midfreq_notperiodic_dutycycle_twoslopeSFMS.png'
            )
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; SFR(M,t); Periodic SF Duty Cycle",  
            sfrevol_prop = {'name': 'squarewave', 'freq_range': [2.*np.pi, 20.*np.pi], 'phase_range': [0, 1]},
            massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.01},
            fig_file = 'test_euler_massevol_midfreq_periodic_dutycycle_twoslopeSFMS.png'
            )
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; SFR(M,t); NewAmp Periodic SF Duty Cycle",  
            sfrevol_prop = {'name': 'newamp_squarewave', 'freq_range': [1.*np.pi, 10.*np.pi], 'phase_range': [0,1], 'sigma': 0.3},
            massevol_prop = {'name': 'integrated', 'type': 'euler', 'f_retain': 0.6, 't_step': 0.01},
            fig_file = 'test_euler_massevol_midfreq_newamp_periodic_dutycycle_twoslopeSFMS.png'
            )


