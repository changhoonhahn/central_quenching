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

def test_mass_evol(mass0, t0, tf, sfrevol_param={'name': 'notperiodic'},
        massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.025}): 

    # spline between z and t_cosmic
    z, t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
    z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 
    z0 = z_of_t(t0)
    zf = z_of_t(tf)

    logsfr_mstar_z, sig_logsfr_mstar_z = get_param_sfr_mstar_z()

    def logsfr_notquenched(logmass, t_input, **kwargs): 
        # this is a simplification. logsfr should depend on the new mass at each time step 
        # but that complicates the SF duty cycle 
        zi = kwargs['z0']
        mass_i = kwargs['mass0']
        #avglogsfr = logsfr_mstar_z(mass_i, zi)
        #print 'SFMS(', mass0[0], ',', str(zi), ') ', avglogsfr[0]
        avglogsfr = logsfr_mstar_z(logmass, z0)

        logsfr_sfms = sfr_evol.logsfr_sfms_evol(zi, z_of_t(t_input))

        logsfr_sfduty = sfr_evol.logsfr_sfduty_fluct(
                t_ancestor, 
                t_input, 
                delta_sfr=ancestor.delta_sfr[sf_ancestors[notquenched]],
                sfrevol_param=sfrevol_param_nq, 
                **sfrevol_prop
                )
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        ####properly implement SF DUTY CYCLE
        logsfr_f = avglogsfr + logsfr_sfms# + logsfr_sfduty
        #print np.max(logsfr_f - logsfr_mstar_z(logmass, z_of_t(t_input))), np.min(logsfr_f - logsfr_mstar_z(logmass, z_of_t(t_input)))

        return logsfr_f
    
    sfrkwargs = {'z0': z0, 'mass0': mass0}
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

def plot_integrated_mass_evol(mass0, t0, tf, massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.025}, title=None, fig_file=None):
    
    prettyplot()
    pretty_colors = prettycolors()

    if title is None: 
        raise ValueError 
    delt = (tf-t0)/10.0
    tfs = np.arange(t0+delt, tf+delt, delt)

    z, t = np.loadtxt('snapshot_table.dat', unpack=True, usecols=[2, 3]) 
    z_of_t = interpolate.interp1d(list(reversed(t)), list(reversed(z)), kind='cubic') 

    logsfr_mstar_z, sig_logsfr_mstar_z = get_param_sfr_mstar_z()

    mass_i_1 = mass0
    t_i_1 = t0

    masses, sfrs, sfmssfrs = [mass0], [logsfr_mstar_z(mass0, z_of_t(t0))], [logsfr_mstar_z(mass0, z_of_t(t0))] 
    for tf_i in tfs: 
        mass_i, sfr_i = test_mass_evol(mass_i_1, t_i_1, tf_i, massevol_prop=massevol_prop)
        #mass_i, sfr_i = test_mass_evol(mass0, t0, tf_i, massevol_prop=massevol_prop)

        masses.append(mass_i)
        sfrs.append(sfr_i)
        sfmssfrs.append(logsfr_mstar_z(mass0, z_of_t(tf_i)))
        
        t_i_1 = tf_i
        mass_i_1 = mass_i

    m_fig = plt.figure(1, figsize=(10,12))
    m_sub = m_fig.add_subplot(211)

    tfs = np.hstack([t0, tfs])
    for i in xrange(len(masses[0])):
        m_sub.plot(tfs, np.array(masses)[:,i], c=pretty_colors[i], lw=3, label=str(mass0[i]))
    
    m_sub.set_xlim([t0-2.5, tf+0.5])
    m_sub.set_ylim([8.5, 12.5])
    m_sub.set_ylabel(r"$\mathtt{log\;M_*}$", fontsize=20)
    m_sub.set_title(title)

    sub = m_fig.add_subplot(212)
    for i in xrange(len(masses[0])):
        sub.plot(tfs, np.array(sfrs)[:,i], c=pretty_colors[i], lw=3, label=str(mass0[i]))
        sub.plot(tfs, np.array(sfmssfrs)[:,i], c='k', lw=2, ls='--')
    
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
        m_fig.savefig(figfile, bbox_inches='tight')

    plt.show()

if __name__=="__main__": 
    mass0 = np.arange(7.5, 12.5, 0.5)
    plot_integrated_mass_evol(mass0, 4.0, 13.2, 
            title="Integrated Mass Evolution; No SF Duty Cycle",  
            massevol_prop = {'name': 'integrated', 'f_retain': 0.6, 't_step': 0.05}
            )
