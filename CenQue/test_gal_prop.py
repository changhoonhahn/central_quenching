'''

Test galactic properties

'''
import numpy as np 
# --- Local --- 
import sfr_evol
import gal_prop as gp 
import util.util as Util
from plotting.plots import PlotSMF
# --- Plotting --- 
import matplotlib.pyplot as plt 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors 

def PlotModelSMF(source='li-march'):
    ''' Plot the different model SMFs for redshift range 0.0 - 1.0
    '''
    zrange = [0.1, 0.3, 0.5, 0.7, 0.9]
    if isinstance(source, str): 
        source_list = [source]
    elif isinstance(source, list): 
        source_list = source
    else: 
        raise ValueError
    lstyles = ['-', '--', '-.', ':']

    smfplot = PlotSMF()
    for i_s, source in enumerate(source_list): 
        for iz, z in enumerate(zrange): 
            lbl_str = None 
            if iz == 0: 
                if i_s ==0: 
                    lbl_str = '$z='+str(round(z,2))+'\;'+source+'$'
                else: 
                    lbl_str = source
            else: 
                if i_s == 0:
                    lbl_str = 'z='+str(round(z,2))

            smfplot.model(z, 
                    source=source, 
                    line_style=lstyles[i_s],
                    line_color=smfplot.pretty_colors[iz], 
                    label=lbl_str)

    smfplot.set_axes()

    fig_file = ''.join([
        'figure/test/', 
        'ModelSMF', 
        '.', '_'.join(source_list), 
        '.png'])
    smfplot.save_fig(fig_file)
    return None

def PlotFq_t(lit='wetzelsmooth'): 
    ''' Plot the quiescent fraction evolution as a function of M* and t
    '''
    prettyplot()
    pretty_colors = prettycolors()
    
    M_arr = np.arange(9.0, 12.5, 0.5)
    t_arr = np.arange(5.0, 14.0, 0.2)

    fig = plt.figure()
    sub = fig.add_subplot(111)
    for Mstar in M_arr:  
        fq_M_t = [] 
        for tt in t_arr: 
            qf = gp.Fq()
            fq_M_t.append(qf.model(np.array([Mstar]), Util.get_zsnap(tt), lit=lit)) 

        sub.plot(t_arr, fq_M_t, label=r"$M_* = "+str(round(Mstar,2))+"$") 
    sub.set_xlim([t_arr.min(), t_arr.max()])
    sub.set_ylim([0., 1.])
    fig_file = ''.join(['figure/test/', 'fq_t.', lit, '.png'])
    fig.savefig(fig_file, bbox_inches='tight') 

    return None

def Plot_dFqdt():  
    ''' Plot the evolution of the derivative of the quiescent fraction as a function of M* and t
    for different parameterized models of f_Q

    d f_Q / d t 

    '''
    prettyplot()
    pretty_colors = prettycolors()
    dt = 0.01
    
    M_arr = np.arange(9.0, 12.5, 0.5)
    t_arr = np.arange(6.0, 13.0, 0.2)

    fig = plt.figure()
    sub = fig.add_subplot(111)

    for iM, Mstar in enumerate(M_arr):  
        qf = gp.Fq()
        lit = 'wetzelsmooth'
        sub.plot(
            t_arr,            
            [gp.dFqdt(Mstar, tt, lit=lit, dt=dt) for tt in t_arr],
            ls='-', lw=3, 
            c=pretty_colors[iM],
            label=r"$\mathtt{M}_* = "+str(round(Mstar,2))+"$"
            ) 

        lit = 'wetzel'
        sub.plot(
            t_arr,            
            [gp.dFqdt(Mstar, tt, lit=lit, dt=dt) for tt in t_arr],
            ls=':', lw=3, 
            c=pretty_colors[iM]
            ) 


        lit = 'cosmosinterp'
        t_blah = np.arange(8.0, 9.6, 0.1)
        sub.plot(
                t_blah, 
                [gp.dFqdt(Mstar, tt, lit=lit, dt=dt) for tt in t_blah],
                ls='--', lw=3, 
                c=pretty_colors[iM]
                ) 

        lit = 'cosmosfit'
        t_blah = np.arange(8.0, 9.6, 0.1)
        sub.plot(
                t_blah, 
                [gp.dFqdt(Mstar, tt, lit=lit, dt=dt) for tt in t_blah],
                ls='-.', lw=3,
                c=pretty_colors[iM]
                ) 

    sub.set_xlim([t_arr.min(), t_arr.max()])
    sub.legend(loc='upper left')
    sub.set_ylim([0., 0.2])
    #plt.show()
    fig_file = ''.join(['figure/test/', 
        'dfqdt.png'])
    fig.savefig(fig_file, bbox_inches='tight') 
    return None

if __name__=='__main__': 
    PlotFq_t(lit='wetzelsmooth')
    PlotFq_t(lit='wetzel')
    Plot_dFqdt()#lit='wetzelsmooth')
    #PlotModelSMF(source=['li-march']) 
    #PlotModelSMF(source=['li-drory-march'])
    #PlotModelSMF(source=['li-drory-march_sameslope'])
