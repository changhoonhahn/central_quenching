'''

Test integrated mass evolution 

'''
import numpy as np
from scipy import interpolate

from smf import SMF
from util.cenque_utility import get_z_nsnap

from defutility.plotting import prettyplot
from defutility.plotting import prettycolors 

def analytic_smf_evol(source='li-drory-march'): 
    '''
    Evolution of the analytic SMF evolution  
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    for i_z, z in enumerate(np.arange(0.0, 1.5, 0.25)):
        smf = SMF()
        analytic_mass, analytic_phi = smf.analytic(z, source=source) 
        sub.plot(analytic_mass, analytic_phi, lw=4, ls='--', c=pretty_colors[i_z], 
                label=r"$ z = "+str(round(z,2))+"$") 

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([8.0, 12.0])
    sub.legend(loc='upper right')
    fig.savefig(
            ''.join([
                'figure/'
                'analytic_smf_', source, '.png'
                ]), 
            bbox_inches='tight'
            )
    plt.close()

if __name__=="__main__": 
    analytic_smf_evol(source='fontana')
    analytic_smf_evol(source='li-march')
