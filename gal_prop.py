'''

Investigate galaxy properties of CenQue project

Author(s): ChangHoon Hahn 


'''

import numpy as np
import matplotlib.pyplot as plt
import mpfit
import scipy.stats as scistat

#----- Local -----
from utility.plotting import prettyplot
from utility.plotting import prettycolors 
import cenque_utility as util 
import cenque as cq 
import cenque_groupcat as cq_group
import sf_mainseq as sfms 

# Stellar Masses ----

def mass_kcorr_isedfit(Mrcut=18): 
    ''' Plot the comparison of mass between k-correct mass and 
    iSEDfit mass 
    
    Parameters
    ----------
    Mrcut : Absolute magnitude cuts to identify SDSS group catalog sample 

    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
        
    centrals = cq_group.central_catalog_match2isedfit(Mrcut=Mrcut, clobber=True) 
    
    fig = plt.figure(1, figsize=(8,8))
    sub = fig.add_subplot(111)

    #sub.scatter(centrals.mass, centrals.kcorrect_mass, s=3, c=pretty_colors[2]) 
    sub.hist2d(centrals.mass, centrals.kcorrect_mass, bins=[30, 30])   # plot (mass, SFR) 

    sub.scatter(np.arange(9.0, 12.0, 0.01), np.arange(9.0, 12.0, 0.01), s=10, c=pretty_colors[4]) 
    sub.set_xlabel('iSEDfit log(Mass)', fontsize=18) 
    sub.set_ylabel('k correct log(Mass)', fontsize=18) 
    sub.set_xlim([9.0, 12.0]) 
    sub.set_ylim([9.0, 12.0]) 
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'mass_kcorr_isedfit_comparison_mrcut', str(Mrcut), '.png']) 
        
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear()

def sfr_mpajhu_isedfit(Mrcut=18): 
    ''' Plot the comparison of SFR between MPA-JHU and 
    iSEDfit SFR 
    
    Parameters
    ----------
    Mrcut : Absolute magnitude cuts to identify SDSS group catalog sample 

    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 
        
    centrals = cq_group.central_catalog_match2isedfit(Mrcut=Mrcut, clobber=True) 
    
    fig = plt.figure(1, figsize=(8,8))
    sub = fig.add_subplot(111)

    #sub.scatter(centrals.sfr, centrals.mpajhu_sfr, s=3, c=pretty_colors[2]) 
    sub.hist2d(centrals.sfr, centrals.mpajhu_sfr, bins=[50, 50])   # plot (mass, SFR) 

    sub.scatter(np.arange(-5.0, 5.0, 0.01), np.arange(-5.0, 5.0, 0.01), s=6, c=pretty_colors[4]) 
    sub.set_xlabel('iSEDfit log(SFR)', fontsize=18) 
    sub.set_ylabel('MPA JHU log(SFR)', fontsize=18) 
    sub.set_xlim([-5.0, 5.0]) 
    sub.set_ylim([-5.0, 5.0]) 
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'sfr_mpajhu_isedfit_comparison_mrcut', str(Mrcut), '.png']) 
        
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear()

def sfr_mass(Mrcut=18): 
    ''' Plot the SFR as a function of mass 
    
    Parameters
    ----------
    Mrcut : Absolute magnitude cuts to identify SDSS group catalog sample 

    '''
    prettyplot()                        #make things pretty 
    pretty_colors = prettycolors() 

    
    fig = plt.figure(1, figsize=(8,8))
    sub = fig.add_subplot(111)
    
    centrals = cq_group.central_catalog(Mrcut=Mrcut, clobber=True) 
    zmid = np.median(centrals.z)
    sub.hist2d(centrals.mass, centrals.sfr, bins=[30, 1000])   # plot (mass, SFR) 
    #sub.scatter(centrals.mass, centrals.sfr, s=4, c=pretty_colors[2])   # plot (mass, SFR) 

    #centrals = cq_group.central_catalog_match2isedfit(Mrcut=Mrcut, clobber=True) 
    #sub.scatter(centrals.mass, centrals.sfr, s=4, c=pretty_colors[3])   # plot (mass, SFR) 

    # plot Salim et al. classification 
    #masses = np.arange(9.0, 12.0, 0.05) 
    #sfrs = (lambda m: -0.49 + 0.65*(m-10.0))(masses)
    #sub.scatter(masses, sfrs, s=6, c=pretty_colors[4]) 

    ##sub.scatter(np.arange(-5.0, 5.0, 0.01), np.arange(-5.0, 5.0, 0.01), s=6, c=pretty_colors[4]) 
    sub.set_ylabel('MPA JHU log(SFR)', fontsize=18) 
    sub.set_xlabel('k correct log(Mass)', fontsize=18) 
    sub.set_xlim([9.0, 12.0]) 
    sub.set_ylim([-5.0, 5.0]) 
    
    fig_file = ''.join(['/home/users/hahn/research/figures/tinker/', 
        'sfr_mass_mrcut', str(Mrcut), '.png']) 
        
    fig.savefig(fig_file, bbox_inches='tight')
    fig.clear()

if __name__=="__main__": 
    #sfr_mass(Mrcut=18)
    #sfr_mass(Mrcut=19)
    #sfr_mass(Mrcut=20)
    #mass_kcorr_isedfit(Mrcut=18)
    #mass_kcorr_isedfit(Mrcut=19)
    #mass_kcorr_isedfit(Mrcut=20)
    sfr_mpajhu_isedfit(Mrcut=18)
    sfr_mpajhu_isedfit(Mrcut=19)
    sfr_mpajhu_isedfit(Mrcut=20)
