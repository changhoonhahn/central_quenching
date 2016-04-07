import numpy as np
import sys
import matplotlib.pyplot as plt
class deltamass: 
    def __init__(self, nsnap): 
        dir = '/data1/hahn/wetzel_tree/'
        file_prefix = 'delta_mass_evolved_sham_snapshot'
        file_ext    = '.dat'
        filename = ''.join([dir, file_prefix, str(nsnap), file_ext])
        self.name = filename
    
        output = np.loadtxt(filename)
        self.dmass = output[:,0]
        self.shammass = output[:,1]
        self.sfrmass = output[:,2]

if __name__=="__main__": 
    figdir = '/home/users/hahn/research/figures/tinker/'
    for i in range(1,13): 
        dmass_nsnap = deltamass(i)
        print 'mean', np.mean(dmass_nsnap.dmass), 'median', np.median(dmass_nsnap.dmass), \
                'min', np.min(dmass_nsnap.dmass), 'max', np.max(dmass_nsnap.dmass)
        print i 
        fig1 = plt.figure(1, figsize=(7, 7)) 
        ax1 = fig1.add_subplot(111)
        ax1.scatter(dmass_nsnap.shammass, dmass_nsnap.sfrmass, color='black')
        ax1.plot([8.5+0.1*np.float(j) for j in range(100)], [8.5+0.1*np.float(j) for j in range(100)], color='red')
        ax1.set_xlim([8.5, 12.0])
        ax1.set_ylim([8.5, 12.0])
        ax1.set_xlabel('SHAM Mass')
        ax1.set_ylabel('SFR Intregal Mass')
        ax1.text(9.0, 11.5, ''.join(['Snapshot ', str(i)]))
        ax1.grid(True)
        figname = ''.join([figdir, 'fig_delta_mass_sham_sfrintegral_snapshot', str(i), '.png'])
        fig1.savefig(figname)
        fig1.clear()
