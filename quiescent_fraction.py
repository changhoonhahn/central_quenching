import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
import sys
import pyfits

class QuiescentFraction: 
    def __init__(self, mass, qf, qf_low, qf_high): 
        self.mass = mass
        self.qf = qf
        self.qf_low = qf_low
        self.qf_high = qf_high
    def interpol_qf(self, mstar): 
        qf_interp = interpolate.interp1d(self.mass,self.qf)
        return qf_interp(mstar) 
    def writetoFile(self, filename):
        f = open(filename, 'w') 
        for i in range(len(self.mass)):
            f.write(str(np.log10(self.mass[i]))+'\t'+str(self.qf[i])+'\t'+str(self.qf_low[i])+'\t'+str(self.qf_high[i])+'\n')
        f.close()
    def bestfit(self, function = 'exp'): 
        if function == 'exp': 
            popt, pcov = curve_fit(expon, np.log10(self.mass), self.qf) 
        elif function == 'linear': 
            popt, pcov = curve_fit(linear, np.log10(self.mass), self.qf)
        elif function == 'gauss': 
            popt, pcov = curve_fit(gauss, np.log10(self.mass), self.qf)
        return popt

def expon(x, sig): 
    return np.exp((x-12.0)/sig) 

def gauss(x, sig):
    return np.exp(-0.5*(x-12.0)**2/sig**2)

def linear(x, x0, a, b): 
    return a+b*(x-x0) 

def upmass(mass): 
    if mass < 9.5: 
        upmass = 9.5
    elif mass < 10.0: 
        upmass = 10.0
    elif mass < 10.5: 
        upmass = 10.5
    elif mass < 11.0: 
        upmass = 11.0
    elif mass < 11.5: 
        upmass = 11.5
    else: 
        print "Mstar is out of range"
    return upmass

def readCosmos(z,type='all'): 
    dir = '/data1/hahn/wetzel_tree/quenched_fractions/'
    # get COSMOS quiescent fraction file name
    # set zlabel
    if z == 0.36: 
        zlabel = 1
    elif z == 0.66:
        zlabel = 2
    elif z == 0.88:
        zlabel = 3
    else: 
        raise NameError
    # set galaxy type
    if (type != 'all' and type != 'cen' and type != 'sat'): 
        print 'not one of the types'
    filename = ''.join(['stats_z', str(zlabel), '.fq_', type])

    cosmos_qf = np.loadtxt(''.join([dir,filename]))
    return cosmos_qf 

if __name__=="__main__": 
    for z in [0.36, 0.66, 0.88]: 
        for type in ['cen']: #,'all', 'sat']: 
            cosmosdata = readCosmos(z, type)
            cosmos = QuiescentFraction(cosmosdata[:,0], cosmosdata[:,1], cosmosdata[:,2], cosmosdata[:,3])
            cosmos.writetoFile(''.join(['/data1/hahn/wetzel_tree/qf_z',str(z),type,'.dat']))
            print 'z ~', z, 'best fit sigma', cosmos.bestfit()
