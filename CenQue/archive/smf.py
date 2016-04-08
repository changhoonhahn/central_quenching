'''

Stellar Mass Function (SMF) of objects

'''
import numpy as np 

from sham_hack import SMFClass 

from util.cenque_utility import get_z_nsnap

class SMF(object): 
    def __init__(self, **kwargs): 
        ''' Class that describes the SMF of class objects
        '''
        self.kwargs = kwargs
        self.mass = None
        self.phi = None

    def Obj(self, obj, dlogm=None, box=None, h=0.7):
        ''' Calculate the SMF for any Class Object with attributes
        mass 

        Parameters
        ----------
        obj : 
            Any object with the attribute 'mass'
        dlogm : 
            log M bin size (default is 0.1)
        box : 
            Box length (default is 250 Mpc/h)
        '''
        return self._smf(obj.mass, dlogm=dlogm, box=box, h=h)

    def analytic(self, redshift, dlogm='', source='li-drory-march'): 
        '''
        Analytic best-fit SMF from Li-Drory-March (Andrew's Treepm Code)

        Return
        ------
        [masses, phi] : 
            array of masses, array of phi (smf number density) 
        '''
        if not dlogm: 
            dlogm = 0.1
        m_arr = np.arange(6.0, 12.1, dlogm)

        if redshift < 0.1:
            redshift = 0.1

        MF = SMFClass(source=source, redshift=redshift)
        
        mass, phi = [], [] 
        for mi in m_arr: 
            mass.append(mi + 0.5 * dlogm)
            phi.append(MF.numden(mi, mi+dlogm)/dlogm)

        #print 'Analytic ', np.sum(np.array(phi))
        return np.array(mass), np.array(phi)

    def _smf(self, masses, dlogm=None, box=None, h=0.7): 
        '''
        Calculate SMF given masses
        '''
        if not dlogm: 
            dlogm = 0.1 # log M bin size
        if not box: 
            box = 250   # 250 Mpc/h Box length

        m_arr = np.arange(0.0, 12.1, dlogm)

        vol = box ** 3  # box volume
            
        mass, phi = [], [] 
        for mi in m_arr: 
            mass_bin = np.where(
                    (masses > mi) & 
                    (masses <= mi+dlogm)
                    )

            ngal_bin = len(masses[mass_bin])
            
            mass.append(mi + 0.5 * dlogm)
            phi.append(np.float(ngal_bin)/vol/dlogm * h**3) 

        return np.array(mass), np.array(phi)

def test_smf(): 
    '''
    Quick code to test the SMF
    '''
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    # central subhalo catalog 
    #censub = Subhalos(scatter = 0.0)
    #censub.build_catalogs(snapshots=[20])
    censub.read(1)

    smf = SMF()
    mass, phi = smf.centralsubhalos(censub)
    sub.plot(mass, phi, lw=4, c='gray') 
    print np.sum(phi)

    analytic_mass, analytic_phi = smf.analytic(get_z_nsnap(1)) 
    sub.plot(analytic_mass, analytic_phi, lw=4, ls='--', c='k') 
    print np.sum(analytic_phi)

    sub.set_yscale('log')
    sub.set_ylim([10**-8, 10**-1])
    sub.set_xlim([7.0, 12.0])

    plt.show()

def test_smf_evol(type='analytic'): 
    '''
    Quick code to test the evolution of SMF
    '''
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    for isnap in range(1, 21): 
    
        # central subhalo catalog 
        #censub = Subhalos(scatter = 0.0)
        ##censub.build_catalogs(snapshots=[20])
        #censub.read(isnap)

        smf = SMF()
        #mass, phi = smf.centralsubhalos(censub)
        #sub.plot(mass, phi, lw=4, c='gray') 
        #print np.sum(phi)
        
        if type == 'analytic': 
            analytic_mass, analytic_phi = smf.analytic(get_z_nsnap(isnap)) 
            sub.plot(analytic_mass, analytic_phi, lw=4, ls='--', c='k') 
            print np.sum(analytic_phi)

    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([8.75, 12.0])

    plt.show()

if __name__ == '__main__': 
    test_smf_evol()
