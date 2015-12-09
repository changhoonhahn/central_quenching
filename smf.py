'''

Stellar Mass Function (SMF) of objects

'''
import numpy as np 

from treepm.sham import SMFClass 

class SMF(object): 

    def __init__(self, **kwargs): 
        '''
        Class that describes the SMF of class objects
        '''
        self.kwargs = kwargs
        self.mass = None
        self.phi = None

    def cenque(self, cq_obj, dlogm='', box=None):
        '''
        Calculate the SMF for given CenQue Class Object

        Parameters
        ----------
        cq_obj : 
            CenQue object
        dlogm : 
            log M bin size (default is 0.1)
        box : 
            Box length (default is 250 Mpc/h)
        '''
        if not dlogm: 
            dlogm = 0.1 # log M bin size
        if not box: 
            box = 250   # 250 Mpc/h Box length

        m_arr = np.arange(8.0, 12.1, dlogm)

        vol = box ** 3
        
        mass, phi = [], [] 
        for m_i in m_arr: 
            mass_bin = np.where(
                    (cq_obj.mass > mi) & 
                    (cq_obj.mass <= mi+dlogm)
                    )

            ngal_bin = len(cq_obj.mass[mass_bin])
            
            mass.append(mi + 0.5 * dlogm)
            phi.append(np.float(ngal_bin)/vol) 

        return np.array(mass), np.array(phi)/dlogm

    def analytic(self, redshift, dlogm=''): 
        '''
        Analytic best-fit SMF from Li-Drory-March (Andrew's Treepm Code)

        Return
        ------
        [masses, phi] : 
            array of masses, array of phi (smf number density) 
        '''
        if not dlogm: 
            dlogm = 0.1
        m_arr = np.arange(8.0, 12.1, dlogm)

        MF = SMFClass(source='li-drory-march', redshift=redshift)
        
        mass, phi = [], [] 
        for mi in m_arr: 
            phi.append(MF.numden(mi, mi+dlogm))
            mass.append(mi + 0.5 * dlogm)

        return np.array(mass), np.array(phi)/dlogm

def test_smf(): 
    '''
    Quick code to test the SMF
    '''
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    smf = SMF()
    mass, phi = smf.analytic(0.1) 

    sub.plot(mass, phi) 
    sub.set_yscale('log')
    sub.set_ylim([10**-5, 10**-1])
    sub.set_xlim([8.75, 12.0])

    plt.show()

if __name__ == '__main__': 
    test_smf()

