# N_body.py

'''
Class that holds info for N-body dynamics

'''
import numpy as np

G = 4*(np.pi)**2

class N_body:
    '''
    Represents an N-body system
    
    Parameters
    ----------
    masses : 
        type = ndarray
        size = (N, 1) 
        units = Msun
    positions : 
        type = ndarray
        size = (N, d) 
        units = AU
    velocities : 
        type = ndarray
        size = (N, d) 
        units = AU / yr
    epsilon : softening factor
        type = int or float (optional)
    radials : (radii, radial densities) 
        type = ndarray (optional)
        size = (N, 2) 
        units = (AU, Msun/AU)
    radius : cluster radius
        type = int or float (optional)
        units = AU
    '''
    
    def __init__(self, masses, positions, velocities, epsilon=None, radials=None, radius=None):
        
        self.ms = np.asarray(masses, dtype=np.float64)
        self.ps = np.asarray(positions, dtype=np.float64)
        self.vs = np.asarray(velocities, dtype=np.float64)
        self.N, self.dimensions = self.ps.shape
        self.radials = radials
        
        if epsilon is not None:
            self.eps = float(epsilon)
        elif radius is not None:
            self.eps = 0.01 * radius / (self.N ** (1/3))
        else:
            self.eps = 0.0
        
        self.kE, self.pE, self.E, self.v_ratio = self.energy()
        self.accels = self.acceleration()
    
    def acceleration(self):
        """
        Gravity equation: G * m_j (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
        Returns acceleration of each body
        
        Parameters
        ----------
        self : N_body object
            - N
            - masses
            - positions
            - epsilon
        """
        
        N = self.N 
        ms = self.ms.ravel()
        ps = self.ps  

        # make a (N,N,d) matrix to hold all position differences
        diff = ps[None, :, :] - ps[:, None, :] 
        
        # multiply by masses
        m_diff = diff * ms[None, :, None]
        
        # radii squared
        r2 = np.sum(diff**2, axis=-1)
        
        # set diagonal to inf to avoid 0/0 self-force errors
        np.fill_diagonal(r2, np.inf)
        
        # denominator
        r2_e = (r2 + self.eps**2)**(-1.5)

        # accelerations, have to reshape denominator for matrix multiplication
        # sums up all accelerations in each component for each body
        accel = G * np.sum(m_diff * r2_e[:, :, None], axis=1)

        return accel
 
    def energy(self):
        """
        Kinetic energy: 1/2 m * v**2
        Potential energy: - G M m / r
        Total energy: E = KE + PE
        Virial ratio: |2KE + PE| / |PE|
        Returns KE, PE, E, v_ratio for entire system
    
        ****** currently set up with loops O(N_body^2) *********
        
        Parameters
        ----------
        self : N_body object
            - masses
            - positions
            - velocities
        """
        
        ms = self.ms.ravel()
        ps = self.ps
        vs = self.vs
    
        v2 = np.sum(vs**2, axis=1)
        
        kE = np.sum(1/2 * ms * v2)
        
        diff = ps[None, :, :] - ps[:, None, :]   
        r2 = np.sum(diff**2, axis=-1)
        r_e = np.sqrt(r2 + self.eps**2)
        
        np.fill_diagonal(r_e, np.inf)
        
        mm = ms[:, None] * ms[None, :] 
        
        pE = -G * np.sum(mm / r_e) / 2 
        
        v_ratio = abs(2*kE + pE) / abs(pE)
        
        return kE, pE, kE+pE, v_ratio
     
    def update_acceleration(self):
        """
        Updated accelerations 
        
        Parameters
        ----------
        self : N_body object
        """
        
        self.accels = self.acceleration()
    
    def update_energies(self):
        """
        Updated energies
        
        Parameters
        ----------
        self : N_body object
        """
        
        self.kE, self.pE, self.E, self.v_ratio = self.energy()