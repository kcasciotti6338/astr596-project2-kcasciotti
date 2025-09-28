# sampling_methods.py

'''
Random, Kroupa IMF, Plummer Sphere

'''

import numpy as np
from n_body import N_body
from constants import MSUN, MEARTH

G = 4*(np.pi)**2

def Earth_Sun(Tt=10, dt=0.1):
    """
    Returns an N-body object of the Sun-Earth system
        
    Parameters
    ----------
    Tt : total time
        type = int or float
        units = yr
        default : 10 yr
    dt : step size
        type = int or float
        units = AU
        default : 0.1 AU
    """
    r = 1.0 #AU
    Msun = 1.0
    Mearth = MEARTH/MSUN
    v_x = 0.0
    v_y = 2*np.pi

    ps = np.array([[0.0,0.0], [r,0.0]])
    vs = np.array([[0.0,0.0], [v_x,v_y]])
    ms = np.array([[Msun], [Mearth]])

    return N_body(ms, ps, vs)

def Earth_Sun_Jupiter(Tt=10, dt=0.1):
    """
    Returns an N-body object of the Sun-Earth-Jupiter system
        
    Parameters
    ----------
    Tt : total time
        type = int or float
        units = yr
        default : 10 yr
    dt : step size
        type = int or float
        units = AU
        default : 0.1 AU
    """
    Rearth = 1.0 #AU
    Rjupiter = 5.2 #AU
    Msun = 1.0
    Mearth = MEARTH/MSUN
    Mjupiter = 9.55e-4
    v_x_earth = 0.0
    v_y_earth = 2*np.pi
    v_x_jupiter = 0.0
    v_y_jupiter = -2*np.pi*Rjupiter / 12
    
    ps = np.array([[0.0,0.0], [Rearth,0.0], [-Rjupiter,0.0]])
    vs = np.array([[0.0,0.0], [v_x_earth,v_y_earth], [v_x_jupiter,v_y_jupiter]])
    ms = np.array([[Msun], [Mearth], [Mjupiter]])

    return N_body(ms, ps, vs)

def Nbody(N, a=100, m_min=0.08, m_max=150, epsilon=0):
    """
    Returns an N-body system
    
    Generates mass with Kroupa_IMF
    Generates positions, velocities with Plummer Sphere
        
    Parameters
    ----------
    N : number of bodies
    a : half-mass radius
        type = float
        units = AU
        default : 100 AU
    m_min : minimum mass
        type = float
        units = Msun
        default : 0.08 Msun
    m_max : maximum mass
        type = float
        units = Msun
        default : 150 Msun 
    epsilon : softening factor
        type = float
        units = AU
        default : 0
    """
    
    def Kroupa_IMF(N, m_min=0.08, m_max=150):
        """
        Kroupa IMF sampling
        Returns an (N, 1) array of masses in Msun
            
        Parameters
        ----------
        N : number of bodies
        m_min : minimum mass
            type = float
            units = Msun
            default : 0.08 Msun
        m_max : maximum mass
            type = float
            units = Msun
            default : 150 Msun 
        """
        
        masses = np.ones([N,1])
        
        alpha1 = 1.3
        alpha2 = 2.3
        
        m_min_alpha = m_min **(1-alpha1)
        m_b_alpha1 = 0.5 **(1-alpha1)
        m_b_alpha2 = 0.5 **(1-alpha2)
        m_max_alpha = m_max **(1-alpha2)
        
        k1 = 1.0
        k2 = k1 * 0.5**(alpha2 - alpha1)
        
        i1 = k1 * (m_b_alpha1 - m_min_alpha) / (1 - alpha1)
        i2 = k2 * (m_max_alpha - m_b_alpha2) / (1 - alpha2)
        
        p1 = i1 / (i1 + i2)
        
        u = np.random.uniform(0, 1, N)
        w = np.random.uniform(0, 1, N)
        low = u < p1
        
        masses[low, 0] = (m_min_alpha + w[low]*(m_b_alpha1 - m_min_alpha))**(1/(1-alpha1))
        masses[~low, 0] = (m_b_alpha2 + w[~low]*(m_max_alpha - m_b_alpha2))**(1/(1-alpha2))
        
        return masses

    def Plummer_sphere(N, ms, a=100):
        """
        Plummer sphere position sampling
        Returns an (N, 3) array of positions in AU
        Returns an (N, 3) array of velocities in AU / yr
            
        Parameters
        ----------
        N : number of bodies
        ms : (N, 1) array of masses
        a : half-mass radius
            type = float
            units = AU
            default : 100 AU
        """
        
        # ---------------- initialization ---------------- #
        
        M = ms.sum()
        ps = np.zeros((N, 3), dtype=float)
        vs = np.zeros((N, 3), dtype=float)
        
        # ---------------- internal functions ---------------- #

        def enclosed_mass(r):
            """
            Enclosed mass for verrification
            
            Parameters
            ----------
            r : radius
            units = AU
            """
            return M * r**3 / (r**2 + a**2)**(3/2)
        
        def center_mass(ms, ps):
            """
            Returns center of mass as an (1, 3) array
            
            Parameters
            ----------
            ms : (N, 1) array of masses
            ps : (N, 3) array of positions
            """
            m = ms.ravel()
            m_ps = m[:, None] * ps
            ps_center = np.sum(m_ps, axis=0) / ms.sum()
            
            return ps_center
        
        def sphere_to_cartesian(r, theta, phi):
            """
            Translates spherical coordinates to cartesian coordinates
                
            Parameters
            ----------
            r : radius
            theta : 
                units = radians
            phi : 
                units = radians
            """
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            return x, y, z
        
        # ---------------- implementation ---------------- #
        
        u = np.random.uniform(1e-12, 1 - 1e-12, N)
        v = np.random.uniform(0, 1, N)
        w = np.random.uniform(0, 1, N)
        
        r = a * (u**(-2/3) - 1)**(-1/2)  
        
        radial_density = ((3*M / (4*np.pi * a**3)) * (1+ (r/a)**2))**(-5/2)
        
        radials = np.array([r, radial_density])
          
        phi = 2 * np.pi * v
        theta = np.arccos(1 - 2*w)
        ps[:,0], ps[:,1], ps[:,2] = sphere_to_cartesian(r, theta, phi)
        
        ps_center = center_mass(ms, ps)
        
        ps_centered = ps - ps_center
        
        sigma = (G*M / (6*a)) * (1 + (r**2 / a**2))**(-1/2)
        
        vs[:, 0] = np.random.uniform(0.0, np.sqrt(sigma))
        vs[:, 1] = np.random.uniform(0.0, np.sqrt(sigma))
        vs[:, 2] = np.random.uniform(0.0, np.sqrt(sigma))
        
        vs_center = center_mass(ms, vs)
        vs_centered = vs - vs_center
   
        return ps_centered, vs_centered, radials
    
    # ---------------- implementation ---------------- #
    
    ms = Kroupa_IMF(N, m_min, m_max)
    ps, vs, radials = Plummer_sphere(N, ms, a)  
    
    return N_body(ms, ps, vs, epsilon, radials, a)
