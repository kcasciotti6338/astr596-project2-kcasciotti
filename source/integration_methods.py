# integration_methods.py

'''
All of the inegration methods needed for n-body dynamics

Euler's, RK2, RK4, Leapfrog

'''
import numpy as np
import copy
from n_body import N_body
from dataclasses import dataclass

@dataclass
class Results:
    '''
    Represents an integrated N-body system
    
    Parameters
    ----------
    ts : times
        type = ndarray
        size = (steps, 1) 
        units = yr
    ps : positions
        type = ndarray
        size = (steps, d) 
        units = AU
    vs : velocities
        type = ndarray
        size = (steps, d) 
        units = AU / yr
    kEs : kinetic energies
        type = ndarray
        size = (steps, 1) 
        units = Msun AU^2 yr^-2
    pEs : potential energies
        type = ndarray
        size = (steps, 1) 
        units = Msun AU^2 yr^-2
    Es : totel energies
        type = ndarray
        size = (steps, 1) 
        units = Msun AU^2 yr^-2
    v_ratios : virial ratios
        type = ndarray
        size = (steps, 1) 
        units = None
    E_errors : energy errors
        type = ndarray
        size = (steps, 1) 
        units = None
    '''
    ts: np.ndarray 
    ps: np.ndarray       
    vs: np.ndarray        
    kEs: np.ndarray     
    pEs: np.ndarray      
    Es: np.ndarray  
    v_ratios: np.ndarray  
    E_errors: np.ndarray 

def integrate(NBody, duration=10, dt=0.1, method='euler', snap_linear=None, snap_log=None):
    """
    Integrates an N-body system
    Returns filled Results class
    
    By default saves all steps.
    Snapshots saves only specified number of frames, either linearly or logarithmicaly spaced
        
    Parameters
    ----------
    NBody : N_body object
    duration : overall time
        type = int or float
        units = yr
        default : 10 yr
    dt : step size
        type = float
        units = yr
        default : 0.1 yr  
    method : name of integration function
        type = str
        options: euler, RK2, RK4, leapfrog 
        default : euler 
    snap_linear : number of saved frames (linearly spaced)(optional)
        type = int
        default : None (all frames saved)
    snap_log : number of saved frames (logarithmicaly spaced)(optional)
        type = int
        default : None (all frames saved)
    """
    # ---------------- initialization ---------------- #
    
    N = NBody.N
    d = NBody.dimensions
    n_steps = int(duration / dt)
    
    if snap_linear is not None:
        saves = np.linspace(0, n_steps, snap_linear, dtype=int)
    elif snap_log is not None:
        saves = np.geomspace(0, n_steps, snap_log, dtype=int)
    else:
        saves = np.arange(0, n_steps+1, 1, dtype=int)
    
    ts  = np.arange(len(saves), dtype=float) * dt
    ps = np.zeros((len(saves), N, d), dtype=float)
    vs = np.zeros((len(saves), N, d), dtype=float)
    kEs = np.zeros(len(saves), dtype=float)
    pEs = np.zeros(len(saves), dtype=float)
    Es  = np.zeros(len(saves), dtype=float)
    v_ratios = np.zeros(len(saves), dtype=float)
    E_errors = np.zeros(len(saves), dtype=float)
    
    # ---------------- integration methods ---------------- #
    
    def euler(NBody, dt):
        """
        Euler's method of integration
        
        Parameters
        ----------
        NBody : N_body object
            - positions
            - velocities
            - accelerations
        dt : step size
            type = float
            units = yr    
        """
    
        system_k = copy.copy(NBody)
        system_k.ps = NBody.ps + dt * NBody.vs
        system_k.vs = NBody.vs + dt * NBody.accels
        system_k.update_acceleration()
    
        return system_k
    
    def RK2(NBody, dt):
        """
        RK2 method of integration
        
        Parameters
        ----------
        NBody : N_body object
            - positions
            - velocities
            - accelerations
        dt : step size
            type = float
            units = yr   
        """
    
        system_k1 = copy.copy(NBody)
        system_k1.ps = NBody.ps + dt/2 * NBody.vs
        system_k1.vs = NBody.vs + dt/2 * NBody.accels
        system_k1.update_acceleration()
        
        system_k2 = copy.copy(system_k1)
        system_k2.ps = NBody.ps + dt * system_k1.vs
        system_k2.vs = NBody.vs + dt * system_k1.accels
        system_k2.update_acceleration()
    
        return system_k2
     
    def RK4(NBody, dt):
        """
        RK4 method of integration
        
        Parameters
        ----------
        NBody : N_body object
            - positions
            - velocities
            - accelerations
        dt : step size
            type = float
            units = yr   
        """
    
        system_k1 = copy.copy(NBody)
        
        system_k2 = copy.copy(NBody)
        system_k2.ps = NBody.ps + dt/2 * system_k1.vs
        system_k2.vs = NBody.vs + dt/2 * system_k1.accels
        system_k2.update_acceleration()
        
        system_k3 = copy.copy(NBody)
        system_k3.ps = NBody.ps + dt/2 * system_k2.vs
        system_k3.vs = NBody.vs + dt/2 * system_k2.accels
        system_k3.update_acceleration()
        
        system_k4 = copy.copy(NBody)
        system_k4.ps = NBody.ps + dt * system_k3.vs
        system_k4.vs = NBody.vs + dt * system_k3.accels
        system_k4.update_acceleration()
        
        system_f = copy.copy(NBody)
        system_f.ps = NBody.ps + dt/6 * (system_k1.vs + 2*system_k2.vs + 2*system_k3.vs + system_k4.vs)
        system_f.vs = NBody.vs + dt/6 * (system_k1.accels + 2*system_k2.accels + 2*system_k3.accels + system_k4.accels)
        system_f.update_acceleration()
        
        return system_f
    
    def leapfrog(NBody, dt):
        """
        symplectic integrators : leapfrog
        
        Parameters
        ----------
        NBody : N_body object
            - positions
            - velocities
            - accelerations
        dt : step size
            type = float
            units = yr   
        """
    
        system_1_2 = copy.copy(NBody)
        system_1_2.vs = NBody.vs + dt/2 * NBody.accels
        
        system_2 = copy.copy(NBody)
        system_2.ps = NBody.ps + dt * system_1_2.vs
        system_2.vs = system_1_2.vs
        system_2.update_acceleration()
        
        system_2.vs = system_2.vs + dt/2 * system_2.accels
        system_2.update_acceleration()
        
        return system_2
    
    # ---------------- integration implementation ---------------- #
    
    methods = {
        'euler': euler,
        'RK2': RK2,
        'RK4': RK4,
        'leapfrog': leapfrog 
    }
    
    step = methods.get(method)
    
    system_j = copy.copy(NBody)
    
    ps[0] = system_j.ps
    vs[0] = system_j.vs
    kEs[0], pEs[0], Es[0], v_ratios[0] = system_j.energy()
    k = 0
    
    for j in range(1, n_steps+1):
        system_j = step(system_j, dt)
        
        if j in saves:
            k += 1
            ps[k] = system_j.ps
            vs[k] = system_j.vs
            system_j.update_energies()
            kEs[k], pEs[k], Es[k], v_ratios[k] = system_j.energy()
            E_errors[k] = (Es[k] - Es[0]) / abs(Es[0])

    return Results(ts=saves*dt, ps=ps, vs=vs, kEs=kEs, pEs=pEs, Es=Es, v_ratios=v_ratios, E_errors=E_errors)