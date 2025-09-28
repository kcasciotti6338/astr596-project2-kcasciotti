# analysis.py

'''
Analysis for Project 2

'''


import numpy as np
from n_body import N_body
from integration_methods import integrate
import plotting as ast
from matplotlib import pyplot as plt
import sampling_methods as sample
import timeit

def body_2(Tt=10, dt=0.01, saves=None, gif=False):
    '''
    Runs, plots, and saves 2 body system with optional animations
    
    Parameters
    ----------
    Tt : total time (yrs)
    dt : step size (yrs)
    saves : number of saved frames
    gif : make True to save an animation 
    '''
    system = sample.Earth_Sun(Tt, dt)

    src = 'Sun, Earth'
    methods = ['RK2', 'RK4', 'leapfrog']

    results_E = integrate(system, Tt, dt, 'euler', snap_linear=saves)
    fig_E, ax_E = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    fig_E, ax_E = ast.plot_energy_diagnosis(results_E, method='euler', methods=['euler'], src=src, ax=ax_E)
    fig_E, _ = ast.plot_orbit(results_E, masses=system.ms, labl='euler', titl=f'Sun-Earth Orbit over {Tt} years (dt={dt})', ax=ax_E[0, 0])
    if gif: ast.make_orbit_gif(results_E, f'Sun_earth_euler_{dt}dt.gif', titl=f'Sun-Earth Orbit over {Tt} years: euler, {dt}dt')
    plt.savefig(f'Sun_earth_euler_{dt}dt.png')

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))

    for i, method in enumerate(methods):
        results = integrate(system, Tt, dt, method, snap_linear=saves)
        fig, ax = ast.plot_energy_diagnosis(results, method=method, methods=methods, src=src, ax=ax, scale='symlog', add_legends=(i == len(methods)-1))
        fig, _ = ast.plot_orbit(results, masses=system.ms, labl=method, titl=f'Sun-Earth Orbit over {Tt} years (dt={dt})', ax=ax[0, 0])
        if gif: ast.make_orbit_gif(results, f'Sun_earth_{method}_{dt}dt.gif', titl=f'Sun-Earth Orbit over {Tt} years: {method}, {dt}dt')
    
    plt.savefig(f'Sun_earth_{dt}dt.png')

def body_3(Tt=10, dt=0.01, saves=None, gif=False):
    '''
    Runs, plots, and saves 3 body system with optional animations
    
    Parameters
    ----------
    Tt : total time (yrs)
    dt : step size (yrs)
    saves : number of saved frames
    gif : make True to save an animation 
    '''

    system = sample.Earth_Sun_Jupiter(Tt, dt)

    src = 'Sun, Earth, Jupiter'
    methods = ['RK4', 'leapfrog']
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))

    for i, method in enumerate(methods):
        results = integrate(system, Tt, dt, method, snap_linear=saves)
        fig, ax = ast.plot_energy_diagnosis(results, method=method, methods=methods, src=src, ax=ax, scale='symlog', add_legends=(i == len(methods)-1))
        fig, _ = ast.plot_orbit(results, masses=system.ms, labl=method, titl=f'Sun-Earth-Jupiter Orbit over {Tt} years (dt={dt})', ax=ax[0, 0])
        if gif: ast.make_orbit_gif(results, f'Sun_earth_jupiter_{method}_{dt}dt.gif', titl=f'Sun-Earth-Jupiter Orbit: {method}, {dt}dt')
    
    plt.savefig(f'Sun_earth_jupiter_{dt}dt.png')

def body_NA(n, Tt=10, dt=0.01, saves=None):
    '''
    Runs simple N-body with random sampled masses and positions
    100 AU, no softening
    
    Parameters
    ----------
    n : number of bodies
    Tt : total time (yrs)
    dt : step size (yrs)
    saves : number of saved frames
    '''
    
    ms = np.ones([n, 1])
    ps = np.random.uniform(0, 100, [n, 3])
    vs = np.zeros([n, 3])

    system = N_body(ms, ps, vs)

    src = f'{n}_Body_System'
    methods = ['RK4', 'leapfrog']
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    ss = ax[0, 0].get_subplotspec()
    fig.delaxes(ax[0, 0])
    ax[0, 0] = fig.add_subplot(ss, projection='3d')

    for i, method in enumerate(methods):
        results = integrate(system, Tt, dt, method, snap_linear=saves)
        fig, ax = ast.plot_energy_diagnosis(results, method=method, methods=methods, src=src, ax=ax, scale='symlog', add_legends=(i == len(methods)-1))
        fig, _ = ast.plot_snapshot(system, 100, titl=f'Simple {src} at t=0', ax=ax[0, 0])
        
    plt.savefig(f'plots/'+f'Simple_{src}_{Tt}yrs.png')

def body_N(n, Tt, dt, a, saves, eps=None):
    '''
    Runs, creates, and saves N-body analysis and plots
    
    Parameters
    ----------
    n : number of bodies
    Tt : total time (yrs)
    dt : step size (yrs)
    a : half-mass radius (AU)
    saves : number of saved frames
    eps : softening factor (optional) 
    '''
    
    src = f'{n}_Body_System'

    system = sample.Nbody(n, a, epsilon=eps)

    fig_1, ax_1 = ast.histogram(system.ms, xlabl=r'log Kroupa IMF (log $M_\odot$)', titl=f'Histogram: {src}', kroupa=True) 
    fig_1.savefig(f'../figures/Good_Plots/Histograms/'+f'{src}_Histogram.png')
    
    fig_2, ax_2 = ast.plot_density_profile(system.radials, a, titl=f'Radial Density Profile: {src}', labl='Radial Density') 
    fig_2.savefig(f'../figures/Good_Plots/Density_Profiles/'+f'{src}_Dennsity_Profile.png')
    
    fig_3, ax_3 = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    ss = ax_3[0, 0].get_subplotspec()
    fig_3.delaxes(ax_3[0, 0])
    ax_3[0, 0] = fig_3.add_subplot(ss, projection='3d')

    results = integrate(system, Tt, dt, 'leapfrog', snap_linear=saves)
    fig_3, ax_3 = ast.plot_energy_diagnosis(results, method='leapfrog', methods=['leapfrog'], src=src, ax=ax_3, scale='symlog', add_legends=True)
    fig_3, _ = ast.plot_snapshot(system, a, titl=f'{src} at t=0', ax=ax_3[0, 0])
    fig_3.savefig(f'../figures/Good_Plots/Energy_Diagnostics/'+f'{src}_Energy_Diagnostics.png')    
        
    results_5 = integrate(system, Tt, dt, 'leapfrog', snap_linear=5)
    fig_5, ax_5 = ast.many_snapshots(results_5, a, system.ms)
    fig_5.savefig(f'../figures/Good_Plots/Evolution/'+f'{src}_Evolution.png')
      
def run_nbody(n, Tt, dt, a, saves, eps=None):
    '''
    Runs very simple n-body for time analysis
    
    Parameters
    ----------
    n : number of bodies
    Tt : total time (yrs)
    dt : step size (yrs)
    a : half-mass radius (AU)
    saves : number of saved frames
    eps : softening factor (optional) 
    '''
    
    system = sample.Nbody(n, a, epsilon=eps)
    integrate(system, Tt, dt, 'leapfrog', snap_linear=saves)
    
def showcase(n, Tt, dt, a, saves, eps=None):
    '''
    Runs, creates, and saves showcase pngs
    
    Parameters
    ----------
    n : number of bodies
    Tt : total time (yrs)
    dt : step size (yrs)
    a : half-mass radius (AU)
    saves : number of saved frames
    eps : softening factor (optional)   
    '''
    
    src = f'{n}-Body System'

    system = sample.Nbody(n, a, epsilon=eps)

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    ss = ax[0, 0].get_subplotspec()
    fig.delaxes(ax[0, 0])
    ax[0, 0] = fig.add_subplot(ss, projection='3d')

    results = integrate(system, Tt, dt, 'leapfrog', snap_linear=saves)
    fig, ax = ast.plot_energy_diagnosis(results, method='leapfrog', methods=['leapfrog'], src=src, ax=ax)
    fig, _ = ast.plot_snapshot(system, a, titl=f'{src}', ax=ax[0, 0])
    fig.savefig(f'{n}_Body_Showcase.png')
    fig_2, ax_2 = ast.many_snapshots(results, a, system.ms)
    fig_2.savefig(f'{n}_Body_Showcase_Evolution.png')
    
    
# ------------------ 2 Body ------------------ # 

dts = [0.005]

for i, dt in enumerate(dts):
    body_2(10, dt, 200, gif=True)

#f.body_2(1, 0.01, 200, gif=True)

print('--- Completed 2 Body ---')

# ------------------ 3 Body ------------------ # 

dts = [0.02, 0.01, 0.005]

for i, dt in enumerate(dts):
    body_3(24, dt, 200, gif=True)

#f.body_3(12, 0.01, 200, gif=True)

print('--- Completed 3 Body ---')

# ------------------ N Body : A ------------------ # 

Tts = [100, 50]

for i, Tt in enumerate(Tts):
    f.body_NA(100, Tt, 0.001, 200)
    
print('--- Completed Simple N Body ---')

# ------------------ N Body : B ------------------ # 


Ns = [500, 1000, 2000]

for i, n in enumerate(Ns):
    body_N(n, 500, 0.01, 1000, 100)
    
print('--- Completed Large N Body ---')
  

# ------------------ N Body : C ------------------ # 

Ns = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
times = np.zeros(len(Ns))
    
for i, n in enumerate(Ns):
    t_temp = timeit.timeit(lambda: run_nbody(n, 2, 0.01, 1000, 5, 5), number=1)
    times[i] = t_temp
    print(f'--- Completed {n}-Body Time Test ---')
            
fig, ax = ast.plot_preformance(Ns, times)
ax.set_title(f'{Ns[0]} - {Ns[-1]} Body Time Analysis')

fig.savefig(f'plots/'+f'{Ns[0]}_{Ns[-1]}_Body_Time_Analysis.png')

showcase(5_000, 200, 0.1, 1000, 10)

