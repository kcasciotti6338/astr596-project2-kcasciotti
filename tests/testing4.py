# testing.py

'''
Testing for Project 2

2 body, 3 body, N-body (A,B,C)

'''

import numpy as np
from constants import MSUN, MEARTH
from N_body import N_body
from integration_methods import integrate
import plotting as ast
from matplotlib import pyplot as plt
from sampling_methods import k

# ------------------ N-Body Testing: A ------------------ # 

N = 10
Tt = 10.0
dt = 0.1

masses = np.ones([N,1])
positions = np.random.uniform(0, 100, [N,2])
velocities = np.zeros([N,2])
times = np.arange(0, Tt, dt)

system = N_body(masses, positions, velocities)

resultsRK4 = integrate(system, Tt, dt, 'RK4')
fig1, ax1 = ast.plot_orbit(resultsRK4, 'N-body: Test A - RK4', masses)
fig2, ax2 = ast.plot_over_time(times, resultsRK4.kEs, 'Kinetic Energy')
fig2, ax2 = ast.plot_over_time(times, resultsRK4.pEs, 'Potenial Energy', ax2)
fig2, ax2 = ast.plot_over_time(times, resultsRK4.Es, 'Total Energy', ax2)
fig2, ax2 = ast.plot_over_time(times, resultsRK4.v_ratios, 'Virial Ratio')

resultsLF = integrate(system, Tt, dt, 'leapfrog')
fig3, ax3 = ast.plot_orbit(resultsLF, 'N-body: Test A - leapfrog', masses)
fig4, ax4 = ast.plot_over_time(times, resultsLF.kEs, 'Kinetic Energy')
fig4, ax4 = ast.plot_over_time(times, resultsLF.pEs, 'Potenial Energy', ax4)
fig4, ax4 = ast.plot_over_time(times, resultsLF.Es, 'Total Energy', ax4)
fig4, ax4 = ast.plot_over_time(times, resultsLF.v_ratios, 'Virial Ratio')

plt.show()

# Bad figs 11 - 14

