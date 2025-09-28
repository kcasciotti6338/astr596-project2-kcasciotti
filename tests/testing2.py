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

# ------------------ 2 Body Testing ------------------ # 

G = 4*(np.pi)**2
r = 1.0 #AU
Msun = 1.0
Mearth = MEARTH/MSUN
v_circular = np.sqrt(G *Msun / r)
v_x = 0.0
v_y = 2*np.pi
Tt = 100.0 #yr
dt = 0.01
times = np.arange(0, Tt, dt)

position_S_E = np.array([[0.0,0.0], [r,0.0]])
velocity_S_E = np.array([[0.0,0.0], [v_x,v_y]])
mass_S_E = np.array([[Msun], [Mearth]])

S_E = N_body(mass_S_E, position_S_E, velocity_S_E)

results1 = integrate(S_E, Tt, dt, 'euler')

fig1, ax1 = ast.plot_orbit(results1, 'Sun Earth: Test 1', mass_S_E)

# Bad fig 6

fig2, ax2 = ast.plot_over_time(times, results1.KEs, 'Kinetic Energy')
fig2, ax2 = ast.plot_over_time(times, results1.PEs, 'Potenial Energy', ax2)
fig2, ax2 = ast.plot_over_time(times, results1.Es, 'Total Energy', ax2)

# Bad fig 7, 8

# Works for 2 bodies for euler

results2 = integrate(S_E, Tt, dt, 'RK2')

fig1, ax1 = ast.plot_orbit(results2, 'Sun Earth: Test 2', mass_S_E)
fig2, ax2 = ast.plot_over_time(times, results2.KEs, 'Kinetic Energy')
fig2, ax2 = ast.plot_over_time(times, results2.PEs, 'Potenial Energy', ax2)
fig2, ax2 = ast.plot_over_time(times, results2.Es, 'Total Energy', ax2)

# Honestly, looks pretty good. RK4 is very good for 2 bodies over a short timescale

results3 = integrate(S_E, Tt, dt, 'RK4')

fig1, ax1 = ast.plot_orbit(results3, 'Sun Earth: Test 3', mass_S_E)
fig2, ax2 = ast.plot_over_time(times, results3.KEs, 'Kinetic Energy')
fig2, ax2 = ast.plot_over_time(times, results3.PEs, 'Potenial Energy', ax2)
fig2, ax2 = ast.plot_over_time(times, results3.Es, 'Total Energy', ax2)

# Hmmmmmm, RK2 did significantly better than RK4, that doesnt seem right
# Fixed:    accidently did 2*K4 as well, fixed and now looks even better

results4 = integrate(S_E, Tt, dt, 'RK4')

fig1, ax1 = ast.plot_orbit(results4, 'Sun Earth: Test 4', mass_S_E)
fig2, ax2 = ast.plot_over_time(times, results4.KEs, 'Kinetic Energy')
fig2, ax2 = ast.plot_over_time(times, results4.PEs, 'Potenial Energy', ax2)
fig2, ax2 = ast.plot_over_time(times, results4.Es, 'Total Energy', ax2)

plt.show()