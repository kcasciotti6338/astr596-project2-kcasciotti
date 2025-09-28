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

# ------------------ 3 Body Testing ------------------ # 

G = 4*(np.pi)**2
Rearth = 1.0 #AU
Rjupiter = 5.2 #AU
Msun = 1.0
Mearth = MEARTH/MSUN
Mjupiter = 9.55e-4
v_x_earth = 0.0
v_y_earth = 2*np.pi
v_x_jupiter = 0.0
v_y_jupiter = -2*np.pi*Rjupiter / 12
Tt = 100.0 #yr
dt = 0.01
times = np.arange(0, Tt, dt)

position_S_E_J = np.array([[0.0,0.0], [Rearth,0.0], [-Rjupiter,0.0]])
velocity_S_E_J = np.array([[0.0,0.0], [v_x_earth,v_y_earth], [v_x_jupiter,v_y_jupiter]])
mass_S_E_J = np.array([[Msun], [Mearth], [Mjupiter]])

S_E_J = N_body(mass_S_E_J, position_S_E_J, velocity_S_E_J)

results1 = integrate(S_E_J, Tt, dt, 'euler')
fig1, ax1 = ast.plot_orbit(results1, 'Sun Earth Jupiter: Test 1', mass_S_E_J)
fig2, ax2 = ast.plot_over_time(times, results1.kEs, 'Kinetic Energy')
fig2, ax2 = ast.plot_over_time(times, results1.pEs, 'Potenial Energy', ax2)
fig2, ax2 = ast.plot_over_time(times, results1.Es, 'Total Energy', ax2)
fig3, ax3 = ast.plot_over_time(times, results1.v_ratios, 'Virial Ratio')

results2 = integrate(S_E_J, Tt, dt, 'RK2')
fig4, ax4 = ast.plot_orbit(results2, 'Sun Earth Jupiter: Test 2', mass_S_E_J)
fig5, ax5 = ast.plot_over_time(times, results2.kEs, 'Kinetic Energy')
fig5, ax5 = ast.plot_over_time(times, results2.pEs, 'Potenial Energy', ax5)
fig5, ax5 = ast.plot_over_time(times, results2.Es, 'Total Energy', ax5)
fig6, ax6 = ast.plot_over_time(times, results2.v_ratios, 'Virial Ratio')

results3 = integrate(S_E_J, Tt, dt, 'RK4')
fig7, ax7 = ast.plot_orbit(results3, 'Sun Earth Jupiter: Test 3', mass_S_E_J)
fig8, ax8 = ast.plot_over_time(times, results3.kEs, 'Kinetic Energy')
fig8, ax8 = ast.plot_over_time(times, results3.pEs, 'Potenial Energy', ax8)
fig8, ax8 = ast.plot_over_time(times, results3.Es, 'Total Energy', ax8)
fig9, ax9 = ast.plot_over_time(times, results3.v_ratios, 'Virial Ratio')

results4 = integrate(S_E_J, Tt, dt, 'leapfrog')
fig10, ax10 = ast.plot_orbit(results4, 'Sun Earth Jupiter: Test 4', mass_S_E_J)
fig11, ax11 = ast.plot_over_time(times, results4.kEs, 'Kinetic Energy')
fig11, ax11 = ast.plot_over_time(times, results4.pEs, 'Potenial Energy', ax11)
fig11, ax11 = ast.plot_over_time(times, results4.Es, 'Total Energy', ax11)
fig12, ax12 = ast.plot_over_time(times, results4.v_ratios, 'Virial Ratio')

plt.show()

# looks pretty snazzy for 3 body