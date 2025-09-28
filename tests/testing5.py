# testing.py

'''
Testing for Project 2

2 body, 3 body, N-body (A,B,C)

'''

import numpy as np
from N_body import N_body
from integration_methods import integrate
import plotting as ast
from matplotlib import pyplot as plt
from sampling_methods import sample_Nbody

# ------------------ N-Body Testing: B,C ------------------ # 

N = 500
Tt = 2000.0
dt = 100
a = 1000

#masses = Kroupa_IMF(N)
#positions = np.random.uniform(0, 100, [N,3])
#velocities = np.zeros([N,3])

system = sample_Nbody(N, a, m_min=0.08, m_max=150)

resultsLF = integrate(system, Tt, dt, 'leapfrog')
fig3, ax3 = ast.plot_orbit(resultsLF, 'N-body: Test C - leapfrog', system.ms)
fig4, ax4 = ast.plot_over_time(resultsLF.ts, resultsLF.kEs, 'Kinetic Energy')
fig4, ax4 = ast.plot_over_time(resultsLF.ts, resultsLF.pEs, 'Potenial Energy', ax4)
fig4, ax4 = ast.plot_over_time(resultsLF.ts, resultsLF.Es, 'Total Energy', ax4)
fig4, ax4 = ast.plot_over_time(resultsLF.ts, resultsLF.v_ratios, 'Virial Ratio')

plt.show()

# Bad fig 15 - 21 for B
# Bad fig 22 - 24 for C

#something is wrong, my radius is getting huge
#forgot to include softening - fixed

