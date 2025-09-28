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
Tt = 10.0 #yr
dt = 0.01
times = np.arange(0, Tt, dt)

position_S_E = np.array([[0.0,0.0], [r,0.0]])
velocity_S_E = np.array([[0.0,0.0], [v_x,v_y]])
mass_S_E = np.array([[Msun], [Mearth]])

S_E = N_body(mass_S_E, position_S_E, velocity_S_E)
#energy_S_E = S_E.energy()

#  Fixed:   divide by actual radius instead of position vector
#           PE and E need to be passed as floats, not 1D arrays
#           works with N^2 loops, should vectorize -> later problem

accel_S_E = S_E.acceleration()

# Fixed:    check for dimentions instead of hard coding 2 or 3
#           reset acceleration to 0 for each body in the loop

S_E.update()

fig1, ax1 = ast.plot_snapshot(S_E, 'SUN_EARTH test1')

S_E_t = integrate(S_E, Tt, dt, 'euler')

# Fixed:    in the integrate function, the system_t array was initialized 
#               to only hold floats, now it holds N-body objects

fig2, ax2 = ast.plot_orbit(S_E_t, 'SUN_EARTH')

# Fixed:    had trouble acceccing and all the positions from array of N-body objects
#               decided to stick with loops since plotting doesn't need to be optimized

# Fixed:    nothing is moving, X and Y are switched
#           Bad fig 1, 2
#           realized i was using G in cgs

# Fixed:    realized i had to use deepcopy instead of just renaming the same thing

# Fixed:    something is definitly wrong with my acceleration function
#           spent a very long time ~7 hrs on vectorizing this
#           works now though, yay

KE, PE, E = S_E.energy()



fig3, ax3 = ast.plot_over_time(times, KE, 'KE')

plt.show()

# Fixed:    had a lot of trouble accessing variables after integration 
#           updated integration to be a dataclass
#           3D arrays are easier than hiding info in the N_body class within an array
# Switching to another test file because i'm making too many stuctural changes