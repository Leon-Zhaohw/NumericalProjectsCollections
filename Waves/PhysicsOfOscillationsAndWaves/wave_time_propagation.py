'''
This program is a numerical computation of the
time propagation of a wave in python.
As of now, this program is not working as it should

The scritp was in 2017 translated by Sebastian G. Winther-Larsen from a matlab script originally written by Arnt Inge Vistnes
'''

from matplotlib import pyplot as plt
import numpy as np
from matplotlib import animation

# Parameters
delta_x = 0.1
delta_t = 0.1
sigma = 2.0
v = 0.5
x = np.arange(-20, 20 + delta_x, delta_x)
N = len(x)

# Generating positions at t=0, Gaussian form
u = np.exp(-(x/(2*sigma)) * (x/(2*sigma)))

# Parameters and time derivative at t=0
factor = (delta_t*v / delta_x)**2
print(factor)
dudt = (v / (2 * sigma * sigma)) * x * u

# Initial values
u_jminus1 = u - delta_t * dudt
u_j = u
u_jplus1 = np.zeros(N)

# Set up figure
fig = plt.figure()
ax = plt.axes(xlim=(0, N+1), ylim=(-1.2, 1.2))
line, = ax.plot(u_j, linewidth = 2)

# Initialisation function
def init():

    # Line
    line.set_ydata([])
    return line,

def animate(i):

    global u_j, u_jminus1, u_jplus1

    u_jplus1[2::] = (2 * (1 - factor)) * u_j[1:-1] - \
        u_jminus1[1:-1] + factor * (u_j[2::] + u_j[0:-2])

    # Line
    temp = np.copy(u_j)
    line.set_ydata(temp)

    # Updating
    u_jminus1 = u_j
    u_j = u_jplus1

    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=1000,\
                                interval=10, blit=True)

plt.show()