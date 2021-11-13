from pprint import pprint
import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

from physics_system import PhysicsSystem
from particle import Particle

# ------------------------------------------------------------
# set up initial state


particle_set = {Particle(1, 1, [np.random.normal(), np.random.normal()]) for _ in range(20)}
system = PhysicsSystem(2, particle_set)

dt = 1. / 30  # 30fps

# ------------------------------------------------------------
# set up figure and animation
x_lower = system.boundaries[0][0]
x_upper = system.boundaries[0][1]
length = x_upper - x_lower
y_lower = system.boundaries[1][0]
y_upper = system.boundaries[1][1]
height = y_upper - y_lower

fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(x_lower - (0.1 * length), x_upper + (0.1 * length)),
                     ylim=(y_lower - (0.1 * height), y_upper + (0.1 * height)))

# particles holds the locations of the particles
particles, = ax.plot([], [], 'bo', ms=6)

# rect is the box edge
rect = plt.Rectangle((x_lower, y_lower), length, height, ec='none', lw=2, fc='none')
ax.add_patch(rect)


def init():
    """initialize animation"""
    global system, rect
    particles.set_data([], [])
    rect.set_edgecolor('none')
    return particles, rect


def animate(i):
    """perform animation step"""
    global system, rect, dt, ax, fig
    system.time_step(dt)
    # print("Kinetic Energy: " + str(system.get_kinetic_energy()))
    # print("Gravitational Potential: " + str(system.get_gravitational_potential()))
    print("Total Energy: " + str(system.get_kinetic_energy() + system.get_gravitational_potential()))

    ms = int(fig.dpi * 2 * 0.14 * fig.get_figwidth()
             / np.diff(ax.get_xbound())[0])

    # update pieces of the animation
    rect.set_edgecolor('k')
    x_vals = [p.position[0] for p in particle_set]
    y_vals = [p.position[1] for p in particle_set]
    particles.set_data(x_vals, y_vals)
    particles.set_markersize(ms)
    return particles, rect


ani = animation.FuncAnimation(fig, animate, interval=10, blit=True, init_func=init)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
# ani.save('particle_box.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
