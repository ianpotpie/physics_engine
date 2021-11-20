import numpy as np
from particle import Particle
from physics_system import PhysicsSystem
from physics_engine import PhysicsEngine
from data_processing import json_to_images

if __name__ == '__main__':
    # generates a random set of particles
    n = 20  # the number of particles to generate
    m = 5000  # the mass of each particle
    r = 10  # the redius of each particle
    x = 400  # the x position around which the particles are initialized
    y = 300  # the y position around which the particles are initialized
    v = 50  # how spread out the particles are at initialization

    # EXAMPLE RUN:
    # initializes a set of particles with the same mass spread out over a normal distribution
    rand = np.random
    particles = {Particle(mass=m, radius=r, position=[rand.normal(x, v), rand.normal(y, v)]) for _ in range(n)}
    system = PhysicsSystem(dimensions=2, particles=particles)
    e = PhysicsEngine(system)
    sample_times = list(range(10))
    timestamp = e.run_interval(0.1, 10, sample_times)
    json_to_images(str(timestamp))

    # r.live_animation(framerate=60, updates_per_frame=10, sample_times=[1, 3, 5, 7, 9, 11])
