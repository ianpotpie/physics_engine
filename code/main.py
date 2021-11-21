import numpy as np
from particle import Particle
from physics_system import PhysicsSystem
from physics_engine import PhysicsEngine
from data_processing import data_to_images

if __name__ == '__main__':
    # # TODO: This creates 2 large particles with initial velocities and run the animation to test collisions
    # particles = {
    #     Particle(mass=50, radius=20, position=[300, 350], velocity=[15, 0]),
    #     Particle(mass=50, radius=20, position=[300, 250],velocity=[-15, 0])
    # }
    # system = PhysicsSystem(dimensions=2, particles=particles)
    # system.boundaries = [[0, 600], [0, 600]]
    # system.boundary_types = [["reflect", "reflect"], ["reflect", "reflect"]]
    # e = PhysicsEngine(system)
    # e.live_animation(framerate=60, updates_per_frame=10, sample_times=None)

    # TODO: generates a random spread of particles with no initial velocities
    n = 10  # the number of particles to generate
    m = 50  # the mass of each particle
    r = 2  # the radius of each particle
    x = 64  # the x position around which the particles are initialized
    y = 64  # the y position around which the particles are initialized
    v = 10  # how spread out the particles are at initialization
    rand = np.random
    particles = {Particle(mass=m, radius=r, position=[rand.normal(x, v), rand.normal(y, v)]) for _ in range(n)}
    system = PhysicsSystem(dimensions=2, particles=particles)
    system.boundaries = [[0, 128], [0, 128]]
    system.boundary_types = [["reflect", "reflect"], ["reflect", "reflect"]]
    e = PhysicsEngine(system)
    # e.live_animation(framerate=60, updates_per_frame=10, sample_times=None)

    sample_times = list(range(10))
    timestamp = e.run_interval(0.1, 10, sample_times)
    data_to_images(str(timestamp))
