import numpy as np
from particle import Particle

from physics_system import PhysicsSystem
from canvas_renderer import Renderer

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    particle_set = {Particle(5000, 10, [np.random.normal(400, 50), np.random.normal(300, 50)]) for _ in range(20)}
    system = PhysicsSystem(2, particle_set)
    r = Renderer(system)
    r.live_animation(framerate=60, updates_per_frame=10, sample_times=[1, 3, 5, 7, 9, 11])
