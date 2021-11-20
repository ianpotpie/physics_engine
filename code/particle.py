import numpy as np


class Particle:
    def __init__(self, position=None, velocity=None, radius=1, mass=1, charge=0):
        default_dimensions = 2
        if position is None:
            position = np.zeros(default_dimensions)
        if velocity is None:
            velocity = np.zeros(default_dimensions)
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.radius = radius
        self.mass = mass
        self.charge = charge

    def set_object(self, particle_object):
        self.position = np.array(particle_object["position"])
        self.velocity = np.array(particle_object["velocity"])
        self.radius = particle_object["radius"]
        self.mass = particle_object["mass"]
        self.charge = particle_object["charge"]

    def get_object(self):
        particle_object = {
            "position": self.position.tolist(),
            "velocity": self.position.tolist(),
            "radius": self.radius,
            "mass": self.mass,
            "charge": self.charge
        }
        return particle_object

    def get_momentum(self):
        """
        Calculates the momentum of the particle from its velocity and mass
        :return: the momentum vector of the particle
        """
        momentum_vector = self.position * self.velocity
        return momentum_vector

    def get_energy(self):
        """
        Calculates the kinetic energy of the particle based on its velocity and mass
        :return: the KE of the particle
        """
        velocity_squared = np.dot(self.velocity, self.velocity)
        energy = 0.5 * self.mass * velocity_squared
        return energy

    def time_step(self, force, delta_t):
        """
        updates the velocity and position of the particles based on a time step and a force applied to the particle of
        the time step
        :param force: a force applied to the particle over the course of the time step
        :param delta_t: the interval over which the time step occurs
        :return: None
        """
        delta_v = force * delta_t / self.mass
        old_velocity = self.velocity
        new_velocity = self.velocity + delta_v
        self.velocity = new_velocity
        # the position is updated with the average velocity throughout its acceleration
        self.position = self.position + (delta_t * (old_velocity + new_velocity) / 2)
