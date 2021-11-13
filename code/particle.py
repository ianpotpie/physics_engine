import numpy as np


class Particle:
    def __init__(self, mass, radius, position, charge=0, velocity=None):
        dimensions = len(position)
        if velocity is None:
            velocity = np.zeros(dimensions)
        self.mass = mass
        self.charge = charge
        self.radius = radius
        self.position = np.array(position)
        self.velocity = np.array(velocity)

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
        self.position = self.position + (delta_t * (old_velocity + new_velocity) / 2)
