from pprint import pprint

import numpy as np

from particle import Particle


class PhysicsSystem:
    def __init__(self, dimensions, particles=None):
        if particles is None:
            particles = set()
        self.dimensions = dimensions  # this will usually just be 2 or 3 dimensions
        self.global_time = 0  # this keeps track of the total amount of time that has passed in the system
        self.particles = particles  # the set of particles in the system
        self.net_force = np.zeros(dimensions)  # a force applied to all particles (may update dynamically)
        self.drag_coefficient = 1000  # force felt by particles against their direction of motion ()
        self.gravitational_constant = 10  # Newton's gravitational constant (or otherwise)
        self.coulomb_constant = 1  # Coulomb's charge attraction constant (or otherwise)
        self.boundaries = [[-20, 20] for _ in range(dimensions)]  # the locations of all boundaries of the system
        self.boundary_type = [["none", "none"] for _ in
                              range(dimensions)]  # this can take values of 'reflect', 'none', 'transport'

    # Note: the energy of the system will only remain constant if the drag_coefficient is set to zero, the net force is
    # zero, and the boundaries are all set to 'none', otherwise energy is energy will be gained or lost

    def get_momentum(self):
        """
        Calculates the net momentum of the system (this should remain the same in an intertial frame of reference)
        :return: the net momentum vector of the system
        """
        net_momentum = np.zeros(self.dimensions)
        for particle in self.particles:
            net_momentum += particle.get_momentum()
        return net_momentum

    def get_kinetic_energy(self):
        """
        Calculates the net kinetic energy of a system by summing the KE of all particles
        :return: the net kinetic energy of the system
        """
        net_energy = 0
        for particle in self.particles:
            net_energy += particle.get_energy()
        return net_energy

    def get_gravitational_potential(self):
        """
        Calculates the net gravitational potential of a system by summing the gravitational potential between all pairs
        of particles
        :return: the net gravitational potential of the system
        """
        checked = set()
        GPE = 0
        for particle1 in self.particles:
            checked.add(particle1)
            for particle2 in self.particles:
                if particle2 not in checked:
                    G = self.gravitational_constant
                    m1 = particle1.mass
                    m2 = particle2.mass
                    p1 = particle1.position
                    p2 = particle2.position
                    r = np.linalg.norm(p2 - p1)
                    GPE -= (G * m1 * m2) / r
        return GPE

    def get_electric_potential(self):
        """
        Calculates the net electrical potential of a system by summing the gravitational potential between all pairs
        of particles
        :return: the net electrical potential of the system
        """
        checked = set()
        EPE = 0
        for particle1 in self.particles:
            checked.add(particle1)
            for particle2 in self.particles:
                if particle2 not in checked:
                    k = self.coulomb_constant
                    q1 = particle1.charge
                    q2 = particle2.charge
                    p1 = particle1.position
                    p2 = particle2.position
                    r = np.linalg.norm(p2 - p1)
                    EPE += (k * q1 * q2) / r
        return EPE

    def time_step(self, delta_t):
        """
        Updates all the forces, velocities, and positions of the particles in the system based on their current positions
        and velocities over the course of a time interval "delta_t"
        :param delta_t: the time interval over which the system makes its update
        :return: None
        """
        self.global_time += delta_t
        particle_forces = {}

        # adds net force of the system to all particles
        for particle in self.particles:
            particle_forces[particle] = self.net_force

        # iterates through pairs of particles and applies intra-particle forces
        for particle1 in self.particles:
            for particle2 in self.particles:
                if particle2 is not particle1:
                    # gets the
                    p1 = particle1.position
                    p2 = particle2.position
                    r = p2 - p1
                    r_norm = np.linalg.norm(r)

                    # this adds the force from gravity
                    G = self.gravitational_constant
                    m1 = particle1.mass
                    m2 = particle2.mass
                    F_gravity = (r / r_norm) * (G * m1 * m2) / ((r_norm * r_norm) + 0.01)
                    particle_forces[particle1] = particle_forces[particle1] + F_gravity

                    # this adds the force from electrical charge
                    k = self.coulomb_constant
                    q1 = particle1.charge
                    q2 = particle2.charge
                    F_charge = (r / r_norm) * (k * q1 * q2) / ((r_norm * r_norm) + 0.01)
                    particle_forces[particle1] = particle_forces[particle1] + F_charge

        # adds the force due to the viscosity of the environment
        for particle in self.particles:
            velocity = particle.velocity
            F_drag = - self.drag_coefficient * velocity
            particle_forces[particle] = particle_forces[particle] + F_drag

        # imparts the force onto the particles
        for particle in self.particles:
            force = particle_forces[particle]
            particle.time_step(force, delta_t)

        # this determines the behavior at the boundaries
        for particle in self.particles:
            for d in range(self.dimensions):
                position = particle.position[d]
                lower_bound = self.boundaries[d][0]
                lower_bound_type = self.boundary_type[d][0]
                upper_bound = self.boundaries[d][1]
                upper_bound_type = self.boundary_type[d][0]
                if position < lower_bound:
                    if lower_bound_type == "reflect":
                        particle.position[d] = lower_bound + (lower_bound - position)
                        particle.velocity[d] = -1 * particle.velocity[d]
                    elif lower_bound_type == "transport":
                        particle.position[d] = upper_bound - (lower_bound - position)
                elif position > upper_bound:
                    if upper_bound_type == "reflect":
                        particle.position[d] = upper_bound - (position - upper_bound)
                        particle.velocity[d] = -1 * particle.velocity[d]
                    elif upper_bound_type == "transport":
                        particle.position[d] = lower_bound + (position - upper_bound)
