from pprint import pprint

import numpy as np
from itertools import combinations
from particle import Particle


class PhysicsSystem:
    def __init__(self, dimensions, particles=None):
        if particles is None:
            particles = set()
        self.dimensions = dimensions  # this will usually just be 2 or 3 dimensions
        self.global_time = 0  # this keeps track of the total amount of time that has passed in the system
        self.particles = particles  # the set of particles in the system
        self.net_force = np.zeros(dimensions)  # a force applied to all particles (may update dynamically)
        self.drag_coefficient = 0  # force felt by particles against their direction of motion ()
        self.gravitational_constant = 1000  # Newton's gravitational constant (or otherwise)
        self.coulomb_constant = 1  # Coulomb's charge attraction constant (or otherwise)
        self.boundaries = [[-20, 20] for _ in range(dimensions)]  # the locations of all boundaries of the system
        self.boundary_types = [["none", "none"] for _ in
                               range(dimensions)]  # this can take values of 'reflect', 'none', 'transport'
        self.collisions = "elastic"

    # Note: the energy of the system will only remain constant if the drag_coefficient is set to zero, the net force is
    # zero, and the boundaries are all set to 'none', otherwise energy is energy will be gained or lost

    def get_momentum(self):
        """
        Calculates the net momentum of the system (this should remain the same in an inertial frame of reference)
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
        gravitational_potential = 0
        for particle1, particle2 in combinations(self.particles, 2):
            G = self.gravitational_constant
            m1 = particle1.mass
            m2 = particle2.mass
            p1 = particle1.position
            p2 = particle2.position
            r = np.linalg.norm(p2 - p1)
            gravitational_potential -= (G * m1 * m2) / r
        return gravitational_potential

    def get_electric_potential(self):
        """
        Calculates the net electrical potential of a system by summing the gravitational potential between all pairs
        of particles
        :return: the net electrical potential of the system
        """
        electrical_potential = 0
        for particle1, particle2 in combinations(self.particles, 2):
            k = self.coulomb_constant
            q1 = particle1.charge
            q2 = particle2.charge
            p1 = particle1.position
            p2 = particle2.position
            r = np.linalg.norm(p2 - p1)
            electrical_potential += (k * q1 * q2) / r
        return electrical_potential

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

        # iterates through pairs of particles and applies inter-particle forces
        for particle1, particle2 in combinations(self.particles, 2):
            # gets the relevant properties of the particle pair
            p1 = particle1.position
            p2 = particle2.position
            p_rel = p2 - p1
            distance = np.linalg.norm(p_rel)

            # this adds the force from gravity
            G = self.gravitational_constant
            m1 = particle1.mass
            m2 = particle2.mass
            F_gravity = (p_rel * G * m1 * m2) / ((distance ** 3) + 0.01)
            particle_forces[particle1] = particle_forces[particle1] + F_gravity
            particle_forces[particle2] = particle_forces[particle2] - F_gravity

            # this adds the force from electrical charge
            k = self.coulomb_constant
            q1 = particle1.charge
            q2 = particle2.charge
            F_charge = (p_rel * k * q1 * q2) / ((distance ** 3) + 0.01)
            particle_forces[particle1] = particle_forces[particle1] + F_charge
            particle_forces[particle2] = particle_forces[particle2] - F_charge

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
                lower_bound_type = self.boundary_types[d][0]
                upper_bound = self.boundaries[d][1]
                upper_bound_type = self.boundary_types[d][0]
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

        # determines the collision behavior
        # it detects when two particles are overlapping - if it finds two which are overlapping then it determines
        # the exact time of the collision by "rewinding" the trajectories of the particle with their current positions
        # and velocities, then it calculates and adjusts the positions and velocities of the particles according to the
        # collision characteristics
        # Assumptions for approximation: the velocity at the time of overlap is the same as the velocity at the time of
        # collision
        if self.collisions == "elastic":
            for particle1, particle2 in combinations(self.particles, 2):
                p1 = particle1.position
                p2 = particle2.position
                r1 = particle1.radius
                r2 = particle2.radius
                p_rel = p2 - p1
                distance = np.linalg.norm(p_rel)

                if distance <= r1 + r2:
                    m1 = particle1.mass
                    m2 = particle2.mass
                    v1 = particle1.velocity
                    v2 = particle2.velocity
                    v_rel = v2 - v1

                    # TODO: fix this block of code
                    # # when two particles are overlapping, their collision actually occurred at some time in the past
                    # # these lines of code determine the time period when that occurred.
                    # # given the velocity and position of the particles, this solves the equation for t:
                    # # ||(v2 - v2)t + (p2 - p1)|| = ||v_rel * t + p_rel|| = (r1 + r2)
                    # a = np.dot(v_rel, v_rel)
                    # b = 2 * ((v_rel[0] * p_rel[0]) + (v_rel[1] * p_rel[1]))
                    # c = r1 + r1 - np.dot(p_rel, p_rel)
                    # # we use the minimum of the quadratic root since we are looking for the time at which the particles
                    # # began overlapping rather than the time at which they will stop overlapping
                    # t_impact = - (b + np.sqrt((b ** 2) - (4 * a * c))) / (2 * a)

                    # the code above is not really working, so we will approximate t_impact to be half the time step
                    t_impact = - delta_t / 2

                    # this gets the particle positions at the time of impact
                    p1_impact = (t_impact * v1) + p1
                    p2_impact = (t_impact * v2) + p2
                    impact_direction = (p2_impact - p1_impact) / np.linalg.norm(p2_impact - p1_impact)

                    # here we project the velocities along the direction of the impact before finding the resulting
                    v1_initial = np.dot(impact_direction, v1) * impact_direction
                    v2_initial = np.dot(impact_direction, v2) * impact_direction

                    v1_final = (v1_initial * (m1 - m2)/(m1 + m2)) + (v2_initial * 2 * m2/(m1 + m2))
                    v2_final = (v1_initial * 2 * m1/(m1 + m2)) + (v2_initial * (m2 - m1)/(m1 + m2))

                    # adjusts the positions and velocities of the particle to account for the impact
                    particle1.velocity = v1 + (v1_final - v1_initial)
                    particle2.velocity = v2 + (v2_final - v2_initial)
                    particle1.position = p1_impact - (t_impact * particle1.velocity)
                    particle2.position = p2_impact - (t_impact * particle2.velocity)






