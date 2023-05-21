import numpy as np


class Particle():
    # __slots__ = ("radius_vector", "previous")

    def __init__(self, start_coordinates, coordinates):
        self.radius_vector = coordinates
        self.previous = start_coordinates
        self.displacement = coordinates - start_coordinates

    def mv_and_get_velocity(self, F, size, dt):
        """Move one particle and return its velocity"""
        delta = self.radius_vector - self.previous
        delta[delta > 0.75*size] -= size
        delta[delta < -0.75*size] += size

        displacement = delta + F*dt**2
        self.displacement += displacement
        r_new = self.radius_vector + displacement

        velocity = displacement/dt
        r_new %= size

        self.previous, self.radius_vector = self.radius_vector, r_new
        return velocity


class Box():
    """Box where particles are located"""
    # __slots__ = ("size", "particles", "dt", "N")

    def __init__(self, N, size, particles, dt):
        self.size = size
        self.particles = particles
        self.dt = dt
        self.N = N

    def random_v(N, size, velocity, dt):
        """generate Box with N particles, which velocities are
        less then velocity * np.sqrt(3) """

        particles = []
        alpha = np.ceil(np.cbrt(N))
        x0 = size / alpha  # size of cell
        velocities = np.random.rand(N, 3)
        velocities -= velocities[::-1, :]
        velocities *= velocity
        kin_energy = np.sum(velocities**2/2)
        n = 0

        for i in np.arange(alpha):
            for j in np.arange(alpha):
                for k in np.arange(alpha):
                    r0 = x0 * (np.array([i, j, k]) + 0.5)
                    vel = velocities[n]
                    r = r0 + vel*dt
                    particles.append(Particle(r0, r))
                    n += 1
                    if n == N:
                        return (kin_energy, Box(N, size, particles, dt))

        return (kin_energy, Box(N, size, particles, dt))

    def move(self):
        forces = np.zeros((self.N, 3))
        pot_energy = 0.
        kin_energy = 0.
        average_displacement = 0.

        for i in np.arange(self.N):
            particle = self.particles[i]

            for k in np.arange(i+1, self.N):
                interacting = self.particles[k]
                delta_r = interacting.radius_vector - particle.radius_vector

                delta_r[delta_r > self.size / 2] -= self.size
                delta_r[delta_r < -(self.size / 2)] += self.size

                # F = Box.simple_force(delta_r)
                pot, F = Box.potential_and_force(delta_r)
                forces[i] += F
                forces[k] -= F

                # pot_energy += Box.potential_energy(delta_r)
                pot_energy += pot

            velocity = particle.mv_and_get_velocity(
                forces[i], self.size, self.dt)
            average_displacement += np.sum(particle.displacement**2)

            kin_energy += (np.linalg.norm(velocity)**2) / 2

        return (average_displacement/self.N, pot_energy, kin_energy)

    def potential_and_force(delta_r):
        """returns touple: (potential energy, force)"""
        module = np.linalg.norm(delta_r)
        attraction = -np.power(module, -6)
        repulsion = np.power(attraction, 2)

        potential = repulsion + attraction
        force = potential + repulsion

        potential *= 4
        force *= (-24) * np.power(module, -2) * delta_r

        return potential, force

    def potential_energy(delta_r):
        """Lennard-Jones potential energy of two particles"""
        module = np.linalg.norm(delta_r)
        return 4*(np.power(module, -12) - np.power(module, -6))

    def simple_force(delta_r):
        """Lennard-Jones interraction force of two particles"""
        module = np.linalg.norm(delta_r)
        module_F = 24 * (np.power(module, -8) - 2 * np.power(module, -14))
        return (module_F * delta_r)

    def get_velocities_x(self):
        "returns np array with x-components of particles's velocities"
        delta_x = np.asarray(
            [(p.radius_vector[0] - p.previous[0]) for p in self.particles]
        )

        delta_x[delta_x > 0.75 * self.size] -= self.size
        delta_x[delta_x < -0.75 * self.size] += self.size

        return (delta_x / self.dt)
