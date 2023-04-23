import numpy as np


class Particle():
    def __init__(self, start_coordinates, coordinates):
        self.radius_vector = coordinates
        self.previous = start_coordinates

    def mv_and_get_velocity(self, F, size, dt):
        delta = self.radius_vector - self.previous
        delta[delta > 0.75*size] -= size
        delta[delta < -0.75*size] += size

        r_new = self.radius_vector + delta + F*dt**2
        velocity = (r_new - self.radius_vector)/dt
        r_new %= size

        self.previous, self.radius_vector = self.radius_vector, r_new
        return velocity


class Box():
    def __init__(self, N, size, particles, dt):
        self.size = size
        self.particles = particles
        self.dt = dt
        self.N = N

    def random_v(N, size, velocity, dt):
        particles = []
        alpha = np.ceil(np.cbrt(N))
        x0 = size / alpha  # size of cell
        velocities = np.random.rand(N, 3) * velocity
        velocities -= velocities[::-1, :]
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

    # TODO: this is not beautifull at all
    def move(self):
        forces = [0.] * self.N
        pot_energy = 0.
        kin_energy = 0.

        for i in np.arange(self.N):
            particle = self.particles[i]

            for k in np.arange(i+1, self.N):
                interacting = self.particles[k]
                delta_r = interacting.radius_vector - particle.radius_vector

                delta_r[delta_r > self.size / 2] -= self.size
                delta_r[delta_r < -(self.size / 2)] += self.size

                F = Box.simple_force(delta_r)
                forces[i] += F
                forces[k] -= F

                pot_energy += Box.potential_energy(delta_r)

            velocity = particle.mv_and_get_velocity(
                forces[i], self.size, self.dt)
            kin_energy += np.sum(np.square(velocity)) / 2

        return (pot_energy, kin_energy)

    def potential_energy(delta_r):
        module = np.linalg.norm(delta_r)
        return (4 * (np.power(module, -12) - np.power(module, -6)))

    def simple_force(delta_r):
        module = np.linalg.norm(delta_r)
        module_F = 24 * (np.power(module, -7) - 2 * np.power(module, -13))
        return (module_F * (delta_r / module))
