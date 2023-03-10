import numpy as np


# TODO: we dont need velocity as atribute, but we have
class Particle():
    def __init__(self, start_coordinates, coordinates):
        self.radius_vector = coordinates
        self.previous = start_coordinates

    def move(self, F, size, dt):
        r_new = 2*self.radius_vector - self.previous + F*dt**2
        r_new %= size
        self.previous, self.radius_vector = self.radius_vector, r_new

    def mv_and_get_velocity(self, F, size, dt):
        r_new = 2*self.radius_vector - self.previous + F*dt**2
        velocity = (r_new - self.radius_vector)
        r_new %= size
        self.previous, self.radius_vector = self.radius_vector, r_new
        return velocity


class Box():
    def __init__(self, size, particles, dt):
        self.size = size
        self.particles = particles
        self.dt = dt

    def make_box(N, size, velocity, dt):
        particles = []
        alpha = np.ceil(np.cbrt(N))
        x0 = size / alpha  # size of cell
        n = 0

        for i in np.arange(alpha):
            for j in np.arange(alpha):
                for k in np.arange(alpha):
                    r0 = x0 * (np.array([i, j, k]) + 0.5)
                    r1 = r0 + velocity*dt
                    particles.append(Particle(r0, r1))
                    n += 1

                    if n == N:
                        return Box(size, particles, dt)

        return Box(size, particles, dt)

    def random_v(N, size, velocity, dt):
        particles = []
        alpha = np.ceil(np.cbrt(N))
        x0 = size / alpha  # size of cell
        n = 0

        for i in np.arange(alpha):
            for j in np.arange(alpha):
                for k in np.arange(alpha):
                    r0 = x0 * (np.array([i, j, k]) + 0.5)
                    vel = np.random.random(3) * velocity
                    # vel - velocity of new particle
                    r1 = r0 + vel*dt
                    particles.append(Particle(r0, r1))
                    n += 1

                    if n == N:
                        return Box(size, particles, dt)

        return Box(size, particles, dt)

    # TODO: this is not beautifull at all
    def move(self):
        forces = [0. for i in np.arange(len(self.particles))]
        for i in np.arange(len(self.particles)):
            particle = self.particles[i]

            for k in np.arange(i+1, len(self.particles)):
                interacting = self.particles[k]
                delta_r = interacting.radius_vector - particle.radius_vector

                delta_r[delta_r > self.size / 2] -= self.size
                delta_r[delta_r < self.size / 2] += self.size

                F = Box.simple_force(delta_r)
                forces[i] += F
                forces[k] -= F
            particle.move(forces[i], self.size, self.dt)

    def simple_force(delta_r):
        module = np.linalg.norm(delta_r)
        module_F = 24 * (np.power(module, -7) - 2 * np.power(module, -13))
        return (module_F * (delta_r / module))
