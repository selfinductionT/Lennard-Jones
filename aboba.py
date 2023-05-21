import ljparser as parser
import ljmodel as model
import numpy as np
import time

steps = 10000  # total number of iterations

kin, a = parser.parse("input.txt")
pot = 0.
displacement = 0.


start = time.time()

# get total potential energy after first step
for i in np.arange(a.N):
    particle = a.particles[i]

    displacement += np.sum(particle.displacement**2)

    for k in np.arange(i+1, a.N):
        interacting = a.particles[k]
        delta_r = interacting.previous - particle.previous

        delta_r[delta_r > a.size / 2] -= a.size
        delta_r[delta_r < -(a.size / 2)] += a.size

        pot += model.Box.potential_energy(delta_r)


kin_energy = np.zeros(steps)
pot_energy = np.zeros(steps)
average_displacement = np.zeros(steps)

kin_energy[0] = kin
pot_energy[0] = pot
average_displacement[0] = displacement/a.N

tics = np.linspace(1, steps, num=steps)

for i in np.arange(1, steps):
    average_displacement[i], pot_energy[i], kin_energy[i] = a.move()
    if (i % 1000) == 0:
        print(i)


# distribution function
velocities_x = a.get_velocities_x()

end = time.time()
print(end - start)


# average_displacement = np.asarray(
#     [np.sum(average_displacement[:i+1], axis=0) for i in np.arange(steps)]
# )

with open("energy.csv", "w") as f:
    for (n, pot, kin) in zip(tics, pot_energy, kin_energy):
        f.write(str(n) + ';' + str(pot) + ';' + str(kin) + '\n')

with open("maxwell.csv", "w") as f:
    for velocity in velocities_x:
        f.write(str(velocity)+'\n')

with open("average_d.csv", "w") as f:
    for d in average_displacement:
        # d = np.linalg.norm(d)
        f.write(str(d) + '\n')
