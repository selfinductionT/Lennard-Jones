import ljparser as parser
import ljmodel as model
import numpy as np
import matplotlib.pyplot as plt

kin, a = parser.parse("input.txt")
pot = 0.

for i in np.arange(a.N):
    particle = a.particles[i]

    for k in np.arange(i+1, a.N):
        interacting = a.particles[k]
        delta_r = interacting.previous - particle.previous

        delta_r[delta_r > a.size / 2] -= a.size
        delta_r[delta_r < a.size / 2] += a.size

        pot += model.Box.potential_energy(delta_r)


kin_energy = np.zeros(11)
pot_energy = np.zeros(11)
sum_energy = np.zeros(11)

kin_energy[0] = kin
pot_energy[0] = pot
sum_energy[0] = kin + pot

print(pot, kin, kin+pot)

tics = np.linspace(1, 11, num=11)

for i in np.arange(10):
    pot_energy[i+1], kin_energy[i+1], sum_energy[i+1] = a.move(True, True)
    print(pot_energy[i+1], kin_energy[i+1], sum_energy[i+1])


figure, axis = plt.subplots(2, 2)

axis[0, 0].plot(tics, pot_energy)
axis[0, 0].set_title("potential energy")

axis[0, 1].plot(tics, kin_energy)
axis[0, 1].set_title("kinetic energy")

axis[1, 0].plot(tics, pot_energy)
axis[1, 0].set_title("mechanical energy")

plt.show()
