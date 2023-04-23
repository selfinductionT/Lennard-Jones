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
        delta_r[delta_r < -(a.size / 2)] += a.size

        pot += model.Box.potential_energy(delta_r)


kin_energy = np.zeros(40)
pot_energy = np.zeros(40)

kin_energy[0] = kin
pot_energy[0] = pot


tics = np.linspace(1, 40, num=40)

for i in np.arange(1, 40):
    pot_energy[i], kin_energy[i] = a.move()


print(tics)

# figure, axis = plt.subplots(2, 2)

# axis[0, 0].scatter(tics, pot_energy)
# axis[0, 0].set_title("potential energy")

# axis[0, 1].scatter(tics, kin_energy)
# axis[0, 1].set_title("kinetic energy")

# axis[1, 0].scatter(tics, pot_energy + kin_energy)
# axis[1, 0].set_title("mechanical energy")

plt.xlabel(r'Время работы программы, тиков', fontsize=14)
plt.ylabel(r'Энергия, у.е.', fontsize=14)

plt.title(r'График зависимости разных видов энергии от времени', fontsize=14)
plt.grid(True)

plt.errorbar(tics, pot_energy + kin_energy, fmt='o',
             color='black', capsize=3, label=r'Полная механическая энергия')

plt.errorbar(tics, kin_energy, fmt='o',
             color='red', capsize=3, label=r'Кинетическая энергия частиц')

plt.errorbar(tics, pot_energy, fmt='o',
             color='green', capsize=3, label=r'Потенциальная энергия частиц')

plt.legend(loc='best', fontsize=12)

plt.show()
