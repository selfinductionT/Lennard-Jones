import numpy as np
import matplotlib.pyplot as plt

tics = []
kin_energy = []
pot_energy = []

with open("energy.csv", "r") as f:
    for line in f:
        line = line.split(';')
        tics.append(float(line[0]))
        pot_energy.append(float(line[1]))
        kin_energy.append(float(line[2]))


tics = np.asarray(tics)
kin_energy = np.asarray(kin_energy)
pot_energy = np.asarray(pot_energy)

# figure, axis = plt.subplots(2, 2)

# axis[0, 0].scatter(tics, pot_energy)
# axis[0, 0].set_title("potential energy")

# axis[0, 1].scatter(tics, kin_energy)
# axis[0, 1].set_title("kinetic energy")

# axis[1, 0].scatter(tics, pot_energy + kin_energy)
# axis[1, 0].set_title("mechanical energy")

plt.xlabel(r'Время работы программы, тиков', fontsize=14)
plt.ylabel(r'Энергия, у.е.', fontsize=14)

plt.title(r'График зависимости разных видов энергии от времени',
          fontsize=14)
plt.grid(True)

plt.errorbar(tics, pot_energy + kin_energy, fmt='o',
             color='black', capsize=3, label=r'Полная механическая энергия')

plt.errorbar(tics, kin_energy, fmt='o',
             color='red', capsize=3, label=r'Кинетическая энергия частиц')

plt.errorbar(tics, pot_energy, fmt='o',
             color='green', capsize=3, label=r'Потенциальная энергия частиц')

plt.legend(loc='best', fontsize=12)

plt.show()
