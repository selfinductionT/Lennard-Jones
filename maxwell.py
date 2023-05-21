import numpy as np
import matplotlib.pyplot as plt

v_x = []

with open("maxwell.csv", "r") as f:
    for line in f:
        v_x.append(float(line))

v_x = np.sort(np.asarray(v_x))
N = len(v_x)

T = np.sum(v_x**2)
plt.hist(v_x, bins=N, density=True)

print(T/N)
plt.show()
