import numpy as np
import matplotlib.pyplot as plt

v_x = []

with open("maxwell.csv", "r") as f:
    for line in f:
        v_x.append(float(line))

v_x = np.sort(np.asarray(v_x))
N = len(v_x)

F, bins = np.histogram(v_x, bins=N, density=True)
print(len(bins))
v_x = np.asarray([(bins[i] + bins[i+1])/2 for i in range(N)])
F = F[:N//2]
v_x = v_x[:N//2]

plt.scatter(v_x**2, np.log(F))
plt.show()
