import numpy as np
import matplotlib.pyplot as plt

d = []
with open("average_d.csv", "r") as f:
    for dr in f:
        dr = float(dr)
        d.append(dr)

d = np.asarray(d)
steps = len(d)
steps = np.linspace(1, steps, num=steps)
steps *= 0.001
plt.scatter(steps, d)
a = np.polyfit(steps, d, deg=1, full=True)

plt.plot(steps, np.polyval(a[0], steps))

print(a)
plt.show()
