import time
import numpy as np

start = time.time()

a = [0 for i in np.arange(10000)]
print(time.time() - start)

start = time.time()

a = [0] * 10000
print(time.time() - start)
