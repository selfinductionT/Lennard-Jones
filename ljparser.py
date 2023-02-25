import numpy as np
import ljmodel as model


def parse(file_name):
    # this function parses input file

    N = 0  # number of particles
    size = np.array([0, 0, 0])
    velocity_random = True  # equals True if velocity will be random
    velocity = np.array([0., 0., 0.])

    with open(file_name, "r") as file:
        for line in file:
            if len(line.rstrip()) == 0 or line[0] == '#':
                continue

            line = line.split()
            k_word = line[0]
            if k_word == 'N':
                N = int(line[2])
            if k_word == 'dt':
                dt = float(line[2])
            if k_word == 'size':
                size = float(line[2])

            if k_word == 'velocity':
                if line[2] == 'random':
                    velocity = float(line[3])
                else:
                    velocity_random = False
                    velocity = np.array([
                        float(line[2]),
                        float(line[3]),
                        float(line[4])
                    ])

    if velocity_random:
        return model.Box.random_v(N, size, velocity, dt)
    else:
        return model.Box.make_box(N, size, velocity, dt)
