#!/usr/local/bin/python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path_txt = os.path.dirname(__file__) + "/../dump_files/"
path_img = os.path.dirname(__file__) + "/../images/"

# plt.figure(figsize=(5,5))
for fname in os.listdir(path_txt):
    if fname.endswith('.txt'):
        with open(path_txt + fname, 'r') as f:
            # size, iteration, grid = [line for line in f]
            grid = [line for line in f][0]
            size = 6
            size = int(size)
            # iteration = int(iteration)
            # print(grid)
            grid = np.array([float(x) for x in grid.split()])

            plt.figure(frameon=False)
            plt.imshow(grid.reshape(size, size))

            fname = fname[:-4]

            print fname

            plt.savefig(path_img + fname, transparent=True, bbox_inches='tight', pad_inches=0)
            plt.close()
