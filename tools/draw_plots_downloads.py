#!/usr/local/bin/python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path_txt = os.path.dirname(__file__) + "/../downloads/"
path_img = os.path.dirname(__file__) + "/../images/"

# plt.figure(figsize=(5,5))
for fname in os.listdir(path_txt):
    with open(path_txt + fname, 'r') as f:
        size, iteration, grid = [line for line in f]
        size = int(size)
        iteration = int(iteration)
        grid = np.array([float(x) for x in grid.split()])
        
        plt.figure(frameon=False)
        plt.imshow(grid.reshape(size, size))

        fname = fname[:-4]

        print fname

        plt.savefig(path_img + fname, transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close()