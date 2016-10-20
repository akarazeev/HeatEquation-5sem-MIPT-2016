#!/usr/local/bin/python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path = "dump_files/"

# plt.figure(figsize=(5,5))
for fname in os.listdir(path):
    with open(path + fname, 'r') as f:
        size, iteration, grid = [line for line in f]
        size = int(size)
        iteration = int(iteration)
        grid = np.array([float(x) for x in grid.split()])
        
        plt.figure(frameon=False)
        plt.imshow(grid.reshape(size, size))
        fname = fname[:-4]
        print 'images/' + fname + '.jpg'
        plt.savefig('images/' + fname, transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close()