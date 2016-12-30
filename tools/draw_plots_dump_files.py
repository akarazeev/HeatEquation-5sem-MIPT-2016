#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path_txt = os.path.dirname(__file__) + "/../dump_files/"
path_img = os.path.dirname(__file__) + "/../images/"

# plt.figure(figsize=(5,5))
for fname in os.listdir(path_txt):
    with open(path_txt + fname, 'r') as f:
        size_x, size_y, iteration, grid = [line for line in f]
        size_x = int(size_x)
        size_y = int(size_y)
        iteration = int(iteration)
        grid = np.array([float(x) for x in grid.split()])

        plt.figure(frameon=False)
        fig = plt.imshow(grid.reshape(size_x, size_y))

        fname = fname[:-4]

        print fname

        plt.axis('off')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)

        plt.savefig(path_img + fname, bbox_inches='tight', pad_inches=0)
        plt.close()
