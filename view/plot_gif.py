import src.read_util as prs
import src.plot_util as mplt
import src.unit_conv as conv

import imageio.v3 as iio
import matplotlib.pyplot as plt
import numpy as np
import sys
import pickle

from pathlib import Path


proj_root = Path(__file__).parents[1]
fname = sys.argv[1]

with open(fname, "rb") as f:
    runs = pickle.load(f)

images = []
for run in runs:
    for state_num in run["saved_states"]:
        state = run["saved_states"][state_num]

        state.values = np.abs(state.values)**2
        state.values *= conv.dens2mum
        xy = prs.Spatial_quantity2D(state, (0, 1))

        fig, ax = plt.subplots()
        ax = mplt.plot2D(xy, ax=ax)
        ax.plot(run["r_right"][0, state_num] - state.dr[0],
                run["r_right"][1, state_num], 'ko')

        ax.plot(run["r_left"][0, state_num],
                run["r_left"][1, state_num], 'ko')

        fig.canvas.draw()
        images.append(np.asarray(fig.canvas.buffer_rgba()))
        plt.close(fig)

iio.imwrite("system_evolutionxy.gif", images, duration=800, loop=0)

images = []
for run in runs:
    for state_num in run["saved_states"]:
        state = run["saved_states"][state_num]

        state.values = np.abs(state.values)**2
        state.values *= conv.dens2mum
        xy = prs.Spatial_quantity2D(state, (0, 2))

        fig, ax = plt.subplots()
        ax = mplt.plot2D(xy, ax=ax)
        fig.canvas.draw()
        images.append(np.asarray(fig.canvas.buffer_rgba()))
        plt.close(fig)

iio.imwrite("system_evolutionxz.gif", images, duration=200, loop=0)
