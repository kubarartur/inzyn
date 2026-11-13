import src.read_util as prs
import src.plot_util as mplt
import src.unit_conv as conv

import imageio.v3 as iio
import matplotlib.pyplot as plt
import numpy as np
import glob

from pathlib import Path


proj_root = Path(__file__).parents[1]
dirname = input("relative path from the project root: ")

fname = proj_root / dirname / "psi_init"
names = {0: fname}

for name in glob.glob(str(proj_root / dirname) + "/real_step*.bin"):
    # remove extension
    name = name.split(".")[-2]

    num = name.split("/")[-1]
    num = num.split("_")[-1]
    num = int(num)

    names[num] = name

images = []

for num in sorted(names):

    psi = prs.Spatial_quantity3D(names[num])
    psi.values = np.abs(psi.values)**2
    psi.values *= conv.dens2mum
    obj = prs.Spatial_quantity2D(psi, (0, 1))

    fig, ax = plt.subplots()

    # 1. Call your existing function untouched
    ax = mplt.plot2D(obj, ax=ax)

    # 2. Capture the plot into memory as an image array
    fig.canvas.draw()
    images.append(np.asarray(fig.canvas.buffer_rgba()))

    # 3. Clean up the figure to save memory
    plt.close(fig)

# Save everything as a GIF (duration is in milliseconds per frame)
iio.imwrite("system_evolutionxy.gif", images, duration=200, loop=0)

images = []
for num in sorted(names):

    psi = prs.Spatial_quantity3D(names[num])
    psi.values = np.abs(psi.values)**2
    psi.values *= conv.dens2mum
    obj = prs.Spatial_quantity2D(psi, (0, 2))

    fig, ax = plt.subplots()

    # 1. Call your existing function untouched
    ax = mplt.plot2D(obj, ax=ax)

    # 2. Capture the plot into memory as an image array
    fig.canvas.draw()
    images.append(np.asarray(fig.canvas.buffer_rgba()))

    # 3. Clean up the figure to save memory
    plt.close(fig)

# Save everything as a GIF (duration is in milliseconds per frame)
iio.imwrite("system_evolutionxz.gif", images, duration=200, loop=0)
