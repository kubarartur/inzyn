import src.read_util as prs
import src.plot_util as mplt
import src.unit_conv as conv
import glob

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

proj_root = Path(__file__).parents[1]
dirname = input("relative path from the project root: ")


def plot_diff(psi1, psi2, title=""):

    xy1_cross = prs.Spatial_quantity2D(psi1, (0, 1))
    xy1_cross.values = np.abs(xy1_cross.values)**2*conv.dens2mum
    xy2_cross = prs.Spatial_quantity2D(psi2, (0, 1))
    xy2_cross.values = np.abs(xy2_cross.values)**2*conv.dens2mum

    xz1_cross = prs.Spatial_quantity2D(psi1, (0, 2))
    xz1_cross.values = np.abs(xz1_cross.values)**2*conv.dens2mum
    xz2_cross = prs.Spatial_quantity2D(psi2, (0, 2))
    xz2_cross.values = np.abs(xz2_cross.values)**2*conv.dens2mum

    x1_cross = prs.Spatial_quantity1D(psi1, (0))
    x1_cross.values = np.abs(x1_cross.values)**2*conv.dens2mum
    x2_cross = prs.Spatial_quantity1D(psi2, (0))
    x2_cross.values = np.abs(x2_cross.values)**2*conv.dens2mum

    xy2_cross.values -= xy1_cross.values
    xz2_cross.values -= xz1_cross.values
    x2_cross.values -= x1_cross.values

    fig, ax = plt.subplots(1, 2)
    ax = ax.ravel()
    mplt.plot2D(xy2_cross, ax=ax[0])
    mplt.plot2D(xz2_cross, ax=ax[1])
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")
    ax[1].set_xlabel("x")
    ax[1].set_ylabel("z")
    plt.suptitle(title)
    plt.show()

    fig, ax = plt.subplots()
    mplt.plotCross(x2_cross, ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("rho")
    plt.suptitle(title)
    plt.show()


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

num = sorted(names)[-1]

psii = prs.Spatial_quantity3D(names[0])
psif = prs.Spatial_quantity3D(names[num])

plot_diff(psii, psif, title=f'step: {num} - init')
