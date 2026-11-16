import src.read_util as prs
import src.plot_util as mplt
import src.unit_conv as conv
from time import sleep


import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys
import pickle

proj_root = Path(__file__).parents[1]
fname = sys.argv[1]


def plot_psi(psi, title=""):
    xy_cross = prs.Spatial_quantity2D(psi, (0, 1))
    xy_cross.values = np.abs(xy_cross.values)**2*conv.dens2mum

    xz_cross = prs.Spatial_quantity2D(psi, (0, 2))
    xz_cross.values = np.abs(xz_cross.values)**2*conv.dens2mum

    x_cross = prs.Spatial_quantity1D(psi, (0))
    x_cross.values = np.abs(x_cross.values)**2*conv.dens2mum

    fig, ax = plt.subplots(1, 2)
    ax = ax.ravel()
    mplt.plot2D(xy_cross, ax=ax[0])
    mplt.plot2D(xz_cross, ax=ax[1])
    ax[0].set_xlabel(r"x $[\mu m]$")
    ax[0].set_ylabel(r"y $[\mu m]$")
    ax[1].set_xlabel(r"x $[\mu m]$")
    ax[1].set_ylabel(r"z $[\mu m]$")
    plt.suptitle(title)
    plt.show()

    fig, ax = plt.subplots()
    mplt.plotCross(x_cross, ax=ax)
    ax.set_xlabel(r"x $[\mu m]$")
    ax.set_ylabel(r"rho $[1/\mu m^3]$")
    plt.suptitle(title)
    plt.show()


def plotAllCrosses(val: prs.Spatial_quantity3D, title="", ax=None,
                   addlabel=""):
    if ax is None:
        fig, ax = plt.subplots()

    x_cross = prs.Spatial_quantity1D(val, 0)
    y_cross = prs.Spatial_quantity1D(val, 1)
    z_cross = prs.Spatial_quantity1D(val, 2)

    mplt.plotCross(x_cross, ax=ax, label="x " + addlabel)
    mplt.plotCross(y_cross, ax=ax, label="y " + addlabel)
    mplt.plotCross(z_cross, ax=ax, label="z " + addlabel)

    plt.legend(["x " + addlabel, "y " + addlabel, "z " + addlabel])
    plt.suptitle(title)


with open(fname, "rb") as f:
    runs = pickle.load(f)

psi_init = runs[0]["psi_init"]
psi_final = runs[-1]["psi_final"]

V_exts = [run["V_ext"] for run in runs]

for V_ext in V_exts:
    V_ext.values = V_ext.values*conv.dens2mum

plot_psi(psi_final, "result")
plot_psi(psi_init, "initial conditions")

fig, ax = plt.subplots()
for i, V_ext in enumerate(V_exts):
    plotAllCrosses(V_ext, title="External potential",
                   ax=ax, addlabel=f"run {i}")
plt.show()

Es = [run["energies"]["total"]*conv.hartree_to_nK for run in runs]
ns = [run["delta_ns"]*conv.dens2mum for run in runs]

n_parts = [run["n_part"] for run in runs]

fig, ax = plt.subplots(1, 2)
ax = ax.ravel()
last_n1 = 0
last_n2 = 0

for i, (E, n, n_part) in enumerate(zip(Es, ns, n_parts)):

    n1 = np.arange(last_n1, last_n1 + E.shape[0])
    n2 = np.arange(last_n2, last_n2 + n.shape[0])
    last_n1 = n1[-1]
    last_n2 = n2[-1]

    ax[0].plot(n1, np.real(E)/n_part, label=f"run {i}")
    ax[1].plot(n2, n, label=f"run {i}")

    ax[0].axvline(x=last_n1)
    ax[1].axvline(x=last_n1)
    plt.legend()

ax[0].set_xlabel("n")
ax[1].set_xlabel("n")
ax[0].set_ylabel("E/n_part [nk]")
ax[1].set_ylabel(r"$\Delta n$")
plt.show()
