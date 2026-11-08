import src.read_util as prs
import src.plot_util as mplt
import src.unit_conv as conv

import f90nml

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

proj_root = Path(__file__).parents[1]
dirname = input("relative path from the project root: ")


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
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")
    ax[1].set_xlabel("x")
    ax[1].set_ylabel("z")
    plt.suptitle(title)
    plt.show()

    fig, ax = plt.subplots()
    mplt.plotCross(x_cross, ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("rho")
    plt.suptitle(title)
    plt.show()


def plotAllCrosses(val: prs.Spatial_quantity3D, title=""):
    val.values = pot.values*conv.hartree_to_nK
    x_cross = prs.Spatial_quantity1D(val, 0)
    y_cross = prs.Spatial_quantity1D(val, 1)
    z_cross = prs.Spatial_quantity1D(val, 2)

    fig, ax = plt.subplots()
    mplt.plotCross(x_cross, ax)
    mplt.plotCross(y_cross, ax)
    mplt.plotCross(z_cross, ax)
    plt.legend(["x", "y", "z"])
    plt.suptitle(title)
    plt.show()


fname = proj_root / dirname / "psi_final"
psi = prs.Spatial_quantity3D(fname)
plot_psi(psi, "result")

fname = proj_root / dirname / "psi_init"
psi = prs.Spatial_quantity3D(fname)
plot_psi(psi, "initial conditions")

fname = proj_root / dirname / "V_ext"
pot = prs.Spatial_quantity3D(fname)
plotAllCrosses(pot, title="External potential")


fname = proj_root / dirname / "E_of_n.npy"
E = np.load(fname)*conv.hartree_to_nK

fname = proj_root / dirname / "delta_n_of_n.npy"
ns = np.load(fname)*conv.dens2mum

fname = proj_root / dirname / "input.nml"
n_part = f90nml.read(fname)["system_params"]["n_part"]
print(f"n_part is: {n_part}")

n1 = np.arange(E.shape[0])
n2 = np.arange(ns.shape[0])
fig, ax = plt.subplots(1, 2)
ax = ax.ravel()
ax[0].plot(n1, np.real(E)/n_part)
ax[1].plot(n2, ns)
ax[0].set_xlabel("n")
ax[1].set_xlabel("n")
ax[0].set_ylabel("E/n_part [nk]")
ax[1].set_ylabel(r"$\Delta n$")
plt.show()
