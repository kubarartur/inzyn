import src.unit_conv as conv

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys
import pickle

proj_root = Path(__file__).parents[1]
fname = sys.argv[1]

with open(fname, "rb") as f:
    runs = pickle.load(f)

print(runs)
fig, ax = plt.subplots()
n_curr = 0
for i, run in enumerate(runs):
    r_right = run["r_right"]
    r_left = run["r_left"]

    n = np.arange(n_curr, n_curr + r_right.shape[1])
    times = n*run["timestep"]*conv.jat_to_ms
    n_curr += r_right.shape[1]

    ax.plot(times, r_right[0, :], label=f"r_right, run: {i}")
    ax.plot(times, r_left[0, :], label=f"r_left, run: {i}")
    ax.axvline(x=times[-1])

plt.legend()
ax.set_xlabel("t [ms]")
ax.set_ylabel(r"r [$\mu$s]")
plt.show()

fig, ax = plt.subplots()
n_curr = 0
for i, run in enumerate(runs):
    ncz = run["n_part"]
    E = np.real(run["energies"]["total"])/ncz*conv.hartree_to_nK
    contact = np.real(run["energies"]["contact"])/ncz*conv.hartree_to_nK
    contact_ext = np.real(
        run["energies"]["contact_ext"])/ncz*conv.hartree_to_nK
    kinetic = np.real(run["energies"]["kinetic"])/ncz*conv.hartree_to_nK
    dipole = np.real(run["energies"]["dipole"])/ncz*conv.hartree_to_nK
    ext = np.real(run["energies"]["ext"])/ncz*conv.hartree_to_nK

    n = np.arange(n_curr, n_curr + E.shape[0])
    n_curr += E.shape[0]
    times = n*run["timestep"]*conv.jat_to_ms
    n_curr += r_right.shape[0]

    ax.plot(times, E, label=f"E in time, run {i}")
    ax.plot(times, contact + contact_ext,
            label=f"tot contact in time, run {i}")
    ax.plot(times, kinetic, label=f"kinetic in time, run {i}")
    ax.plot(times, dipole, label=f"dipole in time, run {i}")
    ax.plot(times, ext, label=f"external in time, run {i}")
    ax.plot(times, E - contact - contact_ext - kinetic - dipole - ext,
            label=f"total - partials, run {i}")

    ax.axvline(x=times[-1])

plt.legend()
ax.set_xlabel("t [ms]")
ax.set_ylabel("E [nK]")
plt.show()


# for num in sorted(names):
# psi = prs.Spatial_quantity3D(names[num])
# psi.values = np.abs(psi.values)**2
# psi.values *= conv.dens2mum
# psixy = prs.Spatial_quantity2D(psi, (0, 1))
# fig, ax = plt.subplots()
# ax = mplt.plot2D(psixy, ax=ax)
# plt.suptitle(f"step {num}")
# plt.show()
