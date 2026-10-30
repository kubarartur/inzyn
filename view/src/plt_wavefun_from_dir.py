import numpy as np

from pathlib import Path
import matplotlib.pyplot as plt
import util.txt_parser as prs
import util.plt_grid as gplt


nm_to_atm = 1/0.0529177
hartree_to_nK = 1.0/3.1668115631*1.0e15

fig, ax = plt.subplots(1, 2)
ax = ax.ravel()

proj_root = Path(__file__).parents[2]
dirname = input("relative path from the project root: ")

print("Reading: " + str(proj_root/dirname/"psi_init.txt"))
grid_init = prs.readTxt(
    proj_root/dirname/"psi_init.txt")

print("Reading: " + str(proj_root/dirname/"V_ext.txt"))
grid_V = prs.readTxt(
    proj_root/dirname/"V_ext.txt")

print("Reading: " + str(proj_root/dirname/"psi_final.txt"))
grid_final = prs.readTxt(
    proj_root/dirname/"psi_final.txt")

gplt.plotN2(grid_init, ax[0])
ax[0].set_title("Initial conditions")
gplt.plotN2(grid_final, ax[1])
ax[1].set_title("Result of minimal energy seek")
plt.show()

fig, ax = plt.subplots(1, 2)
ax = ax.ravel()
gplt.plotN2(grid_V, ax[0])
left, right, top, bottom = grid_V.getExtent((0, 1))
_, _, back, forth = grid_V.getExtent((0, 2))

n = grid_V.vals.shape
print(n)
ax[1].plot(np.linspace(left, right, n[0])/nm_to_atm,
           grid_V.vals[:, n[1]//2, n[2]//2]*hartree_to_nK)
ax[1].plot(np.linspace(top, bottom, n[1])/nm_to_atm,
           grid_V.vals[n[0]//2, :, n[2]//2]*hartree_to_nK)
ax[1].plot(np.linspace(top, bottom, n[2])/nm_to_atm,
           grid_V.vals[n[0]//2, n[1]//2, :]*hartree_to_nK)
plt.title("External potential")
plt.show()

print("Reading: " + str(proj_root/dirname/"E_of_n.txt"))
E_of_psi, extent = prs.readTxt1D(proj_root/dirname/"E_of_n.txt")
extent = [int(x) for x in extent]

n_part = prs.extract_n_part(proj_root/dirname/"input.nml")
print(f"read n_part to be: {n_part}")

skip_first = 0
fig, ax = plt.subplots()
ax.plot(
    np.arange(extent[1] + 1)[skip_first:],
    np.real(E_of_psi)[skip_first:]*hartree_to_nK/n_part,
)
ax.set(
    xlabel="n",
    ylabel="E [nK]"
)
plt.show()
