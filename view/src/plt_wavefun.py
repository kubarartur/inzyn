import os
import numpy as np

from pathlib import Path
import matplotlib.pyplot as plt
import util.txt_parser as prs
import util.plt_grid as gplt


fig, ax = plt.subplots(1, 2)
ax = ax.ravel()

proj_root = Path(__file__).parents[2]
data_dir = proj_root / "data"
dirs = data_dir.glob("run_*")

numbers = [int(os.path.basename(file).split("_")[1]) for file in dirs]

if len(numbers) != 0:
    curr_num = max(numbers)
else:
    curr_num = 1

print("Reading: " + str(proj_root/f"data/run_{curr_num}/psi_init.txt"))
grid_init = prs.readTxt(
    str(proj_root/f"data/run_{curr_num}/psi_init.txt"))

print("Reading: " + str(proj_root/f"data/run_{curr_num}/psi_final.txt"))
grid_final = prs.readTxt(
    str(proj_root/f"data/run_{curr_num}/psi_final.txt"))

gplt.plotN2(grid_init, ax[0])
ax[0].set_title("Initial conditions")
gplt.plotN2(grid_final, ax[1])
ax[1].set_title("Result of minimal energy seek")
plt.show()

print("Reading: " + str(proj_root/f"data/run_{curr_num}/E_of_n.txt"))
E_of_psi, extent = prs.readTxt1D(proj_root/f"data/run_{curr_num}/E_of_n.txt")
extent = [int(x) for x in extent]

skip_first = 0
fig, ax = plt.subplots()
ax.plot(
    np.arange(extent[1])[skip_first:],
    np.real(E_of_psi)[skip_first:],
)
ax.set(
    xlabel="n",
    ylabel="E in atomic units"
)
plt.show()
