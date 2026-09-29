from pathlib import Path
import matplotlib.pyplot as plt
import util.txt_parser as prs
import util.plt_grid as gplt

fig, ax = plt.subplots(1, 2)
ax = ax.ravel()

proj_root = Path(__file__).parents[2]

grid_init = prs.readTxt(str(proj_root/"data/psi_init.txt"))
grid_final = prs.readTxt(str(proj_root/"data/psi_final.txt"))

gplt.plotN2(grid_init, ax[0])
ax[0].set_title("Initial conditions")
gplt.plotN2(grid_final, ax[1])
ax[1].set_title("Result of minimal energy seek")
plt.show()
