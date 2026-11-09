import src.read_util as prs
import src.plot_util as mplt
import src.unit_conv as conv


import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import glob

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


for num in sorted(names):
    psi = prs.Spatial_quantity3D(names[num])
    psi.values = np.abs(psi.values)**2
    psi.values *= conv.dens2mum
    psixy = prs.Spatial_quantity2D(psi, (0, 1))
    fig, ax = plt.subplots()
    ax = mplt.plot2D(psixy, ax=ax)
    plt.suptitle(f"step {num}")
    plt.show()
