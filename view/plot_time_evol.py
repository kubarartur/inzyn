import src.read_util as prs
import src.plot_util as mplt
import src.unit_conv as conv
import src.custom_load as cl


import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import glob
import f90nml

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

dt = f90nml.read(proj_root / dirname /
                 "input.nml")["solver_params"]["CN_real_dt"]
dt *= conv.jat_to_ms

r_right = cl.customLoad(proj_root / dirname / "r_right_of_n")
r_left = cl.customLoad(proj_root / dirname / "r_left_of_n")

r_right = r_right/conv.mum_to_atm
r_left = r_left/conv.mum_to_atm
ns = np.arange(r_right.shape[1])*dt

fig, ax = plt.subplots()
ax.plot(ns, r_right[0, :], label="r_right x")
ax.plot(ns, r_left[0, :], label="r_left x")

ax.plot(ns, r_right[1, :], label="r_right y")
ax.plot(ns, r_left[1, :], label="r_left y")
plt.legend()
plt.show()

ncz = f90nml.read(proj_root / dirname / "input.nml")["system_params"]["n_part"]
Es = cl.customLoad(proj_root / dirname / "E_n_of_n")
Es = np.real(Es)*conv.hartree_to_nK/ncz
ns = np.arange(Es.shape[0])

fig, ax = plt.subplots()
ax.plot(ns, Es, label="E of N")
ax.set_xlabel("n")
ax.set_ylabel("E [nK]")
plt.show()

for num in sorted(names):
    psi = prs.Spatial_quantity3D(names[num])
    psi.values = np.abs(psi.values)**2
    psi.values *= conv.dens2mum
    psixy = prs.Spatial_quantity2D(psi, (0, 1))
    fig, ax = plt.subplots()
    ax = mplt.plot2D(psixy, ax=ax)
    plt.suptitle(f"step {num}")
    plt.show()
