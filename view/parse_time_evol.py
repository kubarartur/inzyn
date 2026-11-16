import src.read_util as ru
import src.custom_load as cl
import src.unit_conv as conv

from pathlib import Path
from copy import deepcopy
import sys
import glob
import f90nml
import pickle

proj_root = Path(__file__).parents[1]
dirnames = sys.argv[1:]


runs = []
saved_states = []

times = []
times.append(0)

for dirnum, dirname in enumerate(dirnames[:-1]):
    print("******************************************************************")
    print(f"reading: {dirname}")
    print("******************************************************************")

    dt = f90nml.read(proj_root / dirname /
                     "input.nml")["solver_params"]["CN_real_dt"]
    nstep = f90nml.read(proj_root / dirname /
                        "input.nml")["solver_params"]["real_steps_todo"]
    n_part = f90nml.read(proj_root / dirname /
                         "input.nml")["system_params"]["n_part"]
    grid_dims = f90nml.read(proj_root / dirname /
                            "input.nml")["grid_params"]["n"]
    box_size = f90nml.read(proj_root / dirname /
                           "input.nml")["grid_params"]["box_size"]

    psi = ru.Spatial_quantity3D(proj_root / dirname / "psi_init")

    energies = cl.customLoad(dirname + "E_of_n")
    contact = cl.customLoad(proj_root / dirname / "contact_of_n")
    contact_ext = cl.customLoad(proj_root / dirname / "contact_ext_of_n")
    ext = cl.customLoad(proj_root / dirname / "ext_of_n")
    kinetic = cl.customLoad(proj_root / dirname / "kinetic_of_n")
    dipole = cl.customLoad(proj_root / dirname / "dipole_of_n")

    rright = cl.customLoad(dirname + "r_right_of_n")
    rleft = cl.customLoad(dirname + "r_left_of_n")
    rright = rright/conv.mum_to_atm
    rleft = rleft/conv.mum_to_atm

    rright_max = cl.customLoad(dirname + "r_right_max_of_n")
    rleft_max = cl.customLoad(dirname + "r_left_max_of_n")
    rright_max = rright_max/conv.mum_to_atm
    rleft_max = rleft_max/conv.mum_to_atm

    runs.append({
        "timestep": dt,
        "n_steps": nstep,
        "n_part": n_part,
        "grid_dims": grid_dims,
        "box_size": box_size,
        "r_right": rright,
        "r_left": rleft,
        "r_right_max": rright_max,
        "r_left_max": rleft_max,
        "energies": {
            "total": energies,
            "contact": contact,
            "contact_ext": contact_ext,
            "ext": ext,
            "kinetic": kinetic,
            "dipole": dipole
        },
        "saved_states": {0: deepcopy(psi)}})

    for name in glob.glob(str(proj_root / dirname) + "/real_step*.bin"):
        # remove extension
        name = name.split(".")[-2]

        num = name.split("/")[-1]
        num = num.split("_")[-1]
        num = int(num)

        psi = ru.Spatial_quantity3D(name)
        runs[-1]["saved_states"][int(num)] = deepcopy(psi)

    psi = ru.Spatial_quantity3D(proj_root / dirname / "psi_final")
    runs[-1]["saved_states"][int(nstep)] = deepcopy(psi)

with open(proj_root / dirnames[-1] / "time_evol_runs.pickle", "wb") as f:
    pickle.dump(runs, f)
