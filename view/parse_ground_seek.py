import src.read_util as ru
import src.custom_load as cl

from pathlib import Path
import sys
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
                     "input.nml")["solver_params"]["Euler_imag_dt"]
    nstep = f90nml.read(proj_root / dirname /
                        "input.nml")["solver_params"]["imag_steps_todo"]
    n_part = f90nml.read(proj_root / dirname /
                         "input.nml")["system_params"]["n_part"]
    grid_dims = f90nml.read(proj_root / dirname /
                            "input.nml")["system_params"]["n_part"]
    box_size = f90nml.read(proj_root / dirname /
                           "input.nml")["grid_params"]["box_size"]

    energies = cl.customLoad(proj_root / dirname / "E_of_n")
    contact = cl.customLoad(proj_root / dirname / "contact_of_n")
    contact_ext = cl.customLoad(proj_root / dirname / "contact_ext_of_n")
    ext = cl.customLoad(proj_root / dirname / "ext_of_n")
    kinetic = cl.customLoad(proj_root / dirname / "kinetic_of_n")
    dipole = cl.customLoad(proj_root / dirname / "dipole_of_n")

    delta_ns = cl.customLoad(proj_root / dirname / "delta_n_of_n")
    psi_init = ru.Spatial_quantity3D(proj_root / dirname / "psi_init")
    psi_final = ru.Spatial_quantity3D(proj_root / dirname / "psi_final")
    V_ext = ru.Spatial_quantity3D(proj_root / dirname / "V_ext")

    runs.append({
        "timestep": dt,
        "n_steps": nstep,
        "n_part": n_part,
        "grid_dims": grid_dims,
        "box_size": box_size,
        "energies": {
            "total": energies,
            "contact": contact,
            "contact_ext": contact,
            "ext": ext,
            "kinetic": kinetic,
            "diople": dipole
        },
        "delta_ns": delta_ns,
        "V_ext": V_ext,
        "psi_init": psi_init,
        "psi_final": psi_final})

with open(proj_root / dirnames[-1] / "ground_seek_runs.pickle", "wb") as f:
    pickle.dump(runs, f)
