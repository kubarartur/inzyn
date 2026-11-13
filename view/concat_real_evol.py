import src.read_util as prs
import src.custom_load as cl

import numpy as np
from pathlib import Path
import glob
import f90nml
import sys
import subprocess as sb

proj_root = Path(__file__).parents[1]
dirnames = sys.argv[1:]


fname = proj_root / dirnames[0] / "psi_init"
names = {0: fname}

curr_tot = 0
for dirname in dirnames[:-1]:
    maxnum = 0
    print("***************************************************************************")
    print(f"reading: {dirname}, curr tot is: {curr_tot}")
    print("***************************************************************************")
    for name in glob.glob(str(proj_root / dirname) + "/real_step*.bin"):
        # remove extension
        name = name.split(".")[-2]

        num = name.split("/")[-1]
        num = num.split("_")[-1]
        num = int(num)

        if num > maxnum:
            maxnum = num

        names[num + curr_tot] = name

    curr_tot += maxnum

for namenum in sorted(names)[1:]:
    print("-------------------------------------")
    for subname in glob.glob(str(names[namenum]) + "*"):
        pthobj = Path(subname)

        print(f"target: {pthobj}")

        split_name = pthobj.name.split("_")
        split_name[-1], extension = split_name[-1].split(".")
        split_name.append(extension)
        split_name[2] = str(namenum).zfill(6)

        to_name = '_'.join(split_name)
        print(f"copied to: {proj_root / dirnames[-1] / to_name}")
        sb.run(["cp", str(pthobj), str(proj_root / dirnames[-1] / to_name)])

