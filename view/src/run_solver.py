import subprocess
import os

from pathlib import Path

proj_root = Path(__file__).parents[2]
data_dir = proj_root / "data"
dirs = data_dir.glob("run_*")

numbers = [int(os.path.basename(file).split("_")[1]) for file in dirs]

if len(numbers) != 0:
    curr_num = max(numbers) + 1
else:
    curr_num = 1
subprocess.run(["mkdir",
               str(proj_root) + f"/data/run_{curr_num}"])

subprocess.run([str(proj_root) + "/build/src/main",
               str(proj_root) + f"/data/run_{curr_num}/"])
