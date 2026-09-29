from pathlib import Path
import util.txt_parser as prs

proj_root = Path(__file__).parents[2]
data_dir = proj_root / "data"
fname1 = data_dir / "psi_final.txt"
fname2 = data_dir / "psi_init.txt"

grid = prs.readTxt(fname1)
print(grid)
