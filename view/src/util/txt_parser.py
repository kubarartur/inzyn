import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import re


class Grid:
    def __init__(self,
                 blc: npt.NDArray[npt.Shape["3"], npt.Float64],
                 dims: npt.NDArray[npt.Shape["3"], npt.Int32],
                 cell_size: npt.NDArray[npt.Shape["3"], npt.Float64],
                 vals: npt.NDArray[npt.Shape["3", "3", "3"]],
                 descr=""):

        self.blc = blc
        self.grid_size = dims
        self.cell_size = cell_size
        self.vals = vals
        self.descr = descr

    def getExtent(self, dim: tuple[int, int]) -> tuple[float, ...]:
        """
        Args:
            dimensions: contains two integers indicating dimensions
            0, 2, 3 -> x, y, z
            rest -> ValueError

        Returns:
            Touple of ranges in each of the specified dimension.
            Compatible with matplotlibs (left, right, bottom, top)
        """
        for di in dim:
            if di not in {0, 1, 2}:
                raise ValueError(
                    f"wrong dimension expected 0, 1 or 2. got: {dim}")

        return (
            self.blc[dim[0]],
            self.blc[dim[0]] + self.cell_size[dim[0]]*self.grid_size[dim[0]],
            self.blc[dim[1]],
            self.blc[dim[1]] + self.cell_size[dim[1]]*self.grid_size[dim[1]]
        )

    def __str__(self):
        return ("parser<Grid>: " + self.descr)


def readTxt(fname: 'str') -> tuple[Grid]:
    blfc = np.loadtxt(fname, skiprows=2, max_rows=1, dtype=np.float64)
    dims = np.loadtxt(fname, skiprows=4, max_rows=1, dtype=np.int32)
    cell_size = np.loadtxt(fname, skiprows=6, max_rows=1)
    with open(fname) as f:
        stri = f.readline()
        stri = stri.split(":")
        descr = ""
        for s in stri[1:]:
            descr += s

    try:
        values = np.loadtxt(fname, skiprows=8, dtype=np.float64)
    except ValueError:
        def s_to_compl128(s):
            s = s.strip("()").split(",")
            return np.complex128(float(s[0]), float(s[1]))

        values = np.loadtxt(fname, skiprows=8,
                            dtype=np.complex128, converters=s_to_compl128)

    print(f"initial shape of values is: {values.shape}")
    if len(values.shape) > 1:
        values = values.reshape(-1)
    print(f"after flatten shape of values is: {values.shape}")
    values = values.reshape(dims, order='F')

    # TODO: make this parse description from the text dump from fortran
    grid = Grid(blfc, dims, cell_size, values, descr=descr)
    return grid


def readTxt1D(fname: str):
    extent = np.loadtxt(fname, skiprows=2, max_rows=1, dtype=np.float64)
    dims = np.loadtxt(fname, skiprows=4, max_rows=1, dtype=np.int32)
    with open(fname) as f:
        stri = f.readline()
        stri = stri.split(":")
        descr = ""
        for s in stri[1:]:
            descr += s

    try:
        values = np.loadtxt(fname, skiprows=8, dtype=np.float64)
    except ValueError:
        def s_to_compl128(s):
            s = s.strip("()").split(",")
            return np.complex128(float(s[0]), float(s[1]))

        values = np.loadtxt(fname, skiprows=8,
                            dtype=np.complex128, converters=s_to_compl128)

    assert (values.shape == dims)

    return values, extent


def plot_real(blc, step, vals, linestyle="-"):
    print(f"shape is {vals.shape}")
    # indices of centers
    xc = vals.shape[0]//2
    yc = vals.shape[1]//2
    zc = vals.shape[2]//2

    print("centers are: ", xc, yc, zc)

    print(f"step is {step} ")
    print(f"blc is {blc} ")

    print(" this should be definitely be zeros ------")
    print(blc[0] + 0.5*step[0]*(vals.shape[0] + 1))
    print(blc[1] + 0.5*step[1]*(vals.shape[1] + 1))
    print(blc[2] + 0.5*step[2]*(vals.shape[2] + 1))

    print(" this should be zeros ------")
    print(blc[0] + step[0]*(xc+1))
    print(blc[1] + step[1]*(yc+1))
    print(blc[2] + step[2]*(zc+1))

    fig, ax = plt.subplots()
    x_coords = blc[0] + (np.arange(vals.shape[0])+1)*step[0]
    ax.plot(x_coords, vals[:, yc, zc], "r"+linestyle)

    y_coords = blc[1] + (np.arange(vals.shape[1]) + 1)*step[1]
    ax.plot(y_coords, vals[xc, :, zc], "g"+linestyle)

    z_coords = blc[2] + (np.arange(vals.shape[2]) + 1)*step[2]
    ax.plot(z_coords, vals[xc, yc, :], "b"+linestyle)
    plt.legend(["x", "y", "z"])

    return fig, ax


def extract_n_part(file_path):
    pattern = re.compile(r"n_part\s*=\s*(\d+)")

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                return int(match.group(1))

    return None
