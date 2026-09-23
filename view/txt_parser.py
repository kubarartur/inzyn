import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def read_real_txt(fname: 'str'):
    blc = np.loadtxt(fname, skiprows=2, max_rows=1)
    dims = np.loadtxt(fname, skiprows=4, max_rows=1, dtype=np.int32)
    cell_size = np.loadtxt(fname, skiprows=6, max_rows=1)
    values = np.loadtxt(fname, skiprows=8)
    values = values.reshape(dims, order='F')  # ugly
    return blc, cell_size, values


def read_compl_txt(fname: 'str'):
    def s_to_compl128(s):
        s = s.strip("()").split(",")
        return np.complex128(float(s[0]), float(s[1]))

    blc = np.loadtxt(fname, skiprows=2, max_rows=1)
    dims = np.loadtxt(fname, skiprows=4, max_rows=1, dtype=np.int32)
    cell_size = np.loadtxt(fname, skiprows=6, max_rows=1)

    values = np.loadtxt(fname, skiprows=8,
                        dtype=np.complex128, converters=s_to_compl128)

    values = values.reshape(dims, order='F')  # ugly
    return blc, cell_size, values


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


def real_2d(blc, step, values):
    extent = (blc[0],  # left
              blc[0] + values.shape[0]*step[0],  # right
              blc[1],  # bottom
              blc[1] + values.shape[1]*step[1])
    print(f" the extent is {extent}")

    fig, ax = plt.subplots()
    ax.imshow(values.T, extent=extent)
    return fig, ax


if __name__ == "__main__":
    proj_root = Path(__file__).parents[1]
    data_dir = proj_root / "data"
    fname1 = data_dir / "V_ext.txt"
    fname2 = data_dir / "psi.txt"

    blc, step, vals = read_real_txt(fname1)
    fig, ax = plot_real(blc, step, vals)
    plt.show()

    blc, step, vals = read_compl_txt(fname2)
    zc = vals.shape[2]//2
    data2D = abs(vals[:, :, zc])  # obnoxious incontinuous slice
    fig, ax = real_2d(blc[0:2], step[0:2], data2D)
    plt.show()
