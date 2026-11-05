import util.txt_parser as prs
import matplotlib.pyplot as plt
import numpy as np

nm_to_atm = 1/0.0529177


def plot_cross(grid: prs.Grid, axes=None, dim=0):
    if axes is None:
        fig, axes = plt.subplots()

    n = np.abs(grid.vals[
        :,
        grid.vals.shape[1]//2,
        grid.vals.shape[2]//2
    ])**2

    left, right, _, _ = grid.getExtent((0, 1))
    xs = np.linspace(left, right, n.shape[0])

    axes.plot(xs, n*1e9*nm_to_atm**3)

    return axes


def plotN2(grid: prs.Grid, axes=None, dims=(0, 1)):
    if axes is None:
        fig, axes = plt.subplots()

    n = np.abs(grid.vals)**2

    dimset = set(dims)
    for dim in range(3):
        if dim not in dimset:
            chosen_dim = dim
            inx = int(n.shape[dim]/2)

    if chosen_dim == 2:
        axes.imshow(n[:, :, inx].T,
                    extent=[x/nm_to_atm for x in grid.getExtent(dims)])
    elif chosen_dim == 1:
        axes.imshow(n[:, inx, :].T,
                    extent=[x/nm_to_atm for x in grid.getExtent(dims)])
    elif chosen_dim == 0:
        axes.imshow(n[inx, :, :].T,
                    extent=[x/nm_to_atm for x in grid.getExtent(dims)])

    dV = grid.cell_size[0]*grid.cell_size[1]*grid.cell_size[2]

    rho_max = n.max()*1e9*nm_to_atm**3

    heights = find_heights(n, 1e-9/nm_to_atm**3, grid.cell_size[2]/nm_to_atm)
    max_height = np.max(heights)

    left, right, bottom, top = grid.getExtent((1, 2))
    back, fwd, _, _ = grid.getExtent((0, 1))
    Vol = np.abs(left - right)*np.abs(top - bottom)*np.abs(fwd - back)
    dV = Vol/(n.shape[0]*n.shape[1]*n.shape[2])

    n_part = np.sum(n*dV)

    return rho_max, max_height, n_part


def find_heights(n, tol, drz):
    heights = np.empty((
        n.shape[0],
        n.shape[1]))

    for i in range(n.shape[0]):
        for j in range(n.shape[1]):
            for k in range(n.shape[2]//2):
                if n[i, j, k] > tol:
                    start = k
            for k in range(n.shape[2]//2, n.shape[2], 1):
                if n[i, j, k] < tol:
                    stop = k

            if (n[i, j, n.shape[2]//2] < tol):
                heights[i, j] = 0
            else:
                try:
                    heights[i, j] = (stop - start)*drz
                except UnboundLocalError:
                    heights[i,j] = 0

    return heights


def int_rhos(n, dV):
    rhos = np.zeros((
        n.shape[0],
        n.shape[1]))

    for i in range(n.shape[0]):
        for j in range(n.shape[1]):
            for k in range(n.shape[2]):
                rhos[i, j] += n[i, j, k]*dV

    return rhos
