import util.txt_parser as prs
import matplotlib.pyplot as plt
import numpy as np


def plotN2(grid: prs.Grid, axes=None):
    if axes is None:
        fig, axes = plt.subplots()

    n = np.real(grid.vals*np.conjugate(grid.vals))
    ih = int(n.shape[2]/2)
    axes.imshow(n[:, :, ih].T, extent=grid.getExtent((0, 1)))
