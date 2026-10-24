import util.txt_parser as prs
import matplotlib.pyplot as plt
import numpy as np

nm_to_atm = 1/0.0529177


def plotN2(grid: prs.Grid, axes=None):
    if axes is None:
        fig, axes = plt.subplots()

    n = np.real(grid.vals*np.conjugate(grid.vals))
    ih = int(n.shape[2]/2)
    axes.imshow(n[:, :, ih].T,
                extent=[x/nm_to_atm for x in grid.getExtent((0, 1))])
