import numpy as np
import matplotlib.pyplot as plt
import src.read_util as prs


def plot2D(quant2D: prs.Spatial_quantity2D,
           ax=None):

    if ax is None:
        _, ax = plt.subplots()

    left = quant2D.blc[0] + quant2D.dr[0]
    right = quant2D.blc[0] + quant2D.dr[0]*quant2D.values.shape[0]
    bottom = quant2D.blc[1] + quant2D.dr[1]
    top = quant2D.blc[1] + quant2D.dr[1]*quant2D.values.shape[1]

    X = np.linspace(left, right, quant2D.values.shape[0])
    X = np.stack([X for _ in range(quant2D.values.shape[1])], axis=0)
    Y = np.linspace(bottom, top, quant2D.values.shape[1])
    Y = np.stack([Y for _ in range(quant2D.values.shape[0])], axis=1)

    mesh = ax.pcolormesh(X, Y, quant2D.values.T, shading='nearest')
    plt.colorbar(mesh, ax=ax)

    return ax


def plotCross(quant1D: prs.Spatial_quantity1D,
              ax=None):

    if ax is None:
        _, ax = plt.subplots()

    left = quant1D.blc[0] + quant1D.dr[0]
    right = quant1D.blc[0] + quant1D.dr[0]*quant1D.values.shape[0]

    xs = np.linspace(left, right, quant1D.values.shape[0])
    ax.plot(xs, quant1D.values)

    return ax
