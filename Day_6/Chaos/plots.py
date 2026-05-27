import matplotlib.pyplot as plt
import numpy as np


def plotFile(filename,xcol=1, ycol=2, xlabel='x',ylabel='y',title='title'):

    data = np.loadtxt(filename)
    fig, ax = plt.subplots()

    ax.set(xlabel=xlabel, ylabel=ylabel,
       title=title)

    ax.plot(data[:, xcol-1], data[:, ycol-1])
    ax.grid()

