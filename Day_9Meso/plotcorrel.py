import matplotlib.pyplot as plt
import numpy as np

def plotCorrel(filename, xcol='time', ycol='pe-total', xlabel='Time', ylabel='Energy', title=''):
    headers = open(filename).readline().split()[1:]
    xc = headers.index(xcol)
    yc = headers.index(ycol)
    data = np.loadtxt(filename, skiprows=1)
    fig, ax = plt.subplots()
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.plot(data[:, xc], data[:, yc])
    ax.grid()

