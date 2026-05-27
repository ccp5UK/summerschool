import matplotlib.pyplot as plt
import numpy as np

def plotStats(filename, xcol=1, ycol=2, xlabel='Time', ylabel='Energy', title=''):
    plt.rcParams['figure.figsize'] = [6.4, 4.8]
    data = np.loadtxt(filename)
    fig, ax = plt.subplots()
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.plot(data[:, xcol-1], data[:, ycol-1])
    ax.grid()

