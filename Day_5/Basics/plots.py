import matplotlib.pyplot as plt
import numpy as np


def plotFile(filename,xlabel='x',ylabel='y',title='title'):

    data = np.loadtxt(filename)
    fig, ax = plt.subplots()

    ax.set(xlabel=xlabel, ylabel=ylabel,
       title=title)

    ax.plot(data[:, 0], data[:, 1])
    ax.grid()

def plotFileErr(filename,xlabel='x',ylabel='y',title='title'):

    data = np.loadtxt(filename)
    fig, ax = plt.subplots()

    ax.set(xlabel=xlabel, ylabel=ylabel,
       title=title)

    ax.errorbar(data[:, 0], data[:, 1],yerr=data[:,2])
    ax.grid()

def plotFileLog(filename,xlabel='x',ylabel='y',title='title'):

    data = np.loadtxt(filename)
    fig, ax = plt.subplots()

    ax.set(xlabel=xlabel, ylabel=ylabel,
       title=title)

    ax.loglog(data[:, 0], data[:, 1])
    ax.grid()

