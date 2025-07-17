import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def plotLBE(filename, title='', solidlines=False):
    # read data from CSV file
    readdata = np.genfromtxt(filename, delimiter=',')
    # remove final column from array due to additional comma at ends of lines
    result = np.delete(readdata, np.s_[-1:], axis=1)
    # find extent of lattice and prepare (x,y) for plots
    row, col = result.shape
    x = np.linspace(0,col,col)
    y = np.linspace(0,row,row)
    X,Y = np.meshgrid(x,y)
    # produce plots of surface and contour
    fig = plt.figure(figsize=[12.0, 4.0])
    fig.suptitle(title)
    ax1 = fig.add_subplot(1, 2, 1, projection='3d') 
    ax1.plot_surface(X, Y, result, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    if solidlines:
        plt.rcParams['contour.negative_linestyle'] = 'solid'
    else:
        plt.rcParams['contour.negative_linestyle'] = 'dashed'
    ax2 = fig.add_subplot(1, 2 ,2)
    cs = ax2.contour(X, Y, result, levels=20, colors='black')
    cf = ax2.contourf(X, Y, result, levels=50, cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax2)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')

def findVortex(filename):
    # read (stream function) data from CSV file, removing final column from array
    # and find extend of lattice to prepare for contour plot
    readdata = np.genfromtxt(filename, delimiter=',')
    result = np.delete(readdata, np.s_[-1:], axis=1)
    row, col = result.shape
    x = np.linspace(0,col,col)
    y = np.linspace(0,row,row)
    X,Y = np.meshgrid(x,y)
    # find location of primary vortex, based on maximum absolute 
    # (i.e. minimum) value of stream function
    ind = np.unravel_index(np.argmin(result, axis=None), result.shape)
    # now generate contour plot and add point for primary vortex centre
    fig, ax = plt.subplots()
    ax.set(xlabel='x', ylabel='y', title='Centre of primary vortex = ({0:d}, {1:d}), '.format(ind[1], ind[0]) + r'$\psi$ = {0:f}'.format(result[ind]))
    plt.rcParams['contour.negative_linestyle'] = 'dashed'
    cs = ax.contour(X, Y, result, levels=20, colors='black')
    cf = ax.contourf(X, Y, result, levels=50, cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax)
    x2 = [ind[1]]
    y2 = [ind[0]]
    ax.scatter(x2, y2, color='black')
    
