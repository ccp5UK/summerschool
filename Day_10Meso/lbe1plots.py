#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# part 1: density

# read comma-separated input files into arrays

readdata1 = np.genfromtxt('data_density.csv', delimiter=',')

# remove last column from array (due to extra comma at ends of lines)

result1 = np.delete(readdata1, np.s_[-1:], axis=1)

# find extent of lattice and prepare (x,y) for plots

row, col = result1.shape

x = np.linspace(0,col,col)
y = np.linspace(0,row,row)

X,Y = np.meshgrid(x,y)

# produce and show surface plot

fig1 = plt.figure('Figure 1: density surface plot')
ax1 = fig1.add_subplot(111, projection='3d')
ax1.plot_surface(X, Y, result1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show contour plot

fig2 = plt.figure('Figure 2: density contour plot')
ax2 = fig2.add_subplot(111)
cs = ax2.contour(X, Y, result1, levels=20, colors='black')
cf = ax2.contourf(X, Y, result1, levels=50, cmap=cm.coolwarm)
fig2.colorbar(cf, ax=ax2)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# part 2: stream functions

# read comma-separated input files into arrays

readdata2 = np.genfromtxt('data_stfn.csv', delimiter=',')

# remove last column from array (due to extra comma at ends of lines)

result2 = np.delete(readdata2, np.s_[-1:], axis=1)

# find extent of lattice and prepare (x,y) for plots

row, col = result2.shape

x = np.linspace(0,col,col)
y = np.linspace(0,row,row)

X,Y = np.meshgrid(x,y)

# produce and show surface plot

fig3 = plt.figure('Figure 3: stream function surface plot')
ax3 = fig3.add_subplot(111, projection='3d')
ax3.plot_surface(X, Y, result2, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show contour plot

plt.rcParams['contour.negative_linestyle'] = 'solid'
fig4 = plt.figure('Figure 4: streamlines (stream function contour plot)')
ax4 = fig4.add_subplot(111)
cs = ax4.contour(X, Y, result2, levels=20, colors='black')
cf = ax4.contourf(X, Y, result2, levels=50, cmap=cm.coolwarm)
fig4.colorbar(cf, ax=ax4)
plt.xlabel("x")
plt.ylabel("y")
plt.show()
