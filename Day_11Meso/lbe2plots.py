#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# read comma-separated input file into arrays

found_density = False
found_rhoN = False
found_K = False
found_psi = False
densitylines = []
rhoNlines = []
Klines = []
psilines = []
with open('density.csv', 'rt') as myfile:
    for myline in myfile:
        if 'density' in myline:
            found_density = True
            found_rhoN = False
            found_K = False
            found_psi = False
        if 'rhoN' in myline:
            found_density = False
            found_rhoN = True
            found_K = False
            found_psi = False
        if 'K' in myline:
            found_density = False
            found_rhoN = False
            found_K = True
            found_psi = False
        if 'psi' in myline:
            found_density = False
            found_rhoN = False
            found_K = False
            found_psi = True
        if found_density:
            densitylines.append(myline.rstrip('\n'))
        if found_rhoN:
            rhoNlines.append(myline.rstrip('\n'))
        if found_K:
            Klines.append(myline.rstrip('\n'))
        if found_psi:
            psilines.append(myline.rstrip('\n'))

# strip out empty entries

densitylines = [i for i in densitylines if i]
rhoNlines = [i for i in rhoNlines if i]
Klines = [i for i in Klines if i]
psilines = [i for i in psilines if i]

# convert string arrays into numerical ones for plotting,
# stripping out last columns due to extra commas at ends of lines

density = np.genfromtxt(densitylines[1:], delimiter=',')
density = np.delete(density, np.s_[-1:], axis=1)
rhoN = np.genfromtxt(rhoNlines[1:], delimiter=',')
rhoN = np.delete(rhoN, np.s_[-1:], axis=1)
Kurv = np.genfromtxt(Klines[1:], delimiter=',')
Kurv = np.delete(Kurv, np.s_[-1:], axis=1)
psi = np.genfromtxt(psilines[1:], delimiter=',')
psi = np.delete(psi, np.s_[-1:], axis=1)

# find extent of lattice and prepare (x,y) for plots
# (should all be identical in size!)

row, col = density.shape

x = np.linspace(0,col,col)
y = np.linspace(0,row,row)

X,Y = np.meshgrid(x,y)

# produce and show surface plot for density

fig1a = plt.figure('Figure 1a: density surface')
ax1a = fig1a.add_subplot(111, projection='3d')
ax1a.plot_surface(X, Y, density, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show contour plot for density

fig1b = plt.figure('Figure 1b: density map')
ax1b = fig1b.add_subplot(111)
cf = ax1b.contourf(X, Y, density, cmap=cm.coolwarm)
fig1b.colorbar(cf, ax=ax1b)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show surface plot for phase index

fig2a = plt.figure('Figure 2a: phase index surface')
ax2a = fig2a.add_subplot(111, projection='3d')
ax2a.plot_surface(X, Y, rhoN, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show contour plot for phase index

fig2b = plt.figure('Figure 2b: phase index map')
ax2b = fig2b.add_subplot(111)
cf = ax2b.contourf(X, Y, rhoN, cmap=cm.coolwarm)
fig2b.colorbar(cf, ax=ax2b)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show surface plot for curvature

fig3a = plt.figure('Figure 3a: curvature surface')
ax3a = fig3a.add_subplot(111, projection='3d')
ax3a.plot_surface(X, Y, Kurv, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show contour plot for curvature

fig3b = plt.figure('Figure 3b: curvature map')
ax3b = fig3b.add_subplot(111)
cf = ax3b.contourf(X, Y, Kurv, cmap=cm.coolwarm)
fig3b.colorbar(cf, ax=ax3b)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# produce and show surface plot for velocity modulus

fig4a = plt.figure('Figure 4a: velocity modulus surface')
ax4a = fig4a.add_subplot(111, projection='3d')
ax4a.plot_surface(X, Y, psi, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel("x")
plt.ylabel("y")
plt.show()


# produce and show contour plot for velocity modulus

fig4b = plt.figure('Figure 4b: velocity modulus surface')
ax4b = fig4b.add_subplot(111)
cf = ax4b.contourf(X, Y, psi, cmap=cm.coolwarm)
fig4b.colorbar(cf, ax=ax4b)
plt.xlabel("x")
plt.ylabel("y")
plt.show()
