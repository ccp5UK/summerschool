import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import skimage as ski

def plotMCLB(filename, title=''):
    # read all data from CSV file
    found_density = False
    found_rhoN = False
    found_K = False
    found_psi = False
    densitylines = []
    rhoNlines = []
    Klines = []
    psilines = []
    with open(filename, 'rt') as myfile:
        for myline in myfile:
            if 'density' in myline:
                found_density = True
                found_rhoN = found_K = found_psi = False
            if 'rhoN' in myline:
                found_rhoN = True
                found_density = found_K = found_psi = False
            if 'K' in myline:
                found_K = True
                found_density = found_rhoN = found_psi = False
            if 'psi' in myline:
                found_psi = True
                found_density = found_rhoN = found_K = False
            if found_density:
                densitylines.append(myline.rstrip('\n'))
            if found_rhoN:
                rhoNlines.append(myline.rstrip('\n'))
            if found_K:
                Klines.append(myline.rstrip('\n'))
            if found_psi:
                psilines.append(myline.rstrip('\n'))
    
    # strip out empty entries and convert string arrays to numerical ones,
    # stripping out last columns due to additional comma at ends of lines

    densitylines = [i for i in densitylines if i]
    density = np.genfromtxt(densitylines[1:], delimiter=',')
    density = np.delete(density, np.s_[-1:], axis=1)
    
    rhoNlines = [i for i in rhoNlines if i]
    rhoN = np.genfromtxt(rhoNlines[1:], delimiter=',')
    rhoN = np.delete(rhoN, np.s_[-1:], axis=1)

    Klines = [i for i in Klines if i]
    Kurv = np.genfromtxt(Klines[1:], delimiter=',')
    Kurv = np.delete(Kurv, np.s_[-1:], axis=1)
    
    psilines = [i for i in psilines if i]
    psi = np.genfromtxt(psilines[1:], delimiter=',')
    psi = np.delete(psi, np.s_[-1:], axis=1)

    # find extent of lattice and prepare (x,y) for plots
    row, col = density.shape
    x = np.linspace(0,col,col)
    y = np.linspace(0,row,row)
    X,Y = np.meshgrid(x,y)

    fig = plt.figure(figsize=[12.0, 16.0])
    fig.suptitle(title)
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.1, hspace=0.3)
    # produce plots of surface and contour for density
    ax1 = fig.add_subplot(4, 2, 1, projection='3d')
    ax1.plot_surface(X, Y, density, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax1.set_title('density surface')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax2 = fig.add_subplot(4, 2, 2)
    cf = ax2.contourf(X, Y, density, cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax2)
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_title('density map')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    # produce plots of surface and contour for phase index
    ax3 = fig.add_subplot(4, 2, 3, projection='3d')
    ax3.plot_surface(X, Y, rhoN, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax3.set_title('phase index surface')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax4 = fig.add_subplot(4, 2, 4)
    cf = ax4.contourf(X, Y, rhoN, cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax4)
    ax4.set_aspect('equal', adjustable='box')
    ax4.set_title('phase index map')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    # produce plots of surface and contour for curvature
    ax5 = fig.add_subplot(4, 2, 5, projection='3d')
    ax5.plot_surface(X, Y, Kurv, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax5.set_title('curvature surface')
    ax5.set_xlabel('x')
    ax5.set_ylabel('y')
    ax6 = fig.add_subplot(4, 2, 6)
    cf = ax6.contourf(X, Y, Kurv, cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax6)
    ax6.set_aspect('equal', adjustable='box')
    ax6.set_title('curvature map')
    ax6.set_xlabel('x')
    ax6.set_ylabel('y')
    # produce plots of surface and contour for velocity modulus (rescaled to actual values)
    ax7 = fig.add_subplot(4, 2, 7, projection='3d')
    ax7.plot_surface(X, Y, psi*0.001, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax7.set_title('velocity modulus surface')
    ax7.set_xlabel('x')
    ax7.set_ylabel('y')
    ax8 = fig.add_subplot(4, 2, 8)
    cf = ax8.contourf(X, Y, psi*0.001, cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax8)
    ax8.set_aspect('equal', adjustable='box')
    ax8.set_title('velocity modulus map')
    ax8.set_xlabel('x')
    ax8.set_ylabel('y')

def findCOMradius(filename, threshold=0.95):
    # read densities, phase indices and velocity moduli from CSV file
    found_density = False
    found_rhoN = False
    found_psi = False
    densitylines = []
    rhoNlines = []
    psilines = []
    with open(filename, 'rt') as myfile:
        for myline in myfile:
            if 'density' in myline:
                found_density = True
                found_rhoN = found_psi = False
            if 'rhoN' in myline:
                found_rhoN = True
                found_density = found_Psi = False
            if 'psi' in myline:
                found_psi = True
                found_density = found_rhoN = False
            if 'K' in myline:
                found_density = found_rhoN = found_psi = False
            if found_density:
                densitylines.append(myline.rstrip('\n'))
            if found_rhoN:
                rhoNlines.append(myline.rstrip('\n'))
            if found_psi:
                psilines.append(myline.rstrip('\n'))
    
    # strip out empty entries and convert string arrays to numerical ones,
    # stripping out last columns due to additional comma at ends of lines

    densitylines = [i for i in densitylines if i]
    density = np.genfromtxt(densitylines[1:], delimiter=',')
    density = np.delete(density, np.s_[-1:], axis=1)
    
    rhoNlines = [i for i in rhoNlines if i]
    rhoN = np.genfromtxt(rhoNlines[1:], delimiter=',')
    rhoN = np.delete(rhoN, np.s_[-1:], axis=1)

    psilines = [i for i in psilines if i]
    psi = np.genfromtxt(psilines[1:], delimiter=',')
    psi = np.delete(psi, np.s_[-1:], axis=1)
    
    # find extent of lattice
    row, col = density.shape

    # use phase indices to work out location of drop centre
    # (dealing with any crossing of periodic boundaries)
    zeta_x = xi_x = 0.0
    zeta_y = xi_y = 0.0
    mass = 0.0
    for j in range(row):
        omega_y = 2.0 * np.pi * float(j) / row
        for i in range(col):
            omega_x = 2.0 * np.pi * float(i) / col
            zeta_x += (1.0+rhoN[j,i]) * np.sin(omega_x)
            xi_x   += (1.0+rhoN[j,i]) * np.cos(omega_x)
            zeta_y += (1.0+rhoN[j,i]) * np.sin(omega_y)
            xi_y   += (1.0+rhoN[j,i]) * np.cos(omega_y)
            mass   += (1.0+rhoN[j,i])
    zeta_x /= mass
    xi_x   /= mass
    zeta_y /= mass
    xi_y   /= mass
    omega_x = np.arctan2(-zeta_x, -xi_x) + np.pi
    omega_y = np.arctan2(-zeta_y, -xi_y) + np.pi
    com_x = 0.5 * omega_x * col / np.pi
    com_y = 0.5 * omega_y * row / np.pi
    print("Centre of mass of drop = ({0:f}, {1:f})".format(com_x, com_y))

    # now find contour corresponding to phase index = 0 (boundary between phases)
    # and fit function for ellipse to find radii (radius)

    verts = ski.measure.find_contours(rhoN[:,:], 0.0)
    verts_periodic = []
    for contour in verts:
        for i in range(len(contour)):
            xx = contour[i, 1] - com_x
            yy = contour[i, 0] - com_y
            xx = xx - round(xx/col) * col
            yy = yy - round(yy/row) * row
            verts_periodic.append(np.asarray([xx,yy]))
    verts_periodic = np.asarray(verts_periodic)
    N = len(verts_periodic)
    x = verts_periodic[:, 0]
    y = verts_periodic[:, 1]
    U, S, V = np.linalg.svd(np.stack((x, y)))
    radii = np.sqrt(2/N)*S
    transform = np.sqrt(2/N) * U.dot(np.diag(S))
    print("Radii (semi-axes) of liquid drop: a = {0:f}, b = {1:f}".format(radii[0], radii[1]))
    radiusmean = (radii[0]*radii[1])**0.5
    print("Mean radius of liquid drop: {0:f}".format(radiusmean))

    # find contours for phase indices close to pure components (+/-threshold, default = 0.95)
    # and fit to functions for ellipse to find interfacial width

    if threshold < 0.0 or threshold > 1.0:
        threshold = 0.95
    verts1 = ski.measure.find_contours(rhoN[:,:], -threshold)
    verts2 = ski.measure.find_contours(rhoN[:,:],  threshold)
    verts1_periodic = []
    verts2_periodic = []
    for contour in verts1:
        for i in range(len(contour)):
            xx = contour[i, 1] - com_x
            yy = contour[i, 0] - com_y
            xx = xx - round(xx/col) * col
            yy = yy - round(yy/row) * row
            verts1_periodic.append(np.asarray([xx,yy]))
    verts1_periodic = np.asarray(verts1_periodic)
    N = len(verts1_periodic)
    x = verts1_periodic[:, 0]
    y = verts1_periodic[:, 1]
    U, S, V = np.linalg.svd(np.stack((x, y)))
    radii1 = np.sqrt(2/N)*S
    for contour in verts2:
        for i in range(len(contour)):
            xx = contour[i, 1] - com_x
            yy = contour[i, 0] - com_y
            xx = xx - round(xx/col) * col
            yy = yy - round(yy/row) * row
            verts2_periodic.append(np.asarray([xx,yy]))
    verts2_periodic = np.asarray(verts2_periodic)
    N = len(verts2_periodic)
    x = verts2_periodic[:, 0]
    y = verts2_periodic[:, 1]
    U, S, V = np.linalg.svd(np.stack((x, y)))
    radii2 = np.sqrt(2/N)*S
    width = ((radii1[0]*radii1[1])**0.5 - (radii2[0]*radii2[1])**0.5)
    print("Estimated interfacial width = {0:f}".format(width))
    
    # find densities inside and outside of drop, representing bulk values
    # of red and blue fluids, to allow us to work out interfacial tension

    y_drop = int(com_y+0.5)
    x_drop = int(com_x+0.5)
    bulk_blue = np.unravel_index(rhoN.argmin(), rhoN.shape)
    rho_red = density[y_drop,x_drop]
    rho_blue = density[bulk_blue[0],bulk_blue[1]]
    print("Density at drop centre = {0:f}, density in background fluid = {1:f}".format(rho_red, rho_blue))
    print("Estimated interfacial tension = {0:f}".format(radiusmean*(rho_red-rho_blue)/3.0))

    # find maximum velocity modulus (characteristic) and report

    umax = psi[:,:].max() * 0.001
    print("Maximum velocity modulus (characteristic microcurrent) = {0:e}".format(umax))

    # use nearest y-value to COM for density plot along x-axis (through drop)

    y_drop = int(com_y+0.5)
    xline = np.linspace(0,col,col)
    fig, ax = plt.subplots()
    ax.set(xlabel='x', ylabel='Density', title='Density plot through drop')
    ax.plot(xline, density[y_drop, :])
    
