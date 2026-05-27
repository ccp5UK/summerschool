import vtk
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import skimage as ski

def plotVTK(filename, title=''):
    # read all data from structured grid VTK file
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    # dimensions of box
    numX = data.GetDimensions()[0]
    numY = data.GetDimensions()[1]
    numZ = data.GetDimensions()[2] # should be 1
    # available data
    numdata = data.GetPointData().GetNumberOfArrays()
    density = []
    for i in range(numdata):
        namedata = data.GetPointData().GetArrayName(i)
        dataset = VN.vtk_to_numpy(data.GetPointData().GetArray(i))
        if namedata == 'velocity':
            dataset = np.hsplit(dataset,3)
            vx = dataset[0].reshape(numY, numX)
            vy = dataset[1].reshape(numY, numX)
        elif namedata.startswith('density'):
            density.append(dataset.reshape(numY, numX))
    # work out total densities, velocity modulus and phase indices
    numfluid = len(density)
    psi = np.sqrt(vx*vx + vy*vy)
    if numfluid>2:
        totalrho = density[0] + density[1] + density[2]
        rhodelim = density[0]+density[1]
        rhoN_1 = np.where(rhodelim/totalrho>=0.5, (density[1]-density[0])/rhodelim, -1.0)
        rhodelim = density[0]+density[2]
        rhoN_2 = np.where(rhodelim/totalrho>=0.5, (density[2]-density[0])/rhodelim, -1.0)
    else:
        totalrho = density[0] + density[1]
        rhoN_1 = (density[1] - density[0])/totalrho
        rhoN_2 = rhoN_1

    # prepare (x,y) for plots
    x = np.linspace(0,numX,numX)
    y = np.linspace(0,numY,numY)
    X,Y = np.meshgrid(x,y)

    #print(np.ptp(X), np.ptp(Y), np.ptp(totalrho))
    
    fig = plt.figure(figsize=[18.0, 8.0])
    fig.suptitle(title)
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.1, hspace=0.3)
    # produce contour plot for density
    ax1 = fig.add_subplot(2, 3, 1)
    cf = ax1.contourf(X, Y, totalrho[:,:], cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax1)
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_title('density map')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    # produce contour plot for phase index 1
    ax2 = fig.add_subplot(2, 3, 2)
    cf = ax2.contourf(X, Y, rhoN_1[:,:], cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax2)
    ax2.set_aspect('equal', adjustable='box')
    if numfluid>2:
        ax2.set_title('phase index map (drop 1)')
    else:
        ax2.set_title('phase index map')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    # produce contour plot for phase index 2
    if numfluid>2:
        ax3 = fig.add_subplot(2, 3, 3)
        cf = ax3.contourf(X, Y, rhoN_2[:,:], cmap=cm.coolwarm)
        fig.colorbar(cf, ax=ax3)
        ax3.set_aspect('equal', adjustable='box')
        ax3.set_title('phase index map (drop 2)')
        ax3.set_xlabel('x')
        ax3.set_ylabel('y')

    # produce contour plots for velocity modulus (speed) and x- and y-components 
    ax4 = fig.add_subplot(2, 3, 4)
    cs = ax4.contour(X, Y, rhoN_1[:,:], levels=[0.0], colors='black')
    if numfluid>2:
        cs1 = ax4.contour(X, Y, rhoN_2[:,:], levels=[0.0], colors='black')
    cf = ax4.contourf(X, Y, psi[:,:], cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax4)
    ax4.set_aspect('equal', adjustable='box')
    ax4.set_title('velocity modulus map')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax5 = fig.add_subplot(2, 3, 5)
    cs = ax5.contour(X, Y, rhoN_1[:,:], levels=[0.0], colors='black')
    if numfluid>2:
        cs1 = ax5.contour(X, Y, rhoN_2[:,:], levels=[0.0], colors='black')
    cf = ax5.contourf(X, Y, vx[:,:], cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax5)
    ax5.set_aspect('equal', adjustable='box')
    ax5.set_title('velocity x-component surface')
    ax5.set_xlabel('x')
    ax5.set_ylabel('y')
    ax6 = fig.add_subplot(2, 3, 6)
    cs = ax6.contour(X, Y, rhoN_1[:,:], levels=[0.0], colors='black')
    if numfluid>2:
        cs1 = ax6.contour(X, Y, rhoN_2[:,:], levels=[0.0], colors='black')
    cf = ax6.contourf(X, Y, vy[:,:], cmap=cm.coolwarm)
    fig.colorbar(cf, ax=ax6)
    ax6.set_aspect('equal', adjustable='box')
    ax6.set_title('velocity y-component map')
    ax6.set_xlabel('x')
    ax6.set_ylabel('y')

def findCOMradius(filename, threshold=0.95):
    # read all data from structured grid VTK file
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    # dimensions of box
    numX = data.GetDimensions()[0]
    numY = data.GetDimensions()[1]
    numZ = data.GetDimensions()[2] # should be 1
    # available data
    numdata = data.GetPointData().GetNumberOfArrays()
    density = []
    for i in range(numdata):
        namedata = data.GetPointData().GetArrayName(i)
        dataset = VN.vtk_to_numpy(data.GetPointData().GetArray(i))
        if namedata == 'velocity':
            dataset = np.hsplit(dataset,3)
            vx = dataset[0].reshape(numY, numX)
            vy = dataset[1].reshape(numY, numX)
        elif namedata.startswith('density'):
            density.append(dataset.reshape(numY, numX))
    # work out total densities, velocity modulus and phase indices
    numfluid = len(density)
    psi = np.sqrt(vx*vx + vy*vy)
    if numfluid>2:
        totalrho = density[0] + density[1] + density[2]
        rhoN_1 = (density[1] - density[0])/(density[0]+density[1])
        rhoN_2 = (density[2] - density[0])/(density[0]+density[2])
    else:
        totalrho = density[0] + density[1]
        rhoN_1 = (density[1] - density[0])/totalrho
        rhoN_2 = rhoN_1

    # use phase indices to work out location of centre for drop 1
    # (dealing with any crossing of periodic boundaries)
    zeta_x = xi_x = 0.0
    zeta_y = xi_y = 0.0
    mass = 0.0
    for j in range(numY):
        omega_y = 2.0 * np.pi * float(j) / numY
        for i in range(numX):
            omega_x = 2.0 * np.pi * float(i) / numX
            zeta_x += (1.0+rhoN_1[j,i]) * np.sin(omega_x)
            xi_x   += (1.0+rhoN_1[j,i]) * np.cos(omega_x)
            zeta_y += (1.0+rhoN_1[j,i]) * np.sin(omega_y)
            xi_y   += (1.0+rhoN_1[j,i]) * np.cos(omega_y)
            mass   += (1.0+rhoN_1[j,i])
    zeta_x /= mass
    xi_x   /= mass
    zeta_y /= mass
    xi_y   /= mass
    omega_x = np.arctan2(-zeta_x, -xi_x) + np.pi
    omega_y = np.arctan2(-zeta_y, -xi_y) + np.pi
    com1_x = 0.5 * omega_x * numX / np.pi
    com1_y = 0.5 * omega_y * numY / np.pi
    print("Centre of mass of drop 1 = ({0:f}, {1:f})".format(com1_x, com1_y))

    # check if fluid 1 actually is at centre of drop: if not, assume drop has broken up

    broken1 = False
    y_drop = int(com1_y+0.5)
    x_drop = int(com1_x+0.5)
    if rhoN_1[y_drop, x_drop]<0.0:
        print("Cannot find fluid 1 at centre of mass for drop 1: assumed to have broken up")
        broken1 = True

    # use phase indices to work out location of centre for drop 2 (if available)
    # (dealing with any crossing of periodic boundaries)
    broken2 = False
    if numfluid>2:
        zeta_x = xi_x = 0.0
        zeta_y = xi_y = 0.0
        mass = 0.0
        for j in range(numY):
            omega_y = 2.0 * np.pi * float(j) / numY
            for i in range(numX):
                omega_x = 2.0 * np.pi * float(i) / numX
                zeta_x += (1.0+rhoN_2[j,i]) * np.sin(omega_x)
                xi_x   += (1.0+rhoN_2[j,i]) * np.cos(omega_x)
                zeta_y += (1.0+rhoN_2[j,i]) * np.sin(omega_y)
                xi_y   += (1.0+rhoN_2[j,i]) * np.cos(omega_y)
                mass   += (1.0+rhoN_2[j,i])
        zeta_x /= mass
        xi_x   /= mass
        zeta_y /= mass
        xi_y   /= mass
        omega_x = np.arctan2(-zeta_x, -xi_x) + np.pi
        omega_y = np.arctan2(-zeta_y, -xi_y) + np.pi
        com2_x = 0.5 * omega_x * numX / np.pi
        com2_y = 0.5 * omega_y * numY / np.pi
        print("Centre of mass of drop 2 = ({0:f}, {1:f})".format(com2_x, com2_y))
    else:
        com2_x = com2_y = 0.0

        # check if fluid 1 actually is at centre of drop: if not, assume drop has broken up

        y_drop = int(com1_y+0.5)
        x_drop = int(com1_x+0.5)
        if rhoN_2[y_drop, x_drop]<0.0:
            print("Cannot find fluid 2 at centre of mass for drop 2: assumed to have broken up")
            broken2 = True


    # now find contour corresponding to phase index = 0 (boundary between phases)
    # and fit function for ellipse to find radii (radius) of drop 1

    if not broken1:
        verts = ski.measure.find_contours(rhoN_1[:,:], 0.0)
        verts_periodic = []
        for contour in verts:
            for i in range(len(contour)):
                xx = contour[i, 1] - com1_x
                yy = contour[i, 0] - com1_y
                xx = xx - round(xx/numX) * numX
                yy = yy - round(yy/numY) * numY
                verts_periodic.append(np.asarray([xx,yy]))
        verts_periodic = np.asarray(verts_periodic)
        N = len(verts_periodic)
        x = verts_periodic[:, 0]
        y = verts_periodic[:, 1]
        U, S, V = np.linalg.svd(np.stack((x, y)))
        radii = np.sqrt(2/N)*S
        transform = np.sqrt(2/N) * U.dot(np.diag(S))
        print("Radii (semi-axes) of liquid drop 1: a = {0:f}, b = {1:f}".format(radii[0], radii[1]))
        print("Mean radius of liquid drop 1: {0:f}".format((radii[0]*radii[1])**0.5))
        radius1 = (radii[0]*radii[1])**0.5

    # now find contour corresponding to phase index = 0 (boundary between phases)
    # and fit function for ellipse to find radii (radius) of drop 2 (if available)

    if numfluid>2 and not broken2:
        verts = ski.measure.find_contours(rhoN_2[:,:], 0.0)
        verts_periodic = []
        for contour in verts:
            for i in range(len(contour)):
                xx = contour[i, 1] - com2_x
                yy = contour[i, 0] - com2_y
                xx = xx - round(xx/numX) * numX
                yy = yy - round(yy/numY) * numY
                verts_periodic.append(np.asarray([xx,yy]))
        verts_periodic = np.asarray(verts_periodic)
        N = len(verts_periodic)
        x = verts_periodic[:, 0]
        y = verts_periodic[:, 1]
        U, S, V = np.linalg.svd(np.stack((x, y)))
        radii = np.sqrt(2/N)*S
        transform = np.sqrt(2/N) * U.dot(np.diag(S))
        print("Radii (semi-axes) of liquid drop 2: a = {0:f}, b = {1:f}".format(radii[0], radii[1]))
        print("Mean radius of liquid drop 2: {0:f}".format((radii[0]*radii[1])**0.5))
        radius2 = (radii[0]*radii[1])**0.5

    # find contours for phase indices close to pure components (+/-threshold, default = 0.95)
    # and fit to functions for ellipse to find interfacial width for drop 1 (and drop 2 if available)

    if threshold < 0.0 or threshold > 1.0:
        threshold = 0.95
    if not broken1:
        verts1 = ski.measure.find_contours(rhoN_1[:,:], -threshold)
        verts2 = ski.measure.find_contours(rhoN_1[:,:],  threshold)
        verts1_periodic = []
        verts2_periodic = []
        for contour in verts1:
            for i in range(len(contour)):
                xx = contour[i, 1] - com1_x
                yy = contour[i, 0] - com1_y
                xx = xx - round(xx/numX) * numX
                yy = yy - round(yy/numY) * numY
                verts1_periodic.append(np.asarray([xx,yy]))
        verts1_periodic = np.asarray(verts1_periodic)
        N = len(verts1_periodic)
        x = verts1_periodic[:, 0]
        y = verts1_periodic[:, 1]
        U, S, V = np.linalg.svd(np.stack((x, y)))
        radii1 = np.sqrt(2/N)*S
        for contour in verts2:
            for i in range(len(contour)):
                xx = contour[i, 1] - com1_x
                yy = contour[i, 0] - com1_y
                xx = xx - round(xx/numX) * numX
                yy = yy - round(yy/numY) * numY
                verts2_periodic.append(np.asarray([xx,yy]))
        verts2_periodic = np.asarray(verts2_periodic)
        N = len(verts2_periodic)
        x = verts2_periodic[:, 0]
        y = verts2_periodic[:, 1]
        U, S, V = np.linalg.svd(np.stack((x, y)))
        radii2 = np.sqrt(2/N)*S
        width = ((radii1[0]*radii1[1])**0.5 - (radii2[0]*radii2[1])**0.5)
        print("Estimated interfacial width for drop 1 = {0:f}".format(width))
    if numfluid>2 and not broken2:
        verts1 = ski.measure.find_contours(rhoN_2[:,:], -threshold)
        verts2 = ski.measure.find_contours(rhoN_2[:,:],  threshold)
        verts1_periodic = []
        verts2_periodic = []
        for contour in verts1:
            for i in range(len(contour)):
                xx = contour[i, 1] - com2_x
                yy = contour[i, 0] - com2_y
                xx = xx - round(xx/numX) * numX
                yy = yy - round(yy/numY) * numY
                verts1_periodic.append(np.asarray([xx,yy]))
        verts1_periodic = np.asarray(verts1_periodic)
        N = len(verts1_periodic)
        x = verts1_periodic[:, 0]
        y = verts1_periodic[:, 1]
        U, S, V = np.linalg.svd(np.stack((x, y)))
        radii1 = np.sqrt(2/N)*S
        for contour in verts2:
            for i in range(len(contour)):
                xx = contour[i, 1] - com2_x
                yy = contour[i, 0] - com2_y
                xx = xx - round(xx/numX) * numX
                yy = yy - round(yy/numY) * numY
                verts2_periodic.append(np.asarray([xx,yy]))
        verts2_periodic = np.asarray(verts2_periodic)
        N = len(verts2_periodic)
        x = verts2_periodic[:, 0]
        y = verts2_periodic[:, 1]
        U, S, V = np.linalg.svd(np.stack((x, y)))
        radii2 = np.sqrt(2/N)*S
        width = ((radii1[0]*radii1[1])**0.5 - (radii2[0]*radii2[1])**0.5)
        print("Estimated interfacial width for drop 2 = {0:f}".format(width))
        
    # only for two-fluid systems: find densities representing bulk values
    # of red and blue fluids, to allow us to work out interfacial tension

    if not broken1 and numfluid==2:
        y_drop = int(com1_y+0.5)
        x_drop = int(com1_x+0.5)
        bulk_blue = np.unravel_index(rhoN_1.argmin(), rhoN_1.shape)
        rho_red = totalrho[y_drop,x_drop]
        rho_blue = totalrho[bulk_blue[0],bulk_blue[1]]
        print("Density at centre of drop = {0:f}, density in background fluid = {1:f}".format(rho_red, rho_blue))
        print("Estimated interfacial tension = {0:f}".format(radius1*(rho_red-rho_blue)/3.0))
        
    # find maximum velocity modulus (characteristic) and report

    umax = psi[:,:,].max()
    print("Maximum velocity modulus = {0:e}".format(umax))
    