import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def plotConcentration(filename, dataset=1):
    numset,lx,concdata,datumNames,chi,chi_legend = readFloryHuggins(filename)
    if dataset<1 or dataset>numset:
        print("Selected dataset out of range (1 - {0:d}) for file {1:s}".format(numset, filename))
        return
    set = dataset-1
    fig, ax = plt.subplots()
    ax.set(xlabel='x [DPD length units]', ylabel=r'Bead concentration, $\phi$', title=datumNames[set])
    ax.set_xlim(left=min(concdata[2*set]), right=max(concdata[2*set]))
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.plot(concdata[2*set], concdata[2*set+1], 'r-')
    ax.legend([r'$\chi$ = {0:f}$\pm${1:f}'.format(chi[set][2], chi[set][3])], loc='upper right')
    ax.grid()

def plotChi(filename, dAmin=10.0, dAmax=100.0):
    numset,lx,concdata,datumNames,chi,chi_legend = readFloryHuggins(filename)
    nAii = len(chi_legend)
    fullset = []
    colourmap = mpl.colormaps['Dark2']
    fig, ax = plt.subplots()
    for i in range(nAii):
        fullset.append([x for x in chi if x[0]==chi_legend[i]])
    colours = colourmap(np.linspace(0, 1, nAii))
    for i in range(nAii):
        tmpx = [(x[1]-x[0]) for x in fullset[i]]
        tmpy = [x[2] for x in fullset[i]]
        tmpe = [x[3] for x in fullset[i]]
        ax.scatter(tmpx, tmpy, color=colours[i], marker='s', label='Aii = {0:s}'.format(str(chi_legend[i])), zorder=1)
        ax.errorbar(tmpx, tmpy, yerr=tmpe, linestyle='None', fmt='None', ecolor='gray', capsize=5, zorder=1)
    xchi = [(x[1]-x[0]) for x in chi]
    chi = [x[2] for x in chi]
    ax.set_xlabel(r'$\Delta$A')
    ax.set_ylabel(r'$\chi$')
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    plotx=[0.0]
    xchidata=[]
    chidata=[]
    for i in range(len(chi)):
        x = xchi[i]
        plotx.append(x)
        if(xchi[i]>=dAmin and xchi[i]<=dAmax):
            xchidata.append(xchi[i])
            chidata.append(chi[i]/x)
    plotx = sorted(plotx)
    a = sum(chidata)/len(chidata)
    ra = 1.0/a
    amax = max(chidata)
    ramin = 1.0/amax
    amin = min(chidata)
    ramax = 1.0/amin
    aerr = max(amax-a, a-amin)
    raerr = max(ramax-ra, ra-ramin)
    fit_fn = np.poly1d([a, 0.0])
    ax.plot(plotx, fit_fn(plotx), '--k')
    ax.set_title(r'$\chi$ = ({0:f}$\pm${1:f})$\Delta$A'.format(a, aerr)+'\n'+r'$\Delta$A = ({0:f}$\pm${1:f})$\chi$'.format(ra, raerr))
    ax.legend()

def readFloryHuggins(filename):
    s = open(filename).read().split('\n')
    i = 0
    chidata = []
    datumNames = []
    Aiiset = []
    d = []
    l = []
    while i<len(s):
        line = s[i].split(',')
        if len(line)!=5:
            break
        Aii = float(line[0])
        Aij = float(line[1])
        chi = float(line[2])
        chisd = float(line[3])
        datumNames.append("Aii = {0:s}, Aij = {1:s}".format(str(Aii), str(Aij)))
        chidata.append([Aii, Aij, chi, chisd])
        datalength = int(line[4])
        Aiiset.append(Aii)
        pointdata = []
        ss = s[i+1:i+datalength+1]
        for line in ss:
            dataelement = [float(j) for j in line.split()]
            pointdata.append(dataelement)
        d.append([x[0] for x in pointdata])
        d.append([x[1] for x in pointdata])
        i += (datalength+1)
    legend = list(sorted(set(Aiiset)))
    n = len(datumNames)
    for i in range(n):
        dx = d[2*i][1]-d[2*i][0]
        maxlx = max(d[2*i][:])
        l.append(maxlx+0.5*dx)
    return n,l,d,datumNames,chidata,legend