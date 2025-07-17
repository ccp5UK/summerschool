from pathlib import Path
import subprocess
import shlex
import matplotlib.pyplot as plt
import numpy as np
import dlmhistoryread as dlmh
from tqdm.auto import tqdm

def calculateIsovolume(rundir='DPD2Ex1/Phase1', species='B', isosurfacesexe='dl_meso/WORK/isosurfaces.exe'):
    outdir = str(Path(rundir).resolve())
    command = str(Path(isosurfacesexe).resolve()) + " -b {0:s}".format(species)
    # quickly check HISTORY file for number of trajectory frames
    _, _, _, _, _, _, _, numframe, _ = dlmh.read_prepare(outdir+"/HISTORY")
    outfile = outdir+"/isoout"
    outpipe = open(outfile, 'w')
    # run isosurfaces.exe in directory and set up progress barout
    isorun = subprocess.Popen(shlex.split(command), cwd=outdir, stdout=outpipe)
    stepnum0 = 0
    pbar = tqdm(total=numframe, desc='Calculating isosurfaces in '+rundir)
    # keep checking isosurfaces.exe is running until it stops (according to our piped output file)
    if isorun.poll() is None:
        while stepnum0<numframe:
            with open(outfile, 'r') as file:
                lines = len(file.readlines())
                stepnum = lines//2
                if stepnum != stepnum0:
                    pbar.update(stepnum-stepnum0)
                    stepnum0 = stepnum
    # close down progress bar
    pbar.close()
    

def plotMoment(filename, name=''):
    data = np.loadtxt(filename)
    fig, ax = plt.subplots()
    plottitle = 'Plot of order parameters vs. time'
    if name!='':
        plottitle += ' for {0:s}'.format(name)
    ax.set(xlabel='Time', ylabel=r'Eigenvalues of second moments, $\mu_i$', title=plottitle)
    ax.plot(data[:, 0], data[:, 1], label = r'$\mu_1$')
    ax.plot(data[:, 0], data[:, 2], label = r'$\mu_2$')
    ax.plot(data[:, 0], data[:, 3], label = r'$\mu_3$')
    ax.legend(loc = "upper left")
    ax.set_ylim(bottom=0.0,top=1.0)
    ax.grid()

