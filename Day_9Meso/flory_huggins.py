#!/usr/bin/env python3
"""Usage:
    flory_huggins.py [--rho <rho>] [--Aii <Aii>] [--Aijmin <Aijmin>] 
                     [--Aijmax <Aijmax>] [--dA <dA>] [--dz <dz>] 
                     [--L <L>] [--W <W>] [--dlmeso <dlmeso>] [--out <out>] 
                     [--nproc <nproc>]

Carries out DL_MESO_DPD calculations to determine relationship between
Flory-Huggins chi parameters and conservative force parameters, analyse
data and plot results

Options:
    -h --help           Print this message
    --rho <rho>         Particle density [default: 3.0]
    --Aii <Aii>         Conservative force parameter for like-like particle
                        interactions [default: 25.0]
    --Aijmin <Aijmin>   Mininum conservative force parameter between particle
                        species [default: 33.0]
    --Aijmax <Aijmax>   Maximum conservative force parameter between particle
                        species [default: 43.0]
    --dA <dA>           Steps between values of conservative force parameter
                        between particle species for each run [default: 1.0]
    --dz <dz>           Bin size for concentration profile in z-direction
                        [default: 0.05]
    --L <L>             Length of box in z-direction [default: 10.0]
    --W <W>             Width of box in x- and y-directions [default: 6.0]
    --dlmeso <dlmeso>   Location of DL_MESO_DPD executable [default: ./dpd.exe]
    --out <out>         Folder for running DL_MESO_DPD calculations [default: out]
    --nproc <nproc>     Number of processor cores to run DL_MESO_DPD calculations [default: 1]

michael.seaton@stfc.ac.uk, 15/07/26
"""
from docopt import docopt
from pathlib import Path
from tqdm.auto import tqdm
import dlmoutputread as dlmo
import subprocess
import shlex
import sys
import numpy as np
import statistics
import math
import os

args = docopt(__doc__)
rho = float(args["--rho"])
Aii = float(args["--Aii"])
Aijmin = float(args["--Aijmin"])
Aijmax = float(args["--Aijmax"])
dA = float(args["--dA"])
dz = float(args["--dz"])
L = float(args["--L"])
W = float(args["--W"])
dlmeso = str(Path(args["--dlmeso"]).resolve())
out = str(Path(args["--out"]).resolve())
nproc = int(args["--nproc"])

bo = sys.byteorder
if(bo == 'big'):
    ri = ">i"
    rd = ">d"
else:
    ri = "<i"
    rd = "<d"

if nproc > 1 :
    invoke = "mpirun -np {0:d} {1:s}".format(nproc, dlmeso)
else :
    invoke = dlmeso

intsize = 4
longintsize = 8

gam = 4.5
Npart = int(round(0.5*rho*L*W*W))
if Npart%2 == 1:
    Npart += 1

factors = []
for i in range(2,Npart):
    if Npart%i == 0:
        factors.append(i)

boxlist = []
errors = []
for i in range(len(factors)):
    Nxtry = Npart//(factors[i]*factors[i])
    Nparttry = Nxtry * factors[i] * factors[i]
    if Npart==Nparttry:
        error = (2.*Nxtry/L - factors[i]/W) ** 2
        listitem = [Nxtry, factors[i]]
        boxlist.append(listitem)
        errors.append(error)

minerror = errors.index(min(errors))
Nz = boxlist[minerror][0]
Nx = boxlist[minerror][1]

print("Box size: {0:.4f} by {1:.4f} by {2:.4f}, density: {3:.4f}".format(W, W, L, rho))
print("Number of particles: {0:d} ({1:d} by {2:d} by {3:d})".format(2*Npart, Nx, Nx, 2*Nz))
print("Like-like interaction parameter (Aii): {0:f}".format(Aii))
print("Interaction parameters between species (Aij): {0:f} to {1:f}".format(Aijmin, Aijmax))
print("Writing data to file: {0:s}/floryhuggins-rho-{1:.3f}.dat".format(out, rho))

# create folder for running simulations (if does not already exist)

os.makedirs(out, exist_ok=True)

# create CONTROL file - same for all simulations

sc = "DL_MESO Flory-Huggins chi-parameter determination\n\n"
sc += "volume {0:.4f} {1:.4f} {2:.4f}\n".format(L, W, W)
sc += "temperature 1.0\n"
sc += "cutoff 1.0\n\n"
sc += "timestep 0.01\n"
sc += "steps 70000\n"
sc += "equilibration steps 20000\n"
sc += "zden sampling every 100\n"
sc += "zden binsize {0:f}\n".format(dz)
sc += "stack size 100\n"
sc += "print every 100\n"
sc += "job time 3600.0\n"
sc += "close time 200.0\n\n"
sc += "ensemble nvt mdvv\n\n"
sc += "l_scr\n\n"
sc += "finish\n"

open(out+"/CONTROL", "w").write(sc)
print("CONTROL file saved.")

# create CONFIG file

wz = 0.5*L/Nz
wx = W/Nx
hz = 0.5*L
hx = 0.5*W
cf = "DL_MESO Flory-Huggins chi-parameter determination\n"
cf += "0\t2\t{0:d}\n".format(2*Npart)
cf += "{0:16.10f}{1:16.10f}{2:16.10f}\n".format(W, 0.0, 0.0)
cf += "{0:16.10f}{1:16.10f}{2:16.10f}\n".format(0.0, W, 0.0)
cf += "{0:16.10f}{1:16.10f}{2:16.10f}\n".format(0.0, 0.0, L)
for i in range(Npart):
    cf += "{0:s}        {1:d}\n".format("A", i+1)
    ix = i%Nx
    iz = i//(Nx*Nx)
    iy = (i%(Nx*Nx))//Nx
    xx = (ix+0.5)*wx-hx
    yy = (iy+0.5)*wx-hx
    zz = (iz+0.5)*wz-hz
    cf += "{0:16.10f}{1:16.10f}{2:16.10f}\n".format(xx, yy, zz)
for i in range(Npart):
    cf += "{0:s}        {1:d}\n".format("B", Npart+i+1)
    ix = i%Nx
    iz = i//(Nx*Nx)
    iy = (i%(Nx*Nx))//Nx
    xx = (ix+0.5)*wx-hx
    yy = (iy+0.5)*wx-hx
    zz = (iz+0.5)*wz
    cf += "{0:16.10f}{1:16.10f}{2:16.10f}\n".format(xx, yy, zz)

open(out+"/CONFIG", "w").write(cf)
print("CONFIG file saved.")

# open file for recording concentration profiles and chi values

filename = out+'/floryhuggins-rho-{0:.3f}.dat'.format(rho)
fw = open(filename, "a+")

# loop through unlike interaction parameters

for Aij in np.arange(Aijmin, Aijmax+0.5*dA, dA):

# create FIELD file

    sf = "DL_MESO Flory-Huggins chi-parameter determination\n\n"
    sf += "SPECIES 2\n"
    sf += "A 1.0 0.0 {0:d} 0\n".format(Npart)
    sf += "B 1.0 0.0 {0:d} 0\n\n".format(Npart)
    sf +="INTERACTIONS 3\n"
    sf += "A A dpd  {0:.3f} 1.0\n".format(Aii)
    sf += "A B dpd  {0:.3f} 1.0\n".format(Aij)
    sf += "B B dpd  {0:.3f} 1.0\n\n".format(Aii)
    sf +="THERMOSTAT 3\n"
    sf += "A A quad {0:.3f} 1.0\n".format(gam)
    sf += "A B quad {0:.3f} 1.0\n".format(gam)
    sf += "B B quad {0:.3f} 1.0\n\n".format(gam)
    sf += "CLOSE\n"

    open(out+"/FIELD", "w").write(sf)
    print("FIELD file saved for Aij = {0:f}.".format(Aij))

# run DL_MESO_DPD and monitor simulation progress

    description = "Aii = {0:f}, Aij = {1:f}".format(Aii, Aij)
    outfile = out+"/OUTPUT-Aii-{0:f}-Aij-{1:f}-rho-{2:f}".format(Aii, Aij, rho)
    outpipe = open(outfile, 'w')
    dlmesorun = subprocess.Popen(shlex.split(invoke), cwd=out, stdout=outpipe)
    stepnum0 = 0
    pbar = tqdm(total=70000, desc='Running '+description)
    terminate = False
    if dlmesorun.poll() is None:
        while not terminate:
            if os.path.getsize(outfile)>0:
                _, startrun, _, _, numstep, terminate, _ = dlmo.read_prepare(outfile)
                if startrun>0:
                    rundata, _, _, _ = dlmo.read_run(outfile, startrun, terminate)
                    stepnum = int(rundata[-1,0])
                    if stepnum != stepnum0:
                        pbar.update(stepnum-stepnum0)
                        stepnum0 = stepnum
    pbar.close()
    
    # open ZDNDAT file and get hold of density profiles
    # for species A and all species to calculate concentration profile

    fz = open(out+"/ZDNDAT", "r")
    content = fz.readlines()
    numlines = len(content)
    # second line contains number of data sets (should be 3) and
    # number of lines per set
    words = content[1].split()
    if (int(words[0])!=3):
        sys.exit("ERROR: ZDNDAT file does not contain data for two species - cannot continue!")
    nz = int(words[1])
    # setup arrays for sampled density profiles
    rhoall = np.zeros(nz)
    rhospec = np.zeros(nz)
    # now look for first available z-density profile for species named "A"
    for line in range(2, numlines):
        if "A" in content[line]:
            for data in range(nz):
                words = content[line+data+1].split()
                rhoz = float(words[1])
                rhospec[data] = rhoz
            break
    # now look for data set starting with "all species":
    # this will be the density profile for all bead species
    for line in range(2, len(content)):
        if "all species" in content[line]:
            for data in range(nz):
                words = content[line+data+1].split()
                rhoz = float(words[1])
                rhoall[data] = rhoz
            break
    # finally convert densities to concentrations (volume fractions)
    # and rename ZDNDAT file for run
    volfrac = rhospec/rhoall
    os.rename(out+'/ZDNDAT', out+'/ZDNDAT-Aii-{0:f}-Aij-{1:f}-rho-{2:f}'.format(Aii, Aij, rho))

    # select region of simulation box where species are likely
    # to be separated out (but hopefully not entirely!) and use
    # concentrations in region to estimate chi value and its error:
    # to better avoid completely separated regions, we will base
    # our search for a region with low concentration of A rather 
    # than high (numerically more likely to be above 0 than below 1)
    
    minchi = int(0.65*nz)
    maxchi = int(0.85*nz)
    meanvolfrac = statistics.mean(volfrac[minchi:maxchi])
    stdvolfrac = statistics.stdev(volfrac[minchi:maxchi])
    if meanvolfrac>0.0:
        chi = math.log((1.0-meanvolfrac)/meanvolfrac)/(1.0-2.0*meanvolfrac)
        minvolfrac = meanvolfrac - stdvolfrac
        maxvolfrac = meanvolfrac + stdvolfrac
        if minvolfrac>0.0:
            chimax = math.log((1.0-minvolfrac)/minvolfrac)/(1.0-2.0*minvolfrac)
        else:
            chimax = chi
        chimin = math.log((1.0-maxvolfrac)/maxvolfrac)/(1.0-2.0*maxvolfrac)
        chierr = max(abs(chimax-chi), abs(chi-chimin))
        # write force parameters, chi, error and number of histogram bins to file,
        # before writing concentration profile (can plot/check this later)
        fw.write("{0:f},{1:f},{2:f},{3:f},{4:d}\n".format(Aii, Aij, chi, chierr, nz))
        for i in range(nz):
            fw.write("{0:11.7f}      {1:11.7f}\n".format((i+0.5)*dz, volfrac[i]))
        fw.flush()
        print("Written data for Aii = {0:f}, Aij = {1:f} to {2:s} - chi = {3:f} +/- {4:f}".format(Aii, Aij, filename, chi, chierr))
    else:
        print("Beads completely separated for Aij = {0:f} - cannot reliably calculate a value of chi!".format(Aij))


