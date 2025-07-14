#!/usr/bin/env python3

"""Convert LJ-units configuration into XYZ format for Ar atoms"""

import numpy as np
import sys

file = input('Enter configuration file name: ')
print('Reading from file name : ',file)

with open(file) as f:
    first_line = f.readline()
    n = int(first_line)
    print('Number of atoms = ',n)

if n <= 0:
    print('No coordinates')
    sys.exit()

# Read in triplets, skipping first row and discarding any velocities
r = np.genfromtxt(file,skip_header=1,usecols=(0,1,2))
box=r[0,:]            # First triplet is actually box
print('Box dimensions  = ',box)
r = r[1:,:]           # Strip off first triplet
print('R array shape   = ',r.shape)

# Apply periodic boundary conditions (may not be necessary)
r = r / box           # Box should be broadcast across first index
r = r - np.rint ( r ) # Periodic boundary correction
r = r * box           # Box should be broadcast across first index

r = r * 3.4           # Multiply by sigma in Angstroms for Ar
header = first_line + "Comment line"
file=file[:-4] + ".xyz"
print('Writing to file name: ',file)
np.savetxt(file,r,header=header,comments='',fmt="Ar   %15.8f %15.8f %15.8f")
