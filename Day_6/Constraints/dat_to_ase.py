#!/usr/bin/env python3

"""Convert LJ-units configuration in pbc into ASE format for Ar atoms"""

import numpy as np
import sys
from ase.atoms import Atoms

def argon(file):

    sigma = 3.4  # Diameter in Angstroms for Ar

    with open(file) as f:
        first_line = f.readline()
        n = int(first_line)
        print('Number of atoms N = ',n)

    # Read in triplets, skipping first row and discarding any velocities
    r = np.genfromtxt(file,skip_header=1,usecols=(0,1,2))
    box=r[0,:] if r.ndim==2 else r[:] # First triplet is actually box
    print('Box dimensions  = ',box)

    # We need to cope with N=0 atoms
    if n<=0:
        box = box*sigma
        return [Atoms(symbols=[],positions=[],pbc=True,cell=box)]

    r = r[1:,:] # Strip off first triplet
    assert r.shape[0]==n, 'N mismatch'

    # Apply periodic boundary conditions (box is broadcast across first index)
    # Shift coordinates to lie in range [0,box]
    r = r / box
    r = r - np.rint ( r )
    r = r + 0.5
    r = r * box

    box   = box*sigma
    r     = r * sigma

    ar=np.full(n,'Ar')
    return [Atoms(symbols=ar,positions=r,pbc=True,cell=box)]
