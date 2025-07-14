#!/usr/bin/env python3
# hdf5_module.py

# DISCLAIMER
# (c) 2022 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
# The author makes no warranties about, and disclaims liability for all uses of, this software.
# The author does not recommend use of this software for any purpose.

"""Simple interface to read a flat HDF5 file."""

def read_file ( file ):
    """Read attributes and datasets from file.

    Argument
    --------
    file : string
       HDF5 file containing attributes and datasets but no groups

    Returns
    -------
    dict
       Dictionary of attributes
    dict
       Dictionary of NumPy arrays containing datasets
    """

    import h5py
    print('Opening file ',file)
    with h5py.File(file,'r') as f:
        print('Attributes ',list(f.attrs.keys()))
        attributes={}
        for key in f.attrs:
            attributes[key]=f.attrs[key]
        print('Datasets   ',list(f.keys()))
        datasets={}
        for key in f:
            datasets[key]=f[key][...]
        print('Done')
        return attributes, datasets



