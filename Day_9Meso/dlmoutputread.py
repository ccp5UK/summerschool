#!/usr/bin/env python3
"""DL_MESO_DPD OUTPUT file readers

Module to read DL_MESO_DPD OUTPUT files:

michael.seaton@stfc.ac.uk, 28/06/23
"""

import sys
import math
import numpy as np
import os

def read_prepare(filename):
    """Scans DL_MESO_DPD OUTPUT file to find essential information for reading further"""
    
    # inputs:
    #   filename        name of OUTPUT file to start reading
    # outputs:
    #   numlines        total number of lines in OUTPUT file
    #   startrun        line number where first timestep in run is recorded (return -1 if does not exist)
    #   startstep       number of first timestep in run (return -1 if does not exist)
    #   termstep        number of timestep at which calculation has been terminated (returns 0 if does not exist)
    #   numstep         number of available timesteps in run (return 0 if cannot find any)
    #   terminate       flag indicating if calculation has terminated properly
    #   datanames       names of columns for data for each timestep

    # open OUTPUT file and split into lines

    try:
        with open(filename) as file:
            content = file.read().splitlines()

        numlines = len(content)
        startrun = -1
        startstep = -1
        numstep = 0
        termstep = 0
        datanames = []
        terminate = False

    # scan through lines to find first line with run results
    # (indicated as a series of dashes)

        if numlines>0:
            for line in range(numlines):
                if "---" in content[line]:
                    startrun = line
                    break

    # look at next line to find column names - all but first one
    # are names of properties being printed - and data line afterwards
    # to find number of first timestep in run

        if startrun>-1 and startrun+3<numlines:
            names = content[startrun+1].split()
            datanames = names[1:]
            words = content[startrun+3].split()
            if len(words)>1:
                startstep = int(words[0])

    # now scan through available data lines to see how many timesteps
    # are available in file - checks for lines with dashes and sees 
    # if third line after also contains dashes, which should skip
    # over any blank lines and column headers - and stop if 
    # "run terminating" is indicated (to avoid reading in averaged
    # properties afterwards)

        if startstep>-1:
            for line in range(startrun, numlines):
                if "---" in content[line] and line+3<numlines:
                    if "---" in content[line+3]:
                        numstep += 1
                elif "run terminating" in content[line]:
                    terminate = True
                    if line+3<=numlines:
                        words = content[line+3].split()
                        termstep = int(words[4])
                    break

    except FileNotFoundError:
        print("ERROR: Cannot open OUTPUT file")
    
    return numlines, startrun, startstep, termstep, numstep, terminate, datanames

def read_run(filename,startrun,terminate):
    """Reads statistical properties written to OUTPUT file during DL_MESO_DPD calculation, including any averaged values"""

    # inputs:
    #   filename        name of OUTPUT file to read
    #   startrun        line number where simulation run starts in OUTPUT file
    #   terminate       flag indicating if simulation terminated properly (and averaged values are available)
    # outputs:
    #   rundata         statistical properties read from OUTPUT file during DL_MESO_DPD calculation:
    #                   each entry includes timestep number, calculation walltime, instantaneous
    #                   values for each property, rolling average values for each property
    #   averages        averaged values for each property, including pressure tensors (conservative, dissipative, 
    #                   random and kinetic contributions, plus overall values)
    #   fluctuations    fluctuations (standard deviations) for each property, including pressure tensors
    #                   (conservative, dissipative, random and kinetic contributions, plus overall values)
    #   datanames       names of data for averaged values and fluctuations (including pressure tensors)

    with open(filename) as file:
        content = file.read().splitlines()

    numlines = len(content)
    rundata = []
    avelines = 0

    # go through all lines in OUTPUT file where data is available,
    # and find all available timesteps - indicated by two lines of
    # dashes separated by two lines of numbers - then read in data
    # and add to list (also work out where data ends and any averages
    # can be found)
 
    for line in range(startrun,numlines):
        if "---" in content[line] and line+3<numlines:
            if "---" in content[line+3]:
                words = content[line+1].split()
                timestep = int(words[0])
                instantdata = list(map(float, words[1:]))
                words = content[line+2].split()
                walltime = float(words[0])
                runningdata = list(map(float, words[1:]))
                data = [timestep, walltime]
                data += instantdata
                data += runningdata
                rundata.append(data)
        elif "run terminating" in content[line]:
            avelines = line
            break

    rundata = np.array(rundata)
    
    # if available, look for final averages and fluctuations for
    # properties and pressure tensors and read these into lists

    averages = []
    fluctuations = []
    datanames = []
    totaltensor = False
    numtensor = 0

    if terminate:
        for line in range(avelines,numlines):
            if "---" in content[line] and line+4<numlines:
                if "---" in content[line+2]:
                    names = content[startrun+1].split()
                    datanames = names[1:]
                    words = content[line+3].split()
                    data = list(map(float, words))
                    averages.extend(data)
                    words = content[line+4].split()
                    data = list(map(float, words))
                    fluctuations.extend(data)
            elif "average conservative" in content[line] or "average dissipative" in content[line] or "average random" in content[line] or "average kinetic" in content[line] and line+4<numlines:
                words = content[line+2].split()
                data = list(map(float, words))
                averages.extend(data[0:3])
                fluctuations.extend(data[3:6])
                words = content[line+3].split()
                data = list(map(float, words))
                averages.extend(data[0:3])
                fluctuations.extend(data[3:6])
                words = content[line+4].split()
                data = list(map(float, words))
                averages.extend(data[0:3])
                fluctuations.extend(data[3:6])
                numtensor += 1
                if "conservative" in content[line]:
                    datanames += ['p_xx^c','p_xy^c','p_xz^c','p_yx^c','p_yy^c','p_yz^c','p_zx^c','p_zy^c','p_zz^c']
                elif "dissipative" in content[line]:
                    datanames += ['p_xx^d','p_xy^d','p_xz^d','p_yx^d','p_yy^d','p_yz^d','p_zx^d','p_zy^d','p_zz^d']
                elif "random" in content[line]:
                    datanames += ['p_xx^r','p_xy^r','p_xz^r','p_yx^r','p_yy^r','p_yz^r','p_zx^r','p_zy^r','p_zz^r']
                elif "kinetic" in content[line]:
                    datanames += ['p_xx^k','p_xy^k','p_xz^k','p_yx^k','p_yy^k','p_yz^k','p_zx^k','p_zy^k','p_zz^k']
            elif "average overall" in content[line] and line+4<numlines:
                totaltensor = True
                words = content[line+2].split()
                data = list(map(float, words))
                averages.extend(data[0:3])
                fluctuations.extend(data[3:6])
                words = content[line+3].split()
                data = list(map(float, words))
                averages.extend(data[0:3])
                fluctuations.extend(data[3:6])
                words = content[line+4].split()
                data = list(map(float, words))
                averages.extend(data[0:3])
                fluctuations.extend(data[3:6])
                datanames += ['p_xx','p_xy','p_xz','p_yx','p_yy','p_yz','p_zx','p_zy','p_zz']

    # if no overall pressure tensor is included, work out average
    # and fluctuations from conservative, dissipative, random and 
    # kinetic contributions (if available)

    if not totaltensor and numtensor==4:
        avetottensor = [0.0] * 9
        flutottensor = [0.0] * 9
        for i in range(9):
            avetottensor[i] += (averages[-36+i]+averages[-27+i]+averages[-18+i]+averages[-9+i])
            flutottensor[i] += math.sqrt(fluctuations[-36+i]*fluctuations[-36+i] + fluctuations[-27+i]*fluctuations[-27+i] + fluctuations[-18+i]*fluctuations[-18+i] + fluctuations[-9+i]*fluctuations[-9+i])
        averages.extend(avetottensor)
        fluctuations.extend(flutottensor)
        datanames += ['p_xx','p_xy','p_xz','p_yx','p_yy','p_yz','p_zx','p_zy','p_zz']

    return rundata,averages,fluctuations,datanames

def read_info(filename,startrun):
    """Reads header of DL_MESO_DPD OUTPUT file to find additional information about simulation (not necessarily given in CONTROL/FIELD files)"""
    
    # inputs:
    #   filename        name of OUTPUT file to start reading
    #   startrun        line in OUTPUT file where simulation data starts (stopping before this)
    # outputs:
    #   version         version number of DL_MESO_DPD used for calculation
    #   procnum         number of processor cores used for calculation
    #   threadnum       number of threads used for calculation
    #   bigendian       flag indicating whether or not big endianness was used for calculation
    #   setuptime       time taken to set up DL_MESO_DPD calculation in seconds
    #   numbeads        total number of beads used in calculation
    #   totstep         total number of timesteps for calculation
    #   equilstep       timestep number when equilibration is due to finish

    with open(filename) as file:
        content = file.read().splitlines()

    # read first lines in OUTPUT file to find version number for DL_MESO_DPD,
    # number of processor cores and threads in use, endianness, total number
    # of beads, total number of timesteps, number of equilibration timesteps, 
    # and the walltime taken to set up the calculation

    threadnum = 1
    for line in range(startrun):
        if "dl_meso version" in content[line]:
            splitstring = content[line].split("version ",1) 
            version = splitstring[1].rstrip()
        elif "RUNNING ON" in content[line]:
            splitstring = content[line].split("RUNNING ON ",1)
            if "ONE" in splitstring[1]:
                procnum = 1
            else:
                procnum = int(splitstring[1].split()[0])
        elif "THREADS/NODE" in content[line]:
            splitstring = content[line].split("WITH ",1)
            if "ONE" in splitstring[1]:
                threadnum = 1
            else:
                threadnum = int(splitstring[1].split()[0])
        elif "LITTLE ENDIAN" in content[line]:
            bigendian = False
        elif "BIG ENDIAN" in content[line]:
            bigendian = True
        elif "time elapsed since job start" in content[line]:
            splitstring = content[line].split("=",1)
            setuptime = float(splitstring[1].split()[0])
        elif "system size (particles)" in content[line]:
            splitstring = content[line].split("=",1)
            numbeads = int(splitstring[1].split()[0])
        elif "number of timesteps" in content[line]:
            splitstring = content[line].split("=",1)
            totstep = int(splitstring[1].split()[0])
        elif "equilibration period" in content[line]:
            splitstring = content[line].split("=",1)
            equilstep = int(splitstring[1].split()[0])
    
    return version,procnum,threadnum,bigendian,setuptime,numbeads,totstep,equilstep
