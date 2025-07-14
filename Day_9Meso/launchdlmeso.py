from pathlib import Path
import subprocess
import shlex
import os
import glob
import dlmoutputread as dlmo
from tqdm.auto import tqdm

def run_DPD(rundir='DPD1Ex2', dlmesoexe='dl_meso/WORK/dpd.exe', numcores=1, deleteold=True, description=''):
    """
    Launches DL_MESO_DPD in a user-selected directory and displays a progress bar as it runs.

    #Parameters:
    #    rundir (str): Path to directory where DL_MESO_DPD will run
    #    dlmesoexe (str): Location of DL_MESO_DPD executable (usually called dpd.exe)
    #    numcores (int): Number of processor cores to run DL_MESO_DPD (default: 1)
    #    deleteold (boolean): Delete any output files from old calculation (default: True)
    #    description (str): Optional description of calculation placed at beginning of line with progress bar
    
    """
    outdir = str(Path(rundir).resolve())
    command = str(Path(dlmesoexe).resolve())
    if numcores>1:
        command = "mpirun -np {0:d} ".format(numcores) + command
    # delete any old output files if they exist (and we want to do this!)
    if deleteold:
        if os.path.exists(outdir+"/CORREL"):
            os.remove(outdir+"/CORREL")
        if os.path.exists(outdir+"/export"):
            os.remove(outdir+"/export")
        if os.path.exists(outdir+"/REVIVE"):
            os.remove(outdir+"/REVIVE")
        if os.path.exists(outdir+"/HISTORY"):
            os.remove(outdir+"/HISTORY")
        if os.path.exists(outdir+"/OUTPUT"):
            os.remove(outdir+"/OUTPUT")
    # get hold of number of timesteps for DL_MESO_DPD calculation
    # from CONTROL file in supplied directory
    with open(outdir+"/CONTROL", 'r') as control:
        lines = control.readlines()
        for line in lines:
            words = line.replace(',',' ').replace('\t',' ').lower().split()
            if len(words)>0:
                if words[0].startswith('steps') and len(words)>1:
                    numstep = int(words[1])
    # prepare OUTPUT file for DL_MESO_DPD to write to
    # (by piping standard output into it)
    outfile = outdir+"/OUTPUT"
    outpipe = open(outfile, 'w')
    # run DL_MESO_DPD in directory and set up progress barout
    dlmesorun = subprocess.Popen(shlex.split(command), cwd=outdir, stdout=outpipe)
    stepnum0 = 0
    pbar = tqdm(total=numstep, desc=description)
    terminate = False
    # keep checking DL_MESO_DPD is running until it stops (according to OUTPUT file)
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
    # close down progress bar
    pbar.close()

def run_LBE(rundir='LBE2Ex4', dlmesoexe='dl_meso/WORK/lbe.exe', numcores=1, deleteold=True, description=''):
    """
    Launches DL_MESO_LBE in a user-selected directory and displays a progress bar as it runs.

    #Parameters:
    #    rundir (str): Path to directory where DL_MESO_LBE will run
    #    dlmesoexe (str): Location of DL_MESO_LBE executable (usually called lbe.exe)
    #    numcores (int): Number of processor cores to run DL_MESO_LBE (default: 1)
    #    deleteold (boolean): Delete any output files from old calculation (default: True)
    #    description (str): Optional description of calculation placed at beginning of line with progress bar
    
    """
    outdir = str(Path(rundir).resolve())
    command = str(Path(dlmesoexe).resolve())
    if numcores>1:
        command = "mpirun -np {0:d} ".format(numcores) + command
    # delete any old output files if they exist (and we want to do this!)
    if deleteold:
        if os.path.exists(outdir+"/lbout.dump"):
            os.remove(outdir+"/lbout.dump")
        if os.path.exists(outdir+"/lbout.screen"):
            os.remove(outdir+"/lbout.screen")
        if os.path.exists(outdir+"/lbout000000.vts"):
            for f in glob.glob(outdir+"/*.vts"):
                os.remove(f)
    # get hold of number of timesteps for DL_MESO_LBE calculation
    # from lbin.sys file in supplied directory
    with open(outdir+"/lbin.sys", 'r') as control:
        lines = control.readlines()
        for line in lines:
            words = line.replace(',',' ').replace('\t',' ').split()
            if len(words)>0:
                if words[0].startswith('total_step') and len(words)>1:
                    numstep = int(words[1])
    # prepare lbout.screen output file for DL_MESO_LBE to write to
    # (by piping standard output into it)
    outfile = outdir+"/lbout.screen"
    outpipe = open(outfile, 'w')
    # run DL_MESO_LBE in directory and set up progress barout
    dlmesorun = subprocess.Popen(shlex.split(command), cwd=outdir, stdout=outpipe)
    stepnum0 = 0
    pbar = tqdm(total=numstep, desc=description)
    terminate = False
    # keep checking DL_MESO_LBE is running until it stops (according to lbout.screen file)
    if dlmesorun.poll() is None:
        while not terminate:
            if os.path.getsize(outfile)>0:
                stepnum = 0
                with open(outfile, 'r') as file:
                    lines = file.readlines()
                    for line in lines:
                        if "MASS: total" in line:
                            words = line.split()
                            stepnum = max(stepnum, int(words[0]))
                        if "Program finished" in line:
                            terminate = True
                if stepnum > stepnum0:
                    pbar.update(stepnum-stepnum0)
                    stepnum0 = stepnum
    # close down progress bar
    pbar.close()
