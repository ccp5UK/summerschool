from pathlib import Path
import subprocess
import shlex
import os
from tqdm.auto import tqdm

def run_MCLB(lbeexe='LBE', numstep=4000, description=''):
    """
    Launches CCP5 multicomponent LBE (MCLB) code and displays a progress bar as it runs.

    #Parameters:
    #    lbeexe (str): Location of LBE code executable (usually called LBE)
    #    numstep (int): Number of timesteps to be carried out by LBE code
    #    description (str): Optional description of calculation placed at beginning of line with progress bar
    
    """
    command = str(Path(lbeexe).resolve())
    # prepare output file to hold standard output and determine current number of 
    outfile = "lbecount"
    outpipe = open(outfile, 'w')
    # run LBE code in directory and set up progress barout
    lberun = subprocess.Popen(shlex.split(command), stdout=outpipe)
    stepnum0 = 0
    pbar = tqdm(total=numstep, desc=description)
    terminate = False
    # keep checking LBE code is running until it stops (according to piped output file)
    if lberun.poll() is None:
        while stepnum0<numstep:
            with open(outfile, 'r') as file:
                lines = file.readlines()
                stepnum = sum(1 for line in lines if "time" in line)
                if stepnum != stepnum0:
                    pbar.update(stepnum-stepnum0)
                    stepnum0 = stepnum
    # close down progress bar
    pbar.close()
    # delete file (we no longer need it!)
    os.remove(outfile)
