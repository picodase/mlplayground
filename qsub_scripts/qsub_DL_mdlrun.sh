#!/bin/sh

# Give the job a name
#$ -N DL_mdlrun

# set the shell
#$ -S /bin/sh

# set working directory on all host to
# directory where the job was started
#$ -cwd

# send all process STDOUT (fd 2) to this file
#$ -o job_output.txt

# send all process STDERR (fd 3) to this file
#$ -e job_output.err

# email information
#$ -m e
 
# Just change the email address. You will be emailed when the job has finished.
#$ -M northj@oregonstate.edu

# generic parallel environment with 2 cores requested
#$ -pe orte 2

############################################<LOAD MODULES>############################################

# Load a module, if needed
module load python/anaconda3-5.0.0.1    # load anaconda

############################################</LOAD MODULES>###########################################

############################################<RUN COMMANDS>############################################

# INIT CONDA
conda init bash

# ACTIVATE CONDA ENV
conda activate ML       # if already setup, no need to run functions below interactively

# INSTALL CONDA PKGS, WILL SKIP IF ALREADY THERE
#conda install -c fastai -c pytorch -c anaconda fastai gh anaconda       # fastai/pytorch
#conda install -c rdkit rdkit 

python mdlrun_gan.py      # run model

conda deactivate ML     # deactivate

############################################</RUN COMMANDS>###########################################

############################################<UNLOAD MODULES>##########################################

module unload python/anaconda3-5.0.0.1  # unload anaconda after done using

############################################</UNLOAD MODULES>#########################################