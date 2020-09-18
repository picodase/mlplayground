#!/bin/sh

# Give the job a name
#$ -N tf_gan_chemimg_grayscale_DL_mdlrun

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
module load cuda/10.1                   # load cuda
#module load cuda91/toolkit/9.1.85       # load cuda toolkit
module load python/anaconda3-5.0.0.1    # load anaconda

############################################</LOAD MODULES>###########################################

############################################<RUN COMMANDS>############################################

# prepare shell: https://stackoverflow.com/questions/55507519/python-activate-conda-env-through-shell-script
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

# INIT CONDA
#conda init bash

# ACTIVATE CONDA ENV
conda activate TF       # if already setup, no need to run functions below interactively

# INSTALL CONDA PKGS, WILL SKIP IF ALREADY THERE
#conda install -c fastai -c pytorch -c anaconda fastai gh anaconda       # fastai/pytorch
#conda install -c rdkit rdkit 

python tf_dcgan_2Dchemim_grayscale.py      # run model

mkdir data

mv image_at_epoch_* job_output.* training_checkpoints/ dcgan.gif data

conda deactivate     # deactivate

############################################</RUN COMMANDS>###########################################

############################################<UNLOAD MODULES>##########################################

module unload python/anaconda3-5.0.0.1  # unload anaconda after done using
#module unload cuda91/toolkit/9.1.85     # unload cuda
module unload cuda/10.1                   # load cuda

############################################</UNLOAD MODULES>#########################################