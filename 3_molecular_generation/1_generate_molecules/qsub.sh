#!/bin/bash
#PBS -N drugex-led3
#PBS -q cheminf
#PBS -l select=1:ncpus=2:ngpus=3:mem=8gb
#PBS -l walltime=720:00:00
#PBS -m ae

set -e

# set important variables
export WORKDIR=${WORKDIR:-`pwd`}
export DRUGEX_SCORERS=${DRUGEX_SCORERS:-"sascore,qsar_cls"}
export DRUGEX_JOB=${DRUGEX_JOB:-"drugex-led3"}

# work begins here
export SCRATCHDIR=/scratch/$USER/$PBS_JOBID # this might be useful -> we can get this variable from python and fetch big data there
mkdir $SCRATCHDIR

# append a line to a file "jobs_info.txt" containing the ID of the job and the current worker hostname
echo "$PBS_JOBID is running on node `hostname -f`." >> $WORKDIR/jobs_info.txt

# go to the working directory and copy over files
cp -r $HOME/projects/led3_score $SCRATCHDIR/led3_score
cd $SCRATCHDIR/led3_score/3_molecular_generation/1_generate_molecules

# activate the conda environment and save version info
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate led3-drugex

python RL_run.py

cp -vr $SCRATCHDIR/led3_score/3_molecular_generation/1_generate_molecules/models/RL/$DRUGEX_JOB $HOME/projects/led3_score/3_molecular_generation/1_generate_molecules/models/RL/$DRUGEX_JOB && rm -r $SCRATCHDIR
