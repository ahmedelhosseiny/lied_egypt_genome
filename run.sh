#!/bin/bash

## If environment.yaml has been changed, the existing environment needs to be removed
## in order to re-generate the environment using:
## source ~/.bashrc; conda env remove -n workflow_lied_egypt_genome

export PATH="$HOME/miniconda3/bin:$PATH"
pwd=$PWD

mkdir -p log

echo "SETTING UP BIOCONDA ENVIRONMENT..."
source activate workflow_lied_egypt_genome
echo "RUNNING SNAKEMAKE WORKFLOW..."
snakemake -k -j 22 --jobname "{jobid}.{rulename}.sh" --cluster "sbatch --mem-per-cpu 16G -o log/%j.{rule}.log --printshellcmds --forceall"
source deactivate

conda list --export > environment_versions.yaml
