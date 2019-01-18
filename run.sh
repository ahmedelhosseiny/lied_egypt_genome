#!/bin/bash

## If environment.yaml has been changed, the existing environment needs to be removed
## in order to re-generate the environment using:
#source ~/.bashrc; conda env remove -n workflow_lied_egypt_genome

mkdir -p log

echo "GENERATING AND ACTIVATING BIOCONDA WORKFLOW ENVIRONMENT..."
export PATH="$HOME/miniconda3/bin:$PATH"
if [ ! -d $HOME/miniconda3/envs/workflow_lied_egypt_genome ]; then
    conda env create -n workflow_lied_egypt_genome --file environment.yaml
fi
source activate workflow_lied_egypt_genome

echo "RUNNING SNAKEMAKE WORKFLOW..."
snakemake -n --touch --reason --rerun-incomplete -k -j 50 --resources io=4 --use-conda --jobname "{jobid}.{rulename}.sh" --cluster "sbatch --mem 16G --partition=shortterm,longterm --time 3-00:00:00 -c 12 -o log/%j.{rule}.log" --printshellcmds  align_assemblies_with_mummer_all comparison_repeatmasker compute_content_and_assembly_numbers dotplots_scaffold_vs_chromosomes_all

source deactivate
conda list -n workflow_lied_egypt_genome --export > environment_versions.yaml
