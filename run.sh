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

snakemake -p -k -j 10 --rerun-incomplete -s Snakefile_assembly --use-conda --jobname "{jobid}.{rulename}.sh" --cluster "sbatch --mem 170G --partition=shortterm,longterm --time 6-00:00:00 -c 48 -o log/%j.{rule}.log" bwa_mapping_10x_against_grch38_all_libs map_pb_GRCh38/EGYPTREF.srt.bam
#snakemake -n --reason --rerun-incomplete -k -j 8 --resources io=4 --use-conda --jobname "{jobid}.{rulename}.sh" --cluster "sbatch --mem 90G --partition=shortterm,longterm --time 3-00:00:00 -c 16 -o log/%j.{rule}.log" --printshellcmds # dotplots_scaffold_vs_chromosomes_all results/repeatmasker_comparison.txt compute_content_and_assembly_numbers comparison_repeatmasker dotplots_scaffold_vs_chromosomes_all

source deactivate
conda list -n workflow_lied_egypt_genome --export > environment_versions.yaml
