#!/bin/bash
#SBATCH --job-name=gcsi_pipeline
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH --time 4-23:59:00
#SBATCH --chdir /cluster/projects/bhklab/projects/gCSI_human/gCSI_PharmacoSet-Pipeline
#SBATCH --output /cluster/projects/bhklab/projects/gCSI_human/gCSI_PharmacoSet-Pipeline/logs/%x_%j.out
source ~/.bashrc
conda activate snakemake 
snakemake --profile workflow/profiles/
