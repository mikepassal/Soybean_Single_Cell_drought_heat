#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=47:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=soybean_realign
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=michael.p@nyu.edu
#SBATCH --account=torch_pr_121_general
#SBATCH --output=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign.out
#SBATCH --error=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign.err


cellranger count \
    --id = run_control_samples \
    --transcriptome = /projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Genomes/soy_Williams82_v3_with_organelles \ 
    --fastqs = /projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Soybean_Single_Cell_drought_heat/Cellranger_output/control \