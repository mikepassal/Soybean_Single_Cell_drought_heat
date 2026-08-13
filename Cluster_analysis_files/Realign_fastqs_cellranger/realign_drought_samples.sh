#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --array=2
#SBATCH --time=1:59:00
#SBATCH --mem=80GB
#SBATCH --job-name=soybean_realign_drought
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=michael.p@nyu.edu
#SBATCH --account=torch_pr_121_general
#SBATCH --output=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign_drought.out
#SBATCH --error=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign_drought.err

stdbuf -o0 echo "Script Started"

SAMPLE_LABELS=(
  "RNA_leaf_D1"
  "RNA-Leaf-D2"
  "RNA-leaf-drought_3"
)

RUN_IDS=(
  "Drought_1"
  "Drought_2"
  "Drought_3"
)

FASTQ_LIST=(
  "/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Raw_fastq_files/GR00_Rep1/ALL_RNA"
  "/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Raw_fastq_files/Rep_2/GR0073_20251117_correct_demultiplex_addedRNA/GR0073/untrimmed_MICHAEL_CREATED_for_confirmation"
  "/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Raw_fastq_files/Rep_3/GR0192_RNA/RNA_fastq_raw"
)

IDX=$((SLURM_ARRAY_TASK_ID-1))
RUN_ID="${RUN_IDS[$IDX]}"
SAMPLE_PREFIX="${SAMPLE_LABELS[$IDX]}"
FASTQS="${FASTQ_LIST[$IDX]}"

mkdir -p /projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Soybean_results/realigned_sc/drought
cd /projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Soybean_results/realigned_sc/drought || exit 1

singularity exec --nv \
    --overlay /home/mp7563/Python_envs/Torch_miniforge_env.ext3:ro /home/mp7563/Python_envs/ubuntu_cuda_image.sif \
    /bin/bash -c "source /ext3/env.sh; conda activate Luke_terrace; \
    cellranger count \
    --id="$RUN_ID" \
    --transcriptome=/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Genomes/soy_Williams82_v3_with_organelles \
    --fastqs="$FASTQS" \
    --sample="$SAMPLE_PREFIX" \
    --chemistry=ARC-v1 \
    --localcores=32 \
    --create-bam=true \
    --localmem=75"