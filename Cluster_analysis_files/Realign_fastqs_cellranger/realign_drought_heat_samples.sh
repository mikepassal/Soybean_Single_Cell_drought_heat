#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --array=1-3
#SBATCH --time=10:00:00
#SBATCH --mem=80GB
#SBATCH --job-name=soybean_realign_HD   
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=michael.p@nyu.edu
#SBATCH --account=torch_pr_121_general
#SBATCH --output=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign_HD.out
#SBATCH --error=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign_HD.err

stdbuf -o0 echo "Script Started"

SAMPLE_LABELS=(
  "RNA_leaf_HD1"
  "RNA-Leaf-HD2"
  "RNA-leaf-HD_3"
)

RUN_IDS=(
  "HD_1"
  "HD_2"
  "HD_3"
)

FASTQ_LIST=(
  "/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Raw_fastq_files/GR00_Rep1/Leaf_HD1/RNA_leaf_HD1"
  "/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Raw_fastq_files/Rep_2/GR0073_RNA_resequenced_20251201_/GR0073/non_trimmed"
  "/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Raw_fastq_files/Rep_3/GR0192_RNA/RNA_fastq_raw"
)

IDX=$((SLURM_ARRAY_TASK_ID-1))
RUN_ID="${RUN_IDS[$IDX]}"
SAMPLE_PREFIX="${SAMPLE_LABELS[$IDX]}"
FASTQS="${FASTQ_LIST[$IDX]}"

mkdir -p /projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Soybean_results/realigned_sc/heat_drought
cd /projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Soybean_results/realigned_sc/heat_drought || exit 1

singularity exec --nv \
    --overlay /home/mp7563/Python_envs/Torch_miniforge_env.ext3:ro /home/mp7563/Python_envs/ubuntu_cuda_image.sif \
    /bin/bash -c "source /ext3/env.sh; conda activate Luke_terrace; \
    cellranger count \
    --id="$RUN_ID" \
    --transcriptome=/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Genomes/soy_Williams82_v3_with_organelles \
    --fastqs="$FASTQS" \
    --sample="$SAMPLE_PREFIX" \
    --chemistry=ARC-v1 \
    --localcores=24 \
    --create-bam=true \
    --localmem=75"