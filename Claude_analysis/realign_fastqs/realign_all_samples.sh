#!/bin/bash
# Realign every library in one array job.
#
# Replaces the four per-condition scripts. The reason for the rewrite is a bug in
# them: for the rep2 libraries the collaborator passed CellRanger *two* --fastqs
# paths (the original 20251117 run and the 20251201 resequencing run), and the old
# scripts passed only the resequencing run. Control_2 therefore aligned 161M reads
# against their 326M, and its count matrix does not match theirs -- while
# Control_1, Control_1A and Control_3, whose inputs did match, reproduce their
# matrices bit for bit. Heat_2, Drought_2 and HD_2 have the same defect.
#
# The collaborator's own _invocation for RNA_leaf_C2_merged_GEX lists:
#   Rep_2/GR0073_20251117_correct_demultiplex_addedRNA/GR0073
#   Rep_2/GR0073_RNA_resequenced_20251201_/GR0073
# (hence the "_merged" in all four rep2 output directory names).
#
# NOTE: our mirror nests the reads one level deeper, under `non_trimmed`. The
# 20251117 path below follows that convention but has not been verified from this
# machine -- the preflight check will abort loudly if it is wrong, and print what
# is actually there so it can be corrected in one place.

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --array=1-13
#SBATCH --time=10:00:00
#SBATCH --mem=80GB
#SBATCH --job-name=soybean_realign
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=michael.p@nyu.edu
#SBATCH --account=torch_pr_121_general
#SBATCH --output=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign.out
#SBATCH --error=/home/mp7563/Jobs_logs/Soybean_jobs/%A_%a_soybean_realign.err

set -euo pipefail
stdbuf -o0 echo "Script Started"

RAW=/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Raw_fastq_files
OUT_BASE=/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Soybean_results/realigned_sc
TRANSCRIPTOME=/projects/rps/cgsb/bergelson/bergelson-lab/Michael_P/Collaborator_data/Genomes/soy_Williams82_v3_with_organelles

# The two rep2 sequencing runs that must be passed together.
REP2_RESEQ="$RAW/Rep_2/GR0073_RNA_resequenced_20251201_/GR0073/non_trimmed"
REP2_ORIG="$RAW/Rep_2/GR0073_20251117_correct_demultiplex_addedRNA/GR0073/non_trimmed"
REP2_BOTH="$REP2_ORIG,$REP2_RESEQ"

REP3="$RAW/Rep_3/GR0192_RNA/RNA_fastq_raw"

#            run_id       condition      sample_prefix          fastqs
SAMPLES=(
  "Control_1   |control      |RNA_leaf_C1          |$RAW/GR00_Rep1/Leaf_Control1/RNA_leaf_C1"
  "Control_1A  |control      |RNA_leaf_C1A         |$RAW/GR00_Rep1/Leaf_control1A/RNA_leaf_C1A"
  "Control_2   |control      |RNA-Leaf-C2          |$REP2_BOTH"
  "Control_3   |control      |RNA-leaf-control_3   |$REP3"
  "Heat_1      |heat         |RNA_leaf_H1          |$RAW/GR00_Rep1/Leaf_heat1/RNA_LeafH1"
  "Heat_2      |heat         |RNA-Leaf-H2          |$REP2_BOTH"
  "Heat_3      |heat         |RNa-leaf-heat_3      |$REP3"
  "Drought_1   |drought      |RNA_leaf_D1          |$RAW/GR00_Rep1/Leaf_drought1/RNA_leaf_D1"
  "Drought_2   |drought      |RNA-Leaf-D2          |$REP2_BOTH"
  "Drought_3   |drought      |RNA-leaf-drought_3   |$REP3"
  "HD_1        |heat_drought |RNA_leaf_HD1         |$RAW/GR00_Rep1/Leaf_HD1/RNA_leaf_HD1"
  "HD_2        |heat_drought |RNA-Leaf-HD2         |$REP2_BOTH"
  "HD_3        |heat_drought |RNA-leaf-HD_3        |$REP3"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
IFS='|' read -r RUN_ID CONDITION SAMPLE_PREFIX FASTQS <<< "${SAMPLES[$IDX]}"
RUN_ID=$(echo "$RUN_ID" | xargs)
CONDITION=$(echo "$CONDITION" | xargs)
SAMPLE_PREFIX=$(echo "$SAMPLE_PREFIX" | xargs)
FASTQS=$(echo "$FASTQS" | xargs)

stdbuf -o0 echo "run_id=$RUN_ID condition=$CONDITION sample=$SAMPLE_PREFIX"
stdbuf -o0 echo "fastqs=$FASTQS"

# Preflight: every --fastqs path must exist and hold reads for this sample prefix.
# Cheap here, and it turns a silently-half-aligned library into an immediate failure.
IFS=',' read -ra FASTQ_DIRS <<< "$FASTQS"
for d in "${FASTQ_DIRS[@]}"; do
  if [[ ! -d "$d" ]]; then
    echo "ERROR: fastq dir does not exist: $d" >&2
    echo "Contents of its parent:" >&2
    ls -la "$(dirname "$d")" >&2 || true
    exit 1
  fi
  n=$(find "$d" -maxdepth 1 -name "${SAMPLE_PREFIX}_*.fastq.gz" | wc -l)
  if [[ "$n" -eq 0 ]]; then
    echo "ERROR: no ${SAMPLE_PREFIX}_*.fastq.gz in $d" >&2
    ls -la "$d" >&2 || true
    exit 1
  fi
  echo "  ok: $n fastq files for $SAMPLE_PREFIX in $d"
done

mkdir -p "$OUT_BASE/$CONDITION"
cd "$OUT_BASE/$CONDITION" || exit 1

singularity exec --nv \
    --overlay /home/mp7563/Python_envs/Torch_miniforge_env.ext3:ro /home/mp7563/Python_envs/ubuntu_cuda_image.sif \
    /bin/bash -c "source /ext3/env.sh; conda activate Luke_terrace; \
    cellranger count \
    --id=$RUN_ID \
    --transcriptome=$TRANSCRIPTOME \
    --fastqs=$FASTQS \
    --sample=$SAMPLE_PREFIX \
    --chemistry=ARC-v1 \
    --localcores=24 \
    --create-bam=true \
    --localmem=75"
