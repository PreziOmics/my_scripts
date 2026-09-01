#!/bin/bash
#PBS -l select=1:ncpus=24:mem=50gb:nodetype=cpu
#PBS -q workq
#PBS -l walltime=18:00:00
#PBS -N STAR_alignment

cd $PBS_O_WORKDIR

# Inizializza Conda
source /apps/miniforge/etc/profile.d/conda.sh
conda activate rna_seq

# Input
sample_folder="/work/gpreziosi/TO_ANALYZE/neurogenesis_projects/MOUSE_aged_hippo_PIWIL2"
genome_dir="/home/gpreziosi/genomes/index_STAR_mouse_GRCm39_mm39_release_M37"
suffix_name="_R1_001_u2_U2.fastq.gz"

# Output
output_dir="${sample_folder}/STAR_output_release_M37"
mkdir -p "$output_dir"

# Loop sui file R1
for file_R1 in "$sample_folder"/*_R1_001_u2_U2.fastq.gz; do
    base=$(basename "$file_R1" "_R1_001_u2_U2.fastq.gz")
    file_R2="${sample_folder}/${base}_R2_001_u2_U2.fastq.gz"
    out_prefix="${output_dir}/${base}_"

    # run STAR
    STAR --runThreadN 24 \
         --runMode alignReads \
         --genomeDir "$genome_dir" \
         --readFilesCommand zcat \
         --quantMode GeneCounts \
         --outSAMtype BAM SortedByCoordinate \
         --readFilesIn "$file_R1" "$file_R2" \
         --outFileNamePrefix "$out_prefix"
done
