#!/bin/bash
#PBS -l select=1:ncpus=24:mem=50gb:nodetype=cpu
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -N Salmon_alignment

cd $PBS_O_WORKDIR

# Inizializza Conda
source /apps/miniforge/etc/profile.d/conda.sh
conda activate salmon_env

# Input
sample_folder="path"
transcriptome_index="path"
suffix_name="_R1_001_u2_U2.fastq.gz"

# Output
output_dir="${sample_folder}/salmon_quant"
mkdir -p "$output_dir"

# Lista campioni
sample_list=(S1 S2 SX)

# Loop
for sample in "${sample_list[@]}"; do
  for file_R1 in "${sample_folder}"/*"${sample}"*"${suffix_name}"; do

    # Controllo di sicurezza: salta se il file non esiste fisicamente
    [ -e "$file_R1" ] || continue

    # Estrazione pulita del nome del campione
    base=$(basename "$file_R1" "$suffix_name")
    
    # Costruzione corretta del file R2 associato
    file_R2="${sample_folder}/${base}_R2_001_u2_U2.fastq.gz"
    out_dir="${output_dir}/${base}"

    mkdir -p "$out_dir"

    echo "--------------------------------------------------"
    echo "Processing sample: $base"
    echo "--------------------------------------------------"

    # Esecuzione Salmon quant con 24 core (pieno utilizzo di ncpus=24)
    salmon quant \
      --libType ISR \
      --index "$transcriptome_index" \
      -1 "$file_R1" \
      -2 "$file_R2" \
      --numBootstraps 10 \
      --threads 24 \
      -o "$out_dir"

  done
done
