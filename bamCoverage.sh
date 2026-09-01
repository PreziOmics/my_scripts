#!/bin/bash

# Define input and output paths
sample_folder="path"

# Effective genome sizes
GRCm39=2654621783   # mouse
GRCh38=2913022398   # human


# Output directory for bigWig tracks
output_dir="${sample_folder}/bigWig_normalized_tracks"

mkdir -p $output_dir 


# Loop through all BAM files
for file in "${sample_folder}"/*sorted.bam; do
    base=$(basename "$file" .bam)
    output_file="${output_dir}/${base}_normalized_BPM.bw"

    echo "Processing ${file}..."
    
    bamCoverage \
        --bam "$file" \
        -o "$output_file" \
        --outFileFormat bigwig \
        --binSize 50 \
        --normalizeUsing BPM \
        --effectiveGenomeSize "$GRCm39" \
        -p 20

    echo "✓ Generated: ${output_file}"
done

echo "All BAM files processed. Normalized bigWig files saved in: $output_dir"

