

#!/bin/bash

# Define input and output paths
sample_folder="path"
mirge_lib="path"

# Output directory for miRge3.0
output_dir="$sample_folder"/mirge_output_folder

# Loop through all FASTQ files
for file in "$sample_folder"/*.fastq.gz; do
    
    base=$(basename "$file" "_cutadapt.fastq.gz")
    out_file="$base"

    echo "Processing ${file}..."

    /home/gpreziosi-iit.local/yes/bin/miRge3.0 \
        -s "$file" \
        -o "$output_dir/$out_file" \
        -on mouse \
        -db mirbase \
        -lib $mirge_lib -dex -ie -gff -cpu 20 \
       

    echo "✓ Generated: $out_file"
done

echo "All miRge output files processed."

