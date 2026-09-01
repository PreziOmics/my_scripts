# Input
adapter3="TGGAATTCTCGGGTGCCAAGG"
adapter5="AGATCGGAAGAGCGTCGTGTAGGGAAAGA"
sample_folder="path"

# Output
outdir="$sample_folder/cutadapt_output"
mkdir -p "$outdir"

# Loop sui file FASTQ
for R1_file in "$sample_folder"/*R1*.fastq.gz; do
    R2_file=${R1_file/R1/R2}  # replace R1 with R2

    R1_base=$(basename "$R1_file" .fastq.gz)
    R2_base=$(basename "$R2_file" .fastq.gz)

    R1_trimmed="${R1_base}_trimmed.fastq.gz"
    R2_trimmed="${R2_base}_trimmed.fastq.gz"

    # Cutadapt paired-end nextflex v4
    cutadapt \
        --cores 16 \
        --pair-adapters \
        --quality-cutoff 20 \
        --adapter "$adapter3" \
        -A "$adapter5" \
        --output "$outdir/$R1_trimmed" \
        --paired-output "$outdir/$R2_trimmed" \
        --minimum-length 16 \
        "$R1_file" \
        "$R2_file"
done
