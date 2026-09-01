# Input
adapter3="TGGAATTCTCGGGTGCCAAGG"
sample_folder="path"

# Output
outdir="$sample_folder/cutadapt_output"
mkdir -p "$outdir"

# Loop sui file FASTQ
for file in "$sample_folder"/fastq_data/*.fastq.gz; do
    base=$(basename "$file" .fastq.gz)X

    # Cutadapt single-end
    cutadapt \
        -a "$adapter3" \
        -e 0.15 \
        --cores 20 \
        -o "$outdir/${base}_cutadapt.fastq.gz" \
        "$file"
done

