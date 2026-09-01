
# Controllo input utente
if [ -z "$1" ]; then
    echo "Uso: qsub script.sh <sample_folder>"
    exit 1
fi

# Input preso da terminale
sample_folder="$1"   # /data_folder
fastq_files="$2"  # /*fastq_prefix*

# Output
outdir="${sample_folder}/fastqc_output"
mkdir -p "$outdir"

# Load FastQC module
module load fastqc

# Esegui FastQC
fastqc --threads 10 -o "$outdir" "${sample_folder}${fastq_files}"

### HOW TO RUN
##
