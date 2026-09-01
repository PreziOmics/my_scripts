
genome_fasta="path"
annotation_gtf="path"
genome_dir="path"
# copy input data to requested scratch space
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $genome_fasta --sjdbGTFfile $annotation_gtf 

# copy back output data to the submission directorye
