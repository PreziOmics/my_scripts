#### FASTQ-DUMP recupera gli SRA di ogni accession e lo converte in formato FASTQ
# Input
sample_folder="path"
file_list=SRR_Acc_List.txt

#output

# Loop sui file 
while read line; do

    fastq-dump --gzip --outdir $sample_folder $line ;
    
done < $sample_folder/$file_list

## n.B split-3 puo essere usato sia per dati single end che paired, questo perche se non trova l'altra coppia biologica non assegna il _1 e _2
