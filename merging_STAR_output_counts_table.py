
# merging STAR output count tables of multiple samples

import pandas as pd
import os 


folder="path"

samples=pd.read_csv(f"{folder}/sample_table.txt",names=["file","sample"],sep="\t")
conta=0
for file in os.listdir(folder): 
    if "ReadsPerGene" in file: 
        sample=samples.loc[samples["file"]==file,"sample"].iloc[0]
        if conta == 0: 
            merged_counts=pd.read_csv(folder+"/"+file,sep="\t",names=["gene_id","unstranded","sforward",sample]).loc[4:,["gene_id",sample]]
            conta+=1
        
        else: 
            cur_counts=pd.read_csv(folder+"/"+file,sep="\t",names=["gene_id","unstranded","sforward",sample]).loc[4:,["gene_id",sample]]
            merged_counts=pd.merge(merged_counts, cur_counts, on=["gene_id"], how="outer")
        
merged_counts.head()

merged_counts.to_csv(folder+"/raw_counts.tsv",sep="\t",index=None)
        
        
    
