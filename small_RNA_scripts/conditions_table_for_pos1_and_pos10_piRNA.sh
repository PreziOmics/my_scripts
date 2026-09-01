import pandas as pd 
from Bio import SeqIO # import Biopython

fastq=open("path/to/file.fastq")

####### import the fastq in a dictionary, only the header and sequence
fastq_dict={}
for record in SeqIO.parse(fastq, "fastq"):
    fastq_dict[record.id]=str(record.seq)
fastq.close()

############### convert dictionary in a dataframe (df)
df = pd.Series(fastq_dict, name="sequence").rename_axis("header").reset_index()
############### initialize a dataframe (matrix) with zeros
conditions=["1A","1B","1C","1D","2A","2B","2C","2D"]
df[conditions]=0        

#### get the position 1 and 10  of each sequence
pos1=df["sequence"].str[0]
pos10=df["sequence"].str[9]
substring=pos1+pos10  
pos1_no_T = substring.str[0] != "T"

#### verify the substring for each sequence:
#### True if sutisfied, otherwise False. astype(int) convert Falso to 0 and True to 1
df['1A'] = (substring == "TA").astype(int) #
df['1B'] = (substring == "TC").astype(int)
df['1C'] = (substring == "TG").astype(int)
df['1D'] = (substring == "TT").astype(int)
df['2A'] = ((pos10 == "A") & pos1_no_T).astype(int)
df['2B'] = ((pos10 == "C") & pos1_no_T).astype(int)
df['2C'] = ((pos10 == "G") & pos1_no_T).astype(int)
df['2D'] = ((pos10 == "T") & pos1_no_T).astype(int)

df.to_csv("/home/gpreziosi-iit.local/Downloads/conditions_table.tsv",sep="\t",index = None)
