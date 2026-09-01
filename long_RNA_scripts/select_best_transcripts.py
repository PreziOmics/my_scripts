'''
A single gene can produce multiple transcript isoforms.
This script processes an Excel sheet of RNA-seq data to select the single "best" (highest-confidence) transcript for each gene.
It extracts quality metrics for every transcript, filters them using a strict priority hierarchy, and exports the final list to a TSV file.

it takes as input the GTF annotation file and parses the geneName field: gene_id " ... "; gene_name " ... "; transcript_id " ... "; ecc

'''

import pandas as pd

gtf_file = pd.read_excel(path/to/file.gtf
)

gtf_file = file["geneName"]

d = {}

for annotation in genes:

    gene_name = annotation.split('gene_name "')[1].split('"')[0]
    transcript_name = annotation.split('transcript_name "')[1].split('"')[0]

    # ---------- TAG ----------
    canonical = "Ensembl_canonical" in annotation
    appris = "appris_principal" in annotation
    ccds = 'tag "CCDS"' in annotation
    basic = 'tag "basic"' in annotation

    # ---------- TSL ----------
    if "transcript_support_level" in annotation:
        tsl = annotation.split('transcript_support_level "')[1].split('"')[0]

        if tsl == "NA":
            tsl = 999
        else:
            tsl = int(tsl)
    else:
        tsl = 999

    if gene_name not in d:
        d[gene_name] = []

    d[gene_name].append(
        (
            transcript_name,
            canonical,
            appris,
            ccds,
            basic,
            tsl
        )
    )
    '''
for i, (key, value) in enumerate(d.items()):
    if i == 5:
        break
    print(key, value, len(d[key]))

'''
best_transcripts = []

for gene, transcripts in d.items(): # transcripts is a list of tuples


    # un solo trascritto
    if len(transcripts) == 1:
        best_transcripts.append(transcripts[0][0])
        continue

    # 1. Canonical
    canonical_tx = [t for t in transcripts if t[1]] #if it is canonical
    if canonical_tx:
        best = min(canonical_tx, key=lambda x: x[5]) # the best tsl (the minimum)
        best_transcripts.append(best[0])
        continue

    # 2. APPRIS
    appris_tx = [t for t in transcripts if t[2]]
    if appris_tx:
        best = min(appris_tx, key=lambda x: x[5])
        best_transcripts.append(best[0])
        continue

    # 3. CCDS
    ccds_tx = [t for t in transcripts if t[3]]
    if ccds_tx:
        best = min(ccds_tx, key=lambda x: x[5])
        best_transcripts.append(best[0])
        continue

    # 4. Basic
    basic_tx = [t for t in transcripts if t[4]]
    if basic_tx:
        best = min(basic_tx, key=lambda x: x[5])
        best_transcripts.append(best[0])
        continue

    # 5. Miglior TSL
    best = min(transcripts, key=lambda x: x[5])
    best_transcripts.append(best[0])

len(best_transcripts)

import pandas as pd

df = pd.DataFrame(best_transcripts, columns=["transcript_id"])

df.to_csv("path/of/output_file", sep="\t", index=False)

    
    
    
