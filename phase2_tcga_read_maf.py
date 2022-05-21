"""
This module works as single file application.
It makes part of Phase2:  process MAF files.
"""

import glob
import pandas as pd
import pickle
import os
import gzip

# Phase II: reading expression data
# Getting all the filenames from TCGA-MAF folder
path = os.path.join('phase2', 'TCGA_MAF','**','*.gz')
all_maf_files = glob.glob(path, recursive=True)

combined_df = pd.DataFrame()
count = 0
for filename in all_maf_files:

    # Reading each file into dataframe
    with gzip.open(filename) as f:
        df = pd.read_csv(filename,sep='\t',index_col=False, header=1, skiprows=6)

    # Extracting the mutation data for LKB1
    filtered_df = df[df['Entrez_Gene_Id']==6794]

    if filtered_df.shape[0]>0:
        # Attaching the mutation to the combined dataframe
        combined_df = pd.concat([combined_df,filtered_df])
        count += 1

# Dumping combined dataframe to pickle
path = os.path.join('phase2', 'TCGA_maf.p')
pickle.dump(combined_df, open(path, 'wb'))
print(f'\nCombined dataframe dumped at {path}.\n')
