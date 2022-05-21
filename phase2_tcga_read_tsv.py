"""
This module works as single file application.
It makes part of Phase2:  process TSV files.
"""

import glob
import pandas as pd
import pickle
import os

# Phase II: reading expression data

# Getting all the filenames from TCGA folder
path = os.path.join('phase2', 'TCGA','**','*.tsv')
all_tsv_files = glob.glob(path, recursive=True)

combined_df = pd.DataFrame()
count = 0

# Reading the pickle with case id dictionary
file_dict = pickle.load(open(os.path.join('phase2', 'file_dict.p'), 'rb'))

for filename in all_tsv_files:

    # Reading each file into dataframe
    df = pd.read_csv(filename,sep='\t',index_col=0, header=1)

    # Extracting the TPM-normalized expression data
    sub_df = pd.DataFrame(df[df['gene_type'] == 'protein_coding'] \
                            ['tpm_unstranded'])

    # Renaming column names with case ids
    name = filename.split(os.path.sep)[3]
    sub_df.rename(columns={'tpm_unstranded': file_dict[name]+ '@' + name},
                  inplace=True)

    # Attaching the expression to the combined dataframe
    combined_df = combined_df.join(sub_df, how='outer')

    count += 1

# Dumping combined dataframe to pickle
path = os.path.join('phase2', 'TCGA_expressions_tpm.p')
pickle.dump(combined_df, open(path, 'wb'))
print(f'\nCombined dataframe dumped at {path}.\n')
