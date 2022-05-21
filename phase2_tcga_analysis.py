"""
This module works as single file application.
It makes Phase2 analysis.
"""

import pickle
from scipy.stats import ttest_ind
from tqdm import tqdm
import pandas as pd
import os


# Phase II. Analysis

# List of genes of interest
LKB1 = 'ENSG00000118046.16', 'LKB1'      # LKB1 - STK11 : ENSG00000118046.16
ONECUT3 = 'ENSG00000205922.5', 'ONECUT3' # ONECUT3 : ENSG00000205922.5
AHR = 'ENSG00000106546.14', 'AHR'        # AHR : ENSG00000106546.14
ERF = 'ENSG00000105722.10', 'ERF'        # ERF : ENSG00000105722.10
NR2F6 = 'ENSG00000160113.5', 'NR2F6'     # NR2F6 : ENSG00000160113.5

genes = [LKB1,ONECUT3,AHR,ERF,NR2F6]

# Loading mutation data from the pickle
path = os.path.join('phase2', 'TCGA_maf.p')
TCGA_maf = pickle.load(open(path, 'rb'))

# Identifying cases with mutated LKB1
mutated_cases = TCGA_maf[TCGA_maf['Variant_Classification'] != 'Silent'] \
                                 ['case_id'].unique()

# Loading normalized expression data from the pickle
path = os.path.join('phase2', 'TCGA_expressions_tpm.p')
TCGA_expr_tpm = pickle.load(open(path, 'rb'))


def confirmed_f2(gene, LKB_minus, LKB_plus, abs):
    """
    Function to compare expression data.
    """
    stats = ttest_ind(LKB_minus[gene], LKB_plus[gene], nan_policy='raise',
                      alternative='greater')
    bigger = stats.pvalue < 0.05 and \
             LKB_minus[gene].mean() - LKB_plus[gene].mean() > abs
    return bigger # Returns True if the difference is statistically significant

# Getting the list of all genes
all_gene_list = TCGA_expr_tpm.index[[True]*len(TCGA_expr_tpm.index) ^ \
                TCGA_expr_tpm.index.str.startswith('N_')]
new_candidates = []
TPM_T = TCGA_expr_tpm.transpose()

# Creating empty array to store mutation status
mutated_bool_index = [False]*len(TPM_T.index)
for case_id in mutated_cases:
    mutated_bool_index = mutated_bool_index | \
                         TPM_T.index.str.startswith(case_id)

# Identifying comparison groups
LKB_minus = TPM_T[mutated_bool_index]
LKB_plus = TPM_T[ [True]*len(TPM_T.index) ^ mutated_bool_index]

# Getting all overexpressed genes
for gene in tqdm(all_gene_list): # Getting the list of overexpressed genes
    # Overexpression threshold is 11 TPM
    cf = confirmed_f2(gene,LKB_minus,LKB_plus, 11)
    if cf:
        new_candidates.append(gene)

# Reading all gene names with Ensembl IDs from one of the expression files
path = os.path.join('phase2', 'TCGA',
                    '0e9c2659-4131-4274-b0b4-3e191239f5f3',
                    '29000c50-6667-4281-a97a-5f02cb24a8f4.' +
                    'rna_seq.augmented_star_gene_counts.tsv')
gene_dict = pd.read_csv(path, sep='\t', index_col=0, header=1)['gene_name']

# Mapping gene names to the list of overexpressed Ensembl IDs
gene_names = gene_dict.transpose()[new_candidates]

# Printing overexpressed genes
print("\nOverexpressed genes: ")
for gname in gene_names:
    print(' ', gname)

print(f'\nNumber of overexpressed genes: {str(len(new_candidates))}.\n')

# Comparing candidate genes expression between groups
for gene in genes:
    cf = confirmed_f2(gene[0], LKB_minus, LKB_plus, 1)
    print(f'  Gene {gene[1]}')
    print(f'    Confirmed:            {cf}')
    print(f'    Expression in LKB1-:  {LKB_minus[gene[0]].mean()}')
    print(f'                    std:  {LKB_minus[gene[0]].std()}')
    print(f'    Expression in LKB1+:  {LKB_plus[gene[0]].mean()}')
    print(f'                    std:  {LKB_plus[gene[0]].std()}')

print("\nPhase 2 is complete.\n")
