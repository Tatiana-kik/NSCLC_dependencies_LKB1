"""
This module works as single file application.
It makes Phase1 DepMap analysis.
"""

import warnings
import time
from tqdm import tqdm
from scipy.stats import ttest_ind
import phase1_dataload_module as dl
import pickle
import os


# Phase I: Analyzing DepMap data

# Constants
LKB1_entrez = 6794             # Entrez gene ID for LKB1. Used in some datasets
LKB1_str = 'STK11 (6794)'      # ID for LKB1. Used in some datasets
SENS_THRESHOLD = 0.0           # Threshold for knockout sensitivity
SENS_LH_THRESHOLD = 0.0        # Likelihood threshold for LKB1- cell lines
SENS_LH_THRESHOLD_OTHER = 0.2  # Likelihood threshold for LKB1+ cell lines
EXPR_THRESHOLD = 3             # LKB1 expression threshold
T_TEST_P_THRESHOLD = 0.05      # P-value threshold

def data_pre_load():
    """
    Function to read the data from the disk and dump it into pickle
    for faster debugging.
    """
    start_time = time.time()

    # Getting all NSCLC cell lines

    # Take cell line IDs for lung cancer
    cell_lines_Lung = dl.get_cell_lines_for_disease('lineage',
                                                    'lung')

    # Take cell line IDs for NSCLC
    cell_lines_NSCLC = dl.get_cell_lines_for_disease('lineage_subtype',
                                                     'NSCLC')

    # Getting NSCLC cell lines with damaging mutations
    badly_mutated_NSCLC_LKB1 = dl.get_mutated_cell_lines(cell_lines_NSCLC,
                                                         LKB1_entrez)

    # Getting LKB1 knockout-insensitive cell lines with thresholds on
    # sensitivity and likelihood
    LKB1_ko_insensitive = dl.get_cl_with_ko_lh(LKB1_str,
                                               SENS_THRESHOLD,
                                               SENS_LH_THRESHOLD)

    # Getting cell lines with underexpressed LKB1 (expression < threshold)
    LKB1_underexpressed = dl.get_cl_expr_under_td(gene_str=LKB1_str,
                                                  threshold=EXPR_THRESHOLD)

    expressions = dl.read_expression()
    mutations = dl.read_mutations()
    gene_dependency = dl.read_gene_dependency()
    gene_effect =  dl.read_gene_effect()

    # Dumping read date into pickle
    data_dump = (cell_lines_Lung,
                 cell_lines_NSCLC,
                 badly_mutated_NSCLC_LKB1,
                 LKB1_ko_insensitive,
                 LKB1_underexpressed,
                 expressions,
                 mutations,
                 gene_dependency,
                 gene_effect)
    path = os.path.join('phase1', 'data_dump.p')
    pickle.dump(data_dump, open(path, 'wb'))
    print("Time for data pre-load: " +
          str(round(time.time()-start_time)) + ' seconds')


if __name__ == '__main__':

    # Comment this line to restart the analysis if pickle data dump
    # is made before
    data_pre_load()

    start_time = time.time()

    path = os.path.join('phase1', 'data_dump.p')

    (cell_lines_Lung,
     cell_lines_NSCLC,
     badly_mutated_NSCLC_LKB1,
     LKB1_ko_insensitive,
     LKB1_underexpressed,
     expressions,
     mutations,
     gene_dependency,
     gene_effect) = pickle.load(open(path, "rb"))

    print("Time for pickle load: " +
          str(round(time.time()-start_time)) + ' seconds')

    # Setting first group of cell lines for comparison (LKB1-)
    LKB1_loss_cell_lines = set(cell_lines_NSCLC) & \
                           set(LKB1_ko_insensitive) & \
                           (set(badly_mutated_NSCLC_LKB1) | \
                            set(LKB1_underexpressed))

    # Setting second group of cell lines for comparison (LKB1+)
    Other_cell_lines = set(cell_lines_Lung) - (set(cell_lines_NSCLC) |
                                               set(LKB1_ko_insensitive) |
                                               set(badly_mutated_NSCLC_LKB1) |
                                               set(LKB1_underexpressed))

    # Setting list of genes to identify candidates

    # Getting full list of genes from column names of the Expressions file
    genes = set(expressions.columns) - set(['Unnamed: 0'])
    count = 0                                       # Candidate counter
    pbar = tqdm(genes)                              # Progress bar
    pbar.set_description(f'Genes found : {count} ') # Progress bar description
    genes_found = []                                # List to store candidates

    for gene in pbar:
        # Get Entrez gene ID from column name of Expressions file
        gene_entrez = int(gene.split("(")[1].split(")")[0])

        # Get cell lines with all mutations
        badly_mutated_this_gene_raw = mutations[mutations['Entrez_Gene_Id'] == gene_entrez]

        # Get cell lines with damaging mutations
        damaging = badly_mutated_this_gene_raw['Variant_annotation'] == 'damaging'
        badly_mutated_this_gene = badly_mutated_this_gene_raw[damaging]['DepMap_ID']

        # Exclude cell lines with damaging mutations from LKB1+ group
        Other_cell_lines_clean = set(Other_cell_lines) - set(badly_mutated_this_gene)

        # Exclude cell lines with damaging mutations from LKB1- group
        LKB1_loss_cell_lines_clean = set(LKB1_loss_cell_lines) - set(badly_mutated_this_gene)

        # Is the gene statistically significantly overexpressed?
        expr_in_lost = expressions[expressions['Unnamed: 0'].
                       isin(LKB1_loss_cell_lines_clean)][gene]
        explr_in_other = expressions[expressions['Unnamed: 0'].
                         isin(Other_cell_lines_clean)][gene]

        # Student's t-test:
        t_expr = ttest_ind(expr_in_lost, explr_in_other, nan_policy='raise',
                           alternative='greater').pvalue < T_TEST_P_THRESHOLD

        # Is the cell line more knockout-sensitive for this gene?
        with warnings.catch_warnings(): # Catch Scipy warnings
            warnings.filterwarnings('error')
            try:
                if gene in gene_effect.columns:
                    # Knockout sensitivity for LKB1+ cell lines
                    sens_in_other = gene_effect[gene_effect['DepMap_ID'].
                                    isin(Other_cell_lines_clean)]. \
                                    loc[:,['DepMap_ID',gene]]
                    # For knockout sensitivity in LKB1- cell lines we filter
                    # the knockout likelihood (adjustment for CRISPR-Cas9
                    # toxicity effect)
                    sens_in_lost = gene_effect[gene_effect['DepMap_ID'].
                                   isin(LKB1_loss_cell_lines_clean)]. \
                                   loc[:,['DepMap_ID',gene]]
                    sens_lh_in_lost = gene_dependency[gene_dependency['DepMap_ID'].
                                      isin(LKB1_loss_cell_lines_clean)]. \
                                      loc[:,['DepMap_ID',gene]]
                    sens_lh_in_lost.rename(columns = {gene: 'lh'})
                    joint_sens_in_lost = sens_in_lost.set_index('DepMap_ID'). \
                                         join(sens_lh_in_lost.set_index('DepMap_ID'),
                                              lsuffix='_sensitivity',
                                              rsuffix='_likelihood')

                    for_tt_lost = joint_sens_in_lost[joint_sens_in_lost[gene+'_likelihood'] > \
                                  SENS_LH_THRESHOLD_OTHER][gene+'_sensitivity']
                    for_tt_other = sens_in_other[gene]

                    if len(for_tt_lost) > 1 and len(for_tt_other) > 1:
                        # Student's t-test for knockout sensitivity
                        t_ko = ttest_ind(for_tt_lost,
                                         for_tt_other,
                                         nan_policy='raise',
                                         alternative='less'). \
                               pvalue < T_TEST_P_THRESHOLD
                        # If cell growth is statistically significantly less
                        #    in the in LKB1- than in LKB1+
                        # AND cell growth is positive in LKB1+
                        # AND cell growth is negative in LKB1-
                        # THEN one of the criteria is met
                        t_ko_f = t_ko and for_tt_lost.mean() < 0 \
                                      and for_tt_other.mean() > 0
                    else:
                        # We don't take the gene if there is not enough data
                        # for Student t-test analysis
                        t_ko_f = False
                else:
                    # We don't take the gene if there is no data of the
                    # knockout sensitivity
                    t_ko_f = False
            except RuntimeWarning as e:
                # Print an error message in case of any Scipy warnings
                print(e)

        if t_ko_f & t_expr:
            # If both cell growth and expression is statistically significantly
            # different, we take the gene as candidate
            count += 1
            pbar.set_description(f'Genes found : ' + str(count) + ' ')
            genes_found.append(gene)
            print("\nKO mean for " + gene +
                  " is " + str(ttest_ind(for_tt_lost, for_tt_other,
                               nan_policy='raise', alternative='less')) +
                  " . In other:" + str(for_tt_other.mean()) +
                  " std: " + str(for_tt_other.std()) +
                  ". In lost:" + str(for_tt_lost.mean()) +
                  " std: " + str(for_tt_lost.std()))
            print("\nExpression mean for " + gene +
                  " is " + str(ttest_ind(expr_in_lost, explr_in_other,
                               nan_policy='raise', alternative='greater')) +
                  " . In other:" + str(explr_in_other.mean()) +
                  " std: " + str(explr_in_other.std()) +
                  ". In lost:" + str(expr_in_lost.mean()) +
                  " std: " + str(expr_in_lost.std()))

    print(f'\nPHASE 1 RESULT:')
    for i in range(len(genes_found)):
        print(f'  {i+1}. {genes_found[i]}')
    print()
