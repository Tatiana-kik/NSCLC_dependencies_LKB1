"""
This file containts helper function for Phase1 analysis.
"""
import pandas as pd
import os


# Phase I: Set of methods to preload DepMap data

def read_cell_lines():
    return pd.read_csv(os.path.join("phase1", "csv_data", "sample_info.csv"))


def read_mutations():
    return pd.read_csv(os.path.join("phase1", "csv_data", "CCLE_mutations.csv"))


def read_expression():
    return pd.read_csv(os.path.join("phase1", "csv_data", "CCLE_expression.csv"))


def read_gene_dependency():
    return pd.read_csv(os.path.join("phase1", "csv_data", "CRISPR_gene_dependency.csv"))


def read_gene_effect():
    return pd.read_csv(os.path.join("phase1", "csv_data", "CRISPR_gene_effect.csv"))


def read_copy_number():
    return pd.read_csv(os.path.join("phase1", "csv_data", "CCLE_gene_cn.csv"))


def read_copy_number_wes():
    return pd.read_csv(os.path.join("phase1", "csv_data", "CCLE_wes_gene_cn.csv"))


def get_cell_lines_for_disease(type, disease):
    cell_lines_raw = read_cell_lines()
    return cell_lines_raw[cell_lines_raw[type] == disease]['DepMap_ID']


def get_mutated_cell_lines(cell_lines,gene):
    mutations_raw = read_mutations()
    mutations_cl = mutations_raw[mutations_raw['DepMap_ID'].isin(cell_lines)]
    mutations_cl_gene = mutations_cl[mutations_cl['Entrez_Gene_Id'] == gene]
    return mutations_cl_gene[mutations_cl_gene['Variant_annotation'] == 'damaging']['DepMap_ID']


def get_cl_with_ko_sens_over_td(gene_str,threshold):
    sensitivity_raw = read_gene_effect()
    cell_lines = sensitivity_raw.loc[:,['DepMap_ID',gene_str]]
    return cell_lines[cell_lines[gene_str] > threshold]['DepMap_ID']


def get_cl_with_lh_for_ko_over_td(gene_str,lh_threshold):
    sensitivity_raw = read_gene_dependency()
    cell_lines = sensitivity_raw.loc[:,['DepMap_ID',gene_str]]
    return cell_lines[cell_lines[gene_str] > lh_threshold]['DepMap_ID']


def get_cl_with_ko_lh(gene_str, th, lhth):
    cl_ge = get_cl_with_ko_sens_over_td(gene_str,th)
    cl_gd = get_cl_with_lh_for_ko_over_td(gene_str,lhth)
    return pd.Series(list(set(cl_ge) & set(cl_gd)))


def get_cl_expr_under_td(gene_str, threshold):
    expr = read_expression()
    cell_lines_onegene_expr = expr.loc[:,['Unnamed: 0',gene_str]]
    return cell_lines_onegene_expr[cell_lines_onegene_expr[gene_str] < threshold]['Unnamed: 0']


def get_cl_with_cn_under_td(gene_str, threshold):
    cn = read_copy_number()
    cell_lines_cn = cn.loc[:, ['Unnamed: 0', gene_str]]
    return cell_lines_cn[cell_lines_cn[gene_str] < threshold]['Unnamed: 0']
