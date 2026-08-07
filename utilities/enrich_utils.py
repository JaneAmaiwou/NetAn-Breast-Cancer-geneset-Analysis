import pandas as pd
import scipy.stats as sts
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

__author__ = 'Rasmus Magnusson, Eduard Ghemes'
__COPYRIGHT__ = 'Copyright (C) 2024 Rasmus Magnusson, Eduard Ghemes'
__contact__ = 'rasma774@gmail.com'

PATH = os.path.dirname(os.path.abspath(__file__)).replace('utilities', '')

# Helper functions
def _load_annotations(annotations):
    if annotations == 'GO':
        return pd.read_csv(os.path.join(PATH, 'data/GO/GO_clean.csv'), index_col=0)
    if annotations == 'GO_BP':
        return pd.read_csv(os.path.join(PATH, 'data/GO/GO_clean_biological_process.csv'), index_col=0)
    if annotations == 'GO_MF':
        return pd.read_csv(os.path.join(PATH, 'data/GO/GO_clean_molecular_function.csv'), index_col=0)
    if annotations == 'GO_CC':
        return pd.read_csv(os.path.join(PATH, 'data/GO/GO_clean_cellular_component.csv'), index_col=0)
    elif annotations == 'KEGG':
        return pd.read_csv(os.path.join(PATH, 'data/KEGG/KEGG_clean.csv'), index_col=0)
    return annotations

def _filter_annotations(annotations, thresh):
    unique_annots, counts = np.unique(annotations.index, return_counts=True)
    return annotations.loc[unique_annots[counts > thresh]]

def _perform_fisher_test(genes, annot_genes, ngenes_background):
    in_both = np.isin(genes, annot_genes).sum()
    only_genes = len(genes) - in_both
    only_annot = len(annot_genes) - in_both
    not_in_any = ngenes_background - (in_both + only_annot + only_genes)
    return sts.fisher_exact([[in_both, only_annot], [only_genes, not_in_any]], alternative='greater')

# Main function
def calc_enrich(
    genes: list, 
    annotations: str = 'GO', 
    ngenes_background: int = 22000, 
    thresh: int = 20
    ) -> pd.DataFrame:
    """
    Calculate the enrichment of a set of genes against a specified annotation set using Fisher's Exact Test.

    Parameters:
    genes: A list of genes for which to calculate enrichment.
    annotations [Default = 'GO']: A string specifying a built-in annotation ('GO' or 'KEGG').
    ngenes_background [Default = 22000]: The total number of genes in the background set.
    thresh [Default = 20]: The minimum number of genes an annotation must be associated with to be considered.

    Returns:
    pandas.DataFrame: Contains the enrichment results, including odds ratios and adjusted p-values, sorted by p-value.
    """
    np.random.seed(0)  # Set a seed for reproducibility

    if isinstance(annotations, str):
        annotations = _load_annotations(annotations)

    annotations = _filter_annotations(annotations, thresh)

    all_res = {}
    all_genes = {}
    print('Calculating Fisher test')

    for term in np.unique(annotations.index):
        annot_genes = np.unique(annotations.loc[term].values.T[0])
        user_genes_in_term = np.intersect1d(genes, annot_genes)
        fishres = _perform_fisher_test(genes, annot_genes, ngenes_background)
        
        all_res[term] = {'OR': fishres[0], 'pval': fishres[1]}
        all_genes[term] = ', '.join(user_genes_in_term)

    df = pd.DataFrame(all_res).transpose()
    df['p_adj'] = multipletests(df.pval.values, method='fdr_bh')[1]
    df['genes'] = df.index.map(all_genes)

    return df.sort_values('pval')
