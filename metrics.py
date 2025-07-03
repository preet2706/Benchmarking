import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import zlib
import base64
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp

def run_pyscenic(path_adata, data_dir, model:str):
    os.chdir(data_dir)
    f_loom_path_scenic = f"{model}data_scenic.loom"
    f_tfs = 'allTFs_hg38.txt'
    f_motif_path = 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
    f_db_names = 'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'

    rna_data = sc.read_h5ad(path_adata)

    # Create loom file for RNA data
    row_attrs = {
    "Gene": np.array(rna_data.var_names) ,
    }
    col_attrs = {
        "CellID": np.array(rna_data.obs_names) ,
        "nGene": np.array( np.sum(rna_data.X.transpose()>0 , axis=0)).flatten() ,
        "nUMI": np.array( np.sum(rna_data.X.transpose() , axis=0)).flatten() ,
        "cell_type": np.array(rna_data.obs['Cell Types']) ,
    }
    lp.create(f_loom_path_scenic, rna_data.X.transpose(), row_attrs, col_attrs)

    # Run pyscenic
    cmd = f'pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 20'
    os.system(cmd)
    cmd = f'pyscenic ctx adj.csv {f_db_names} --annotations_fname {f_motif_path} --expression_mtx_fname {f_loom_path_scenic} --output reg.csv --mask_dropouts --num_workers 20'
    os.system(cmd)

def compute_auc(input_path, data_dir, loom_file_name, output_file_name):
    os.chdir(data_dir)
    _, extension = os.path.splitext(input_path)

    #Convert h5ad to loom, if input is adata object
    if extension.lower() != '.loom':
        adata = sc.read_h5ad(input_path)
        row_attrs = {
            "Gene": np.array(adata.var_names) ,
        }
        col_attrs = {
            "CellID": np.array(adata.obs_names) ,
            "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
            "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
            "cell_type": np.array(adata.obs['Cell Types']) ,
        }
        lp.create(loom_file_name, adata.X.transpose(), row_attrs, col_attrs)
        f_loom_path_scenic = loom_file_name
    else:
        f_loom_path_scenic = input_path

    # Run Aucell
    cmd = f'pyscenic aucell {f_loom_path_scenic} reg.csv --output {output_file_name} --num_workers 20'
    os.system(cmd)

    # Convert output loom file into pd.DataFrame
    lf = lp.connect(output_file_name, mode='r+', validate=False )
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    lf.close()

    return auc_mtx

def compute_grcs(auc_mtx_full, auc_mtx_pred):
    auc_mtx_true = auc_mtx_full.loc[list(auc_mtx_pred.index),:]
    grcs_corr = 0
    grcs_dict = {}
    for i in auc_mtx_true.columns:
        if i in auc_mtx_pred.columns:
            if auc_mtx_true[i].std() > 0 and auc_mtx_pred[i].std() > 0:
                grcs_corr += np.corrcoef(auc_mtx_true[i], auc_mtx_pred[i])[0,1]
                grcs_dict[i] = np.corrcoef(auc_mtx_true[i], auc_mtx_pred[i])[0,1]
            else:
                grcs_dict[i] = 0
    return grcs_corr/len(auc_mtx_full.columns), grcs_dict

def bootstrap_auc(input_file:str, data_dir:str, n_iterations:int = 50, seed:int = 999):
    grcs = np.zeros(n_iterations)
    for i in range(n_iterations):
        auc_mtx = compute_auc(input_path=input_file, data_dir=data_dir, output_file_name='bootstrap_auc.loom', loom_file_name='')
        auc_mtx2 = compute_auc(input_path=input_file, data_dir=data_dir, output_file_name='bootstrap_auc.loom', loom_file_name='')
        grcs[i] = compute_grcs(auc_mtx, auc_mtx2)[0]
    return {'mean': np.mean(grcs), 'std': np.std(grcs)}   

def compute_deg(RNA_data, hv=False, plot=False):
    RNA_data.obs['group'] = RNA_data.obs['Cell Types'].map(
    lambda x: 'B' if x in ['B'] else ('LYM' if x in ['Tumor B', 'Tumor B cycling'] else None)
        )
    data_gsea = RNA_data[RNA_data.obs['group'].isin(['B','LYM'])].copy()

    if hv is True:
        data_gsea = data_gsea[:, data_gsea.var['highly_variable']]

    sc.pp.log1p(data_gsea)
    sc.tl.rank_genes_groups(
    data_gsea,
    groupby='group',
    groups=['LYM'],
    reference='B',
    method='wilcoxon',
    pts=True
    )

    deg = sc.get.rank_genes_groups_df(data_gsea, group='LYM')
    deg['neg_log10_padj'] = -np.log10(deg['pvals_adj'].replace(0, np.nan))

    deg_list = deg[(deg["pvals_adj"] < 0.05) & (abs(deg["logfoldchanges"]) > 1)]['names'].tolist()

    enr = gp.enrichr(
    gene_list=deg_list,
    gene_sets='MSigDB_Hallmark_2020',  # or 'GO_Biological_Process_2021', 'Reactome_2022', etc.
    organism='Human',
    background=list(data_gsea.var_names),          
    outdir=None,     # optional
    cutoff=0.01                   # p-value cutoff
    )

    sig_enrichments = enr.results[enr.results['Adjusted P-value'] < 0.05]

    if plot is True:  
        df = deg.copy()  
        df['color'] = 'grey'
        df.loc[(df['logfoldchanges'] > 1) & (df['pvals_adj'] < 0.001), 'color'] = 'red'
        df.loc[(df['logfoldchanges'] < -1) & (df['pvals_adj'] < 0.001), 'color'] = 'blue'

        plt.figure(figsize=(7, 6))
        plt.scatter(df['logfoldchanges'], df['neg_log10_padj'], s=10, alpha=0.7, c=df['color'], edgecolor='none')
        plt.axhline(-np.log10(0.001), color='green', linestyle='--')
        plt.axvline(1, color='green', linestyle='dotted')
        plt.axvline(-1, color='green', linestyle='dotted')

        sig_df = df[df['pvals_adj'] < 0.001]
        top_up = sig_df.sort_values(by='logfoldchanges', ascending=False).head(10)
        top_down = sig_df.sort_values(by='logfoldchanges', ascending=True).head(10)

        for _, row in pd.concat([top_up, top_down]).iterrows():
            plt.text(row['logfoldchanges'], row['neg_log10_padj'],
                    row['names'], fontsize=7, ha='center', va='bottom', color='black')

        plt.xlabel('Log2 FC')
        plt.ylabel('-Log10 Adj_p')
        plt.title('Volcano Plot: Lymphoma B vs B Cells')
        plt.tight_layout()
        plt.savefig('my_plot.svg', format='svg')

    return sig_enrichments, deg

def compute_pathway_overlap(sig_true, sig_pred):
# Identify common enriched terms
    shared_terms = set(sig_true['Term']).intersection(sig_pred['Term'])

    dict = {}
    jaccards = 0.0
    for term in shared_terms:
        # Get gene sets for the term from both dataframes
        genes_true = set(sig_true[sig_true['Term'] == term]['Genes'].values[0].split(';'))
        genes_pred = set(sig_pred[sig_pred['Term'] == term]['Genes'].values[0].split(';'))

        # Compute Jaccard index (overlap)
        union = genes_true.union(genes_pred)
        if len(union) == 0:
            continue  # avoid division by zero
        jaccard = len(genes_true.intersection(genes_pred)) / len(union)
        jaccards += jaccard
        fp = len(genes_pred - genes_true)
        fn = len(genes_true - genes_pred)
        dict[term] = {
            'Jaccard': jaccard,
            'FP': fp,
            'FN': fn,
        }

    if jaccards > 0:
        return jaccards / len(dict.values()), dict 
    else:
        return 0.0, dict  # or np.nan if you prefer
