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
    return grcs_corr/len(auc_mtx_full.columns), grcs_dict
