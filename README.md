# Benchmarking of Deep-Learning Multi-omics Translation Tools using Biologically Relevant Metrics

Name: Preet Nikhil Shah<br>
Subject: Molekulare Biotechnologie<br>
Department: IPMB, Bioinformatik<br>
Supervisors: Prof. Dr. Carl Herrmann, Dr. Robin Droit and Dr. Daria Doncevic<br>
Contact: preet.shah@stud.uni-heidelberg.de

## Repository structure/notebooks
- `Translation.ipynb`: Processing and plotting of translated RNA and ATAC data for lymphoma and BMMC dataset
- `grcs.ipynb`: Inference of GRNs using pySCENIC from true data and calculation of cellular enrichment with AUCell for true and predicted data; Metric: Gene Regulatory Consistency Score
- `pps.ipynb`: Computation of DEGs for lymphoma dataset and inference of enriched pathways in lymphoma B-cells compared to normal B-cells; Comparison of Pathway and Gene level overlap between true and predicted data; Metric: Pathway Preservation Score
- `dacs.ipynb`: Computation of DAPs for ATAC Lymphoma data and calculation of DAP overlap between true and predicted data, Metric: Differential Accessibility Concordance Score
- `gacs.ipynb`: Inference of Gene activity scores from ATAC data using episcanpy package and comparison of HVG overlap with true RNA data, Metric: Gene Activity Concordance Score
- `metrics.py`: Python file containing all relevant functions for calculation of the metrics  
- `archive`: Folder containing all experimantal/developemental code

## Datasets
- Lymphoma dataset: Download .h5 file conatining both RNA and ATAC data at following [link](https://www.10xgenomics.com/datasets/fresh-frozen-lymph-node-with-b-cell-lymphoma-14-k-sorted-nuclei-1-standard-1-0-0).
Cell-type annotation file is available in the `data` folder.
- BMMC dataset: Download .h5ad file conatining both RNA and ATAC data at following [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122) (cell-type annotation included). 

All preprocessed data (.h5ad files for RNA and ATAC) is available in `data` folder, which is not pushed to GitHub. Any required files can be provided upon request.

## Benchmarked Tools
The implementation of each deep-learning tool with its dependencies in docker environments are available in separate GitHub repositories with relavant code for input data processing, model training and translation output generation. These repos were forked from the original GitHub repos of each tool ([BABEL](https://github.com/wukevin/babel), [Polarbear](https://github.com/Noble-Lab/Polarbear) and [scButterfly](https://github.com/BioX-NKU/scButterfly))

- BABEL (Wu et al., 2021): [GitHub](https://github.com/preet2706/babel)<br> 
Input data format: Single filtered h5 file (Feature-Barcode matrix) with RNA and ATAC data together<br>
Output data format: Translated RNA and ATAC data in separate .h5ad files<br>
- Polarbear (Zhang et al., 2022): [GitHub](https://github.com/preet2706/Polarbear)<br> 
Input data format: Separate .mtx files for RNA and ATAC data, .csv file with all cell barcodes, and .csv file with ATAC peak names<br>
Output data format: Translated RNA and ATAC data in separate .txt files<br>
- scButterfly (Cao et al, 2024): [GitHub](https://github.com/preet2706/scButterfly)<br> 
Input data format: Separate .h5ad files for RNA and ATAC data<br>
Output data format: Translated RNA and ATAC data in separate .h5ad objects<br>

## Packages and Versions
Here are all relavant packages listed for the implementation of each model as well as for the Benchmarking study. Additionally a [`requirements.txt`](requirements.txt) file is provided for a detailed list of all installed packages for the Benchmarking study. In each repository for the tools, separate environment.yaml or requirements.txt files can be found.

- BABEL:<br>
pytorch (version 1.3.1) 
- Polarbear:<br>
tensorflow (version 1.15.0) 
- scButterfly:<br>
pytorch (version 2.4.1)
- Benchmarking:<br>
    General packages:
    - anndata (version: 0.8.0)
    - numpy (version: 1.23.5)
    - pandas (version: 1.5.1)
    - scanpy (version: 1.9.1)
    - matplotlib (version: )
    - seaborn (version: 0.12.1)
    - scikit-learn (version: 1.1.3)

    Benchmarking-related packages:
    - gseapy (version: 1.1.8)
    - pyscenic (version: 0.12.1)
    - episcanpy (version: 0.4.0)
    - loompy (version: 3.0.7)
    - MulticoreTSNE (version: 0.1)
