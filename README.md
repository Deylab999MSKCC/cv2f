## A consensus variant-to-function score to functionally prioritize variants for disease

This repository contains a snakemake workflow to calculate cV2F scores and cV2F metrics method. These two methods can be run indepenedently. 


## Table of Contents

- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)
  - [cv2f score](#cv2fscore)
  - [cv2f metric](#cv2fmetric)

## Installation

Ensure you have Conda installed. If not, you can install it from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).


Create and activate a Conda environment with the necessary dependencies:

`conda create -n cv2f -c bioconda -c conda-forge snakemake r-base=4.2.0`
`conda activate cv2f`

`conda install conda-forge::r-devtools`
`conda install conda-forge::r-matrix`
`conda install conda-forge::r-ggplot2`
`conda install conda-forge::r-ggforce`
`conda install conda-forge::r-ggExtra`
`conda install conda-forge::r-ggpubr`


ERROR: dependencies ‘AnnotationFilter’, ‘BiocGenerics’, ‘ensembldb’, ‘GenomeInfoDb’, ‘GenomicRanges’, ‘IRanges’, ‘rtracklayer’ are not available for package ‘locuszoomr’

Ensure the following R packages are installed in your environment:

`install.packages(c("data.table", "xgboost", "pROC", "PRROC", "optparse", "R.utils", "lightgbm", "zoo"))`

BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v75"))

Warning: 

## Configuration

Edit the config.yaml file to specify your input and output settings:

1-snpcell: Path to the SNP cell directory.

2-featurecell: Path to the feature cell directory.

3-mafpath: Path to the MAF features file.

4-bimpath: Path to the BIMS directory.

5-ldblockspath: Path to the LAVA LD blocks file.

6-outputcell: Path to the output directory.

7-script_path: Path to the R script for calculating cV2F scores.

8-chromosomes: List of chromosomes to process.

9-pos_prefix: Prefix for the positive set file.

10-neg_prefix: Prefix for the negative set file.

11-annotation_prefix: Prefix for the annotation file.

12-output_prefix: Prefix for the output file.

## Usage

There are two seperate Snakemake workflows, one for cv2f score and the other is cv2f method. To run either, ensure that your Conda environment is activated and then simply execute the following command in the respective method folder (score or metric):

`snakemake --cores 1`

If you want to run the workflow in the background:

`nohup snakemake --cores 1 > snakemake.log 2>&1 &`


## Output

### cv2fscore

For the cv2f score, the workflow generates one output file per chromosome with the suffix .cv2f.txt. Each output file will contain the following columns:

1-CHR: Chromosome number

2-BP: Base pair position

3-SNP: SNP identifier

4-CM: Centimorgan distance

5-cV2F: cV2F score


Example output file header:
`CHR BP SNP CM cV2F`

### cv2fmetric

For the cv2f metric, the workflow generates one output file for all chromosomes indicated in the config.yaml with the suffix .metrics. The output file will contain the following columns:


1-AUROC.mean: Mean Area Under the Receiver Operating Characteristic curve. This metric evaluates the ability of the model to discriminate between positive and negative classes. 

2-AUROC.sd: Standard deviation of the AUROC.

3-AUPRC.mean: Mean Area Under the Precision-Recall Curve. This metric is useful for imbalanced datasets and measures the trade-off between precision and recall. 

4-AUPRC.sd: Standard deviation of the AUPRC.

Example output file header:
`AUROC.mean	AUROC.sd	AUPRC.mean	AUPRC.sd`
