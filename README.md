## A consensus variant-to-function score to functionally prioritize variants for disease

This repository contains a workflow for the cV2F method.


## Table of Contents

- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)


## Installation

Ensure you have Conda installed. If not, you can install it from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Create and activate a Conda environment with the necessary dependencies:

```
conda create -n cv2f -c bioconda -c conda-forge snakemake r-base=4.2.0
conda activate cv2f
conda install conda-forge::r-devtools r-matrix r-ggplot2 r-ggforce r-ggExtra r-ggpubr
```

Ensure the following R packages are installed in your environment:

```
install.packages(c("data.table", "xgboost", "pROC", "PRROC", "optparse", "R.utils", "lightgbm", "zoo"))
BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v86"))
devtools::install_github(c("liuyanguu/SHAPforxgboost", "myles-lewis/locuszoomr"))
```

## Usage

For detailed instructions on running the cV2F model check the [Wiki](https://github.com/Deylab999MSKCC/cv2f/wiki).

## cV2F Data

cV2F scores both tissue-agnostic and tissue-specific can be found in https://mskcc.box.com/shared/static/hsrogtr3fddtmd53hyy5ph7dlp20eq72.txt.

cV2F fine-mapping results can be found in https://mskcc.box.com/s/wyl206gnn5tqhgnpuakzik98yrdonilw.

## Citation

If you find this pipeline or the cV2F scores useful for your work please cite our paper which is out on bioRxiv:

Fabiha et al 2024. A consensus variant-to-function score to functionally prioritize variants for disease. bioRxiv.

## Acknowledgements
Thank you to all co-authors of the cV2F manuscript and the ENCODE4 Consortium.
