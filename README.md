## A consensus variant-to-function score to functionally prioritize variants for disease

This repository contains a snakemake workflow to calculate cV2F scores ad cV2F metrics method.



## Table of Contents

- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)
  - [cv2f score](#cv2fscore)
  - [cv2f metric](#cv2fmetric)

## Installation

Ensure you have Conda installed. If not, you can install it from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Configuration

Edit the config.yaml file to specify your input and output settings:

## Usage

To run the Snakemake workflow, use the provided wrapper script. This script ensures that the correct Conda environment is activated and Snakemake is run with the specified profile.

`chmod +x run_snakemake.sh`

`nohup ./run_snakemake.sh > snakemake.log 2>&1 &`

## Output

### cv2fscore

For the cv2f score, the workflow generates one output file per chromosome with the suffix .cv2f.txt. Each output file will contain the following columns:

1-CHR: Chromosome number

2-BP: Base pair position

3-SNP: SNP identifier

4-CM: 

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
