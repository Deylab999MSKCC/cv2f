## A consensus variant-to-function score to functionally prioritize variants for disease

This repository contains a snakemake workflow to calculate cV2F scores. The workflow is configured to run with LSF and is designed to be flexible and easy to modify. There is also the cV2F metrics method that can be run indepenedently from the cV2F scores: https://github.com/Deylab999MSKCC/snakemake_cv2f_metric



## Table of Contents

- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)

## Installation

Ensure you have Conda installed. If not, you can install it from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Configuration

Edit the config.yaml file to specify your input and output settings:

## Usage

To run the Snakemake workflow, use the provided wrapper script. This script ensures that the correct Conda environment is activated and Snakemake is run with the specified profile.

`chmod +x run_snakemake.sh`

`nohup ./run_snakemake.sh > snakemake.log 2>&1 &`

## Output

The workflow generates one output file per chromosome with the suffix .cv2f.txt. Each output file will contain the following columns:

1-CHR: Chromosome number

2-BP: Base pair position

3-SNP: SNP identifier

4-CM: 

5-cV2F: cV2F score


Example output file header:
`CHR BP SNP CM cV2F`
