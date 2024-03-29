# CF_Varsity

[![DOI](https://zenodo.org/badge/513838306.svg)](https://zenodo.org/badge/latestdoi/513838306)

This repository contains code used for analysis of Cystic Fibrosis Varsity project data.

## Contents

The repository contains the following directories:
- `DESeq2_res` differential gene expression results from DESeq2.
- `gene_counts` normalised gene expression matrices containing log2(FPKM) values.
- `sample_tables` tables of sample information.

## Instructions

You need to first install [`miniconda`](https://docs.conda.io/en/latest/miniconda.html).

To run the code, use the shell script (`run_analysis.sh`), which will create the necessary conda environment and then run the R analysis code.
