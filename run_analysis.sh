#!/bin/bash

conda env create -f environment.yml
conda activate cf_rna-1.0.0

Rscript cluster_analysis.R
