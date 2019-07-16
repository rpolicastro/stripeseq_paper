#!/bin/bash

########################
## Rerun Yeast Analysis
########################

## Retrieve TSSs
## ----------

cd ./get_TSSs
Rscript get_TSSs.R
cd ..

## TSS Correlation
## ----------

cd ./TSS_correlation
Rscript TSS_correlation.R
cd ..

## Export TMM Normalized Bedgraphs
## ----------

cd ./make_bedgraphs
Rscript make_bedgraphs.R
cd ..

## Annotate TSSs
## ----------

cd ./annotate_TSSs
Rscript annotate_TSSs.R
cd ..

## TSS Promoter Proximity
## ----------

cd ./TSS_genomic_distribution
Rscript promoter_proximal.R
cd ..

## TSS Average Plots
## ----------

cd ./TSS_average_plots
Rscript TSS_average_plots.R
cd ..
