#!/usr/bin/env Rscript

library("tidyverse")
library("TSRchitect")
library("GenomicRanges")

#############
## Find TSRs
#############

## Prepare Files for Analysis
## ----------

## Get files.

bam.files <- list.files("../../aligned_reads", pattern=".*\\.bam$")

sample.names <- str_replace(bam.files, "_Aligned.out_cleaned.bam", "")
