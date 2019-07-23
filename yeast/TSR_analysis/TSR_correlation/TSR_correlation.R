#!/usr/bin/env Rscript

library("GenomicRanges")
library("tidyverse")
library("gtools")

######################################
## TSR Correlation Between Replicates
######################################

## Loading Data
## ----------

TSRs <- readRDS("../get_TSRs/TSRs.RDS")

## Get Overlapping Regions
## ----------

## Comparison combinations.

sample.combinations <- TSRs %>%
	names %>%
	combinations(6, 2, ., repeats.allowed=FALSE) %>%
	as_tibble %>%
	dplyr::rename("sample_1"=1, "sample_2"=2)
