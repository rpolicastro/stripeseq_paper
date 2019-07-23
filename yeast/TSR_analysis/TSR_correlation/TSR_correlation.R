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

## Jaccard Index
## ----------

intersects <- GenomicRanges::intersect(TSRs[[1]], TSRs[[2]], ignore.strand = FALSE)
intersection <- sum(width(intersects))
union <- sum(width(GenomicRanges::union(TSRs[[1]], TSRs[[2]], ignore.strand = FALSE)))
ans <- DataFrame(
	intersection,
	union,
	jaccard = intersection/union,
	n_intersections = length(intersects)
)
