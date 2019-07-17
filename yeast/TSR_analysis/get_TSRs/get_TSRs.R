#!/usr/bin/env Rscript

library("tidyverse")
library("TSRchitect")
library("GenomicRanges")

#############
## Find TSRs
#############

## Prepare Data
## ----------

## Load TSS GRanges.

TSSs <- readRDS("../../TSS_analysis/get_TSSs/Yeast_TSSs.RDS")

## Load determined thresholds.

source("threshold_values.conf")

## Filter TSSs by Threshold
## ----------

## Filter TSSs.

TSSs.filtered <- map(
	names(TSSs),
	~TSSs[[.]][score(TSSs[[.]]) >= threshold.values[[.]]]
) %>% setNames(names(TSSs))

## Export filtered TSSs.

saveRDS(TSSs.filtered, "Filtered_TSSs.RDS")

## Merge TSSs
## ----------

## Function to merge TSSs and get stats.

merge.TSSs <- function(x) {
	merged.TSSs <- reduce(
		x,
		ignore.strand=FALSE,
		min.gapwidth=26L,
		with.revmap=TRUE
	)
	mcols(merged.TSSs) <- aggregate(
		x,
		merged.TSSs$revmap,
		score=sum(score),
		nTSSs=lengths(seqnames)
	)
	merged.TSSs$grouping <- NULL
	return(merged.TSSs)
}

## Merge TSSs within 25 bases.

TSRs <- map(TSSs.filtered, ~merge.TSSs(.))

## Export TSRs
## ----------

saveRDS(TSRs, "TSRs.RDS")
