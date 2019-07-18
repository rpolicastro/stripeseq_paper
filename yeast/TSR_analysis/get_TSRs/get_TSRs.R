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

TSSs <- readRDS("../../TSS_analysis/get_TSSs/TSRchitect_TSS_Object.RDS")

## Load determined thresholds.

source("thresholds_values.conf")

threshold.values <- threshold.values %>%
	bind_rows %>%
	gather(key="sample", value="threshold")

## Find TSRs
## ----------

TSRs <- TSSs

## Function to Find TSRs.

for (samp in TSSs@sampleNames) {
	threshold <- threshold.values %>%
		filter(sample == samp) %>%
		pull(threshold)

	sample.index <- which(TSSs@sampleNames == samp) %>%
		as.character
	
	TSRs <- determineTSR(
		experimentName=TSRs,
		n.cores=4,
		tssSetType="replicates",
		tssSet=sample.index,
		tagCountThreshold=threshold,
		clustDist=25,
		writeTable=FALSE,
		mixedorder=FALSE
	)
}

## Convert TSRs to GRanges
## ----------

## Convert to GRanges.

TSRs.GRanges <- TSRs@tsrData %>%
	map(
		~dplyr::rename(., "seqnames"=seq) %>%
		makeGRangesFromDataFrame(keep.extra.columns=TRUE)
	) %>%
	setNames(TSRs@sampleNames)

## Save to RDS

saveRDS(TSRs.GRanges, "TSRs.RDS")

## _______________________________________ deprecated code _______________________________________________

## Filter TSSs by Threshold
## ----------

## Filter TSSs.

#TSSs.filtered <- map(
#	names(TSSs),
#	~TSSs[[.]][score(TSSs[[.]]) >= threshold.values[[.]]]
#) %>% setNames(names(TSSs))

## Export filtered TSSs.

#saveRDS(TSSs.filtered, "Filtered_TSSs.RDS")

## Merge TSSs
## ----------

## Function to merge TSSs and get stats.

#merge.TSSs <- function(x) {
#	merged.TSSs <- reduce(
#		x,
#		ignore.strand=FALSE,
#		min.gapwidth=26L,
#		with.revmap=TRUE
#	)
#	mcols(merged.TSSs) <- aggregate(
#		x,
#		merged.TSSs$revmap,
#		score=sum(score),
#		nTSSs=lengths(seqnames)
#	)
#	merged.TSSs$grouping <- NULL
#	return(merged.TSSs)
#}

## Merge TSSs within 25 bases.

#TSRs <- map(TSSs.filtered, ~merge.TSSs(.))

## Export TSRs
## ----------

#saveRDS(TSRs, "TSRs.RDS")
