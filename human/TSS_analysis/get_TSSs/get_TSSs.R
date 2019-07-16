#!/usr/bin/env Rscript

library("TSRchitect")
library("GenomicRanges")
library("rtracklayer")
library("tidyverse")

############
## Get TSSs
############

## TSRchitect to Extract TSSs
## ----------

## Load TSS object.

tss.obj <- loadTSSobj(
	experimentTitle="Human STRIPE-seq",
	inputDir="../aligned_reads",
	n.cores = 4,
	isPairedBAM=TRUE,
	isPairedBED=FALSE,
	sampleSheet="sample_sheet.tsv",
	sampleNames="",
	replicateIDs=""
)

## Get fragment 5` ends.

tss.obj <- inputToTSS(tss.obj)

## Find TSSs.

tss.obj <- processTSS(
	tss.obj,
	n.cores=4,
	tssSet="all",
	writeTable=FALSE
)

## Pull TSSs from TSRchitect Object and Export
## ----------

## Grab TSSs from object.

TSSs <- tss.obj@tssCountData %>%
	setNames(tss.obj@sampleNames)

## Convert TSSs to GRanges object.

TSSs.granges <- map(
	TSSs,
	~ GRanges(
		seqnames=.$seq,
		ranges=IRanges(start=.$TSS, width=1),
		strand=.$strand,
		score=.$nTAGs
	)
)

## Export GRanges objects for further analysis.

saveRDS(TSSs.granges, "Human_TSSs.RDS")
