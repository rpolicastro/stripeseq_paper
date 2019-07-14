#!/usr/bin/env Rscript

library("TSRchitect")
library("GenomicRanges")
library("rtracklayer")

############
## Get TSSs
############

## TSRchitect to Extract TSSs
## ----------

## Load TSS object.

tss.obj <- loadTSSobj(
	experimentTitle="Yeast STRIPE-seq",
	inputDir="../../aligned_reads/",
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
