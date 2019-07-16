#!/usr/bin/env Rscript

library("GenomicRanges")
library("GenomicFeatures")
library("ChIPseeker")
library("tidyverse")

####################################
## Get Promoter Proximal Percentage
####################################

## Annotate TSSs
## ----------

## Load TSS data.

TSSs <- readRDS("../get_TSSs/Human_TSSs.RDS")

## Load genomic annotation.

annotation <- makeTxDbFromGFF(
	"../../genome/Homo_sapiens.GRCh38.97.gtf",
	format="gtf"
)

## Annotate TSSs.

TSSs.annotated <- map(
	TSSs,
	~ annotatePeak(.,
		tssRegion=c(-1000,100),
		TxDb=annotation,
		sameStrand=TRUE
	) %>%
	as_tibble
)

## Export Annotated TSSs
## ----------

##  Export to RDS.

saveRDS(TSSs.annotated, "TSSs_annotated.RDS")

##  Save to file.

dir.create("annotated_TSSs")

map(
	names(TSSs.annotated),
	~ write.table(
		TSSs.annotated[[.]], file.path("annotated_TSSs", paste0("Annotated-TSSs_", ., ".tsv")),
		sep="\t", col.names=T, row.names=F, quote=F, na=""
	)
)
