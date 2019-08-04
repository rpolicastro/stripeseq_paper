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

TSSs <- readRDS("../get_TSSs/Yeast_TSSs.RDS")

## Load genomic annotation.

annotation <- makeTxDbFromGFF(
	"../../genome/Saccharomyces_cerevisiae.R64-1-1.97.gtf",
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
	as_tibble(.name_repair="unique")
)

## Export Annotated TSSs
## ----------

##  Export to RDS.

saveRDS(TSSs.annotated, "TSSs_annotated.RDS")

##  Save to file.

dir.create("annotated_TSSs", showWarnings=FALSE)

map2(
	TSSs.annotated, names(TSSs.annotated),
	~ write.table(
		.x, file.path("annotated_TSSs", paste0("Annotated-TSSs_", .y, ".tsv")),
		sep="\t", col.names=T, row.names=F, quote=F, na=""
	)
)
