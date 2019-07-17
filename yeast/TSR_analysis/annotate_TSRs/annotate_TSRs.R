#!/usr/bin/env Rscript

library("GenomicRanges")
library("GenomicFeatures")
library("tidyverse")
library("ChIPseeker")

#################
## Annotate TSRs
#################

## Load Necessary Data
## ----------

## Load TSRs.

TSRs <- readRDS("../get_TSRs/TSRs.RDS")

## Load GTF as TxDb.

annotations <- makeTxDbFromGFF("../../genome/Saccharomyces_cerevisiae.R64-1-1.97.gtf", "gtf")

## Annotate TSRs.

TSRs.anno <- map(
	TSRs,
	~annotatePeak(.,
		tssRegion=c(-1000,100),
		TxDb=annotations,
		sameStrand=TRUE
	) %>%
	as_tibble
)

## Export Annotated TSRs.
## ----------

## Export to R object.

saveRDS(TSRs.anno, "Annotated_TSRs.RDS")

## Export to tsv.

dir.create("annotated_TSR_files")

map(
	names(TSRs.anno),
	~write.table(
		TSRs.anno[[.]], file.path("annotated_TSR_files", paste0("Annotated-TSRs_", ., ".tsv")),
		sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
	)
)
