#!/bin/bash

library("rtracklayer")
library("tidyverse")
library("GenomicRanges")

######################
## Make TSS Bedgraphs
######################

## Prepare Count Data
## ----------

## Load TMM normalized counts from the TSS correlation analysis.

TMM <- read.delim(
	"../TSS_correlation/TMM_Normalized_Counts.tsv",
	sep="\t", header=T, stringsAsFactors=F
) %>% as_tibble(.name_repair="unique")

## Format the data for export.

grange <- TMM %>%
	gather(key=sample, value=score, -TSS_position) %>%
	separate(TSS_position, into=c("seqnames", "start", "end", "strand"), sep="_") %>%
	group_split(sample, strand) %>%
	map(~makeGRangesFromDataFrame(., keep.extra.columns=TRUE))

## Export Bedgraphs
## ----------

dir.create("bedgraphs", showWarnings=FALSE)

export.bedgraphs <- function(x) {
	sample.name <- x %>% 
		.$sample %>%
		unique %>%
		gsub(., pattern="\\.", replacement="-")
	strand <- x %>% 
		strand %>%
		unique %>%
		as.character
	
	if (strand == "+") {
		file.name <- paste0("pos_",sample.name,".bedgraph")
	} else {
		file.name <- paste0("min_",sample.name,".bedgraph")
	}
	
	export(x, file.path("bedgraphs", file.name), "bedgraph")
}

invisible(map(grange, ~export.bedgraphs(.)))
