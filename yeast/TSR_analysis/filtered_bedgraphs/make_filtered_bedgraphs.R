#!/usr/bin/env Rscript

library("tidyverse")
library("GenomicRanges")
library("rtracklayer")

###############################
## Make Filtered TSS Bedgraphs
###############################

## Load TSSs.

TSSs <- readRDS("../get_TSRs/Filtered_TSSs.RDS")

## Export TSSs.

dir.create("bedgraphs")

export.bedgraphs <- function(x) {
	min.name <- paste0("min_", x, ".bedgraph")
	min <- TSSs[[x]][strand(TSSs[[x]]) == "-"]
	export(min, file.path("bedgraphs", min.name), "bedgraph")

	pos.name <- paste0("pos_", x, ".bedgraph")
	pos <- TSSs[[x]][strand(TSSs[[x]]) == "+"]
	export(pos, file.path("bedgraphs", pos.name), "bedgraph")
}

map(names(TSSs), ~export.bedgraphs(.))
