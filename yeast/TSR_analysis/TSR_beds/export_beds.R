#!/usr/env/bin Rscript

library("tidyverse")
library("GenomicRanges")
library("rtracklayer")

###################
## Export TSR Beds
###################

## Load TSRs.

TSRs <- readRDS("../get_TSRs/TSRs.RDS")

## Export TSRs to bed.

dir.create("bed_files")

map(
	names(TSRs),
	~export(TSRs[[.]], file.path("bed_files", paste0("TSRs_", ., ".bed")), "bed")
)
