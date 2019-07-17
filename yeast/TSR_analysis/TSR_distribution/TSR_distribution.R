#!/usr/bin/env

library("tidyverse")

############################
## TSR Genomic Distribution
############################

## Load and Prepare Data
## ----------

## Loading TSR data.

TSRs <- readRDS("../annotate_TSRs/Annotated_TSRs.RDS")

## Cleaning TSR data.

TSRs <- map(
	TSRs,
	~ dplyr::select(.,
		geneId, transcriptId,
		score, nTSSs,
		annotation, distanceToTSS
	) %>%
	mutate(cleaned.annotations=case_when(
		annotation == "Promoter" ~ "promoter",
		grepl(annotation, pattern="Exon") ~ "exon",
		grepl(annotation, pattern="Intron") ~ "intron",
		grepl(annotation, pattern="Downstream") ~ "downstream",
		annotation == "Distal Intergenic" ~ "intergenic"
	))
)
