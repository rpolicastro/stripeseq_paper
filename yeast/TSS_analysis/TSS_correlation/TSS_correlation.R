#!/usr/bin/env Rscript

library("GenomicRanges")
library("tidyverse")
library("edgeR")

##########################################
## Get TSS Correlation Between Replicates
##########################################

## Prepare Data
## ----------

## Load GRanges object.

TSSs <- readRDS("../get_TSSs/Yeast_TSSs.RDS")

## Prepare counts sheet.

counts <- TSSs %>%
	# Convert GRanges to tibbles.
	map(~as_tibble(.)) %>%
	# Make the TSS name a concatenation of the chromosome, start, end, and strand.
	map(~mutate(., position=paste(seqnames, start, end, strand, sep="_"))) %>%
	# Select the TSS name column and the score.
	map(~dplyr::select(., position, score)) %>%
	# Turn the list of tibbles into one tibble with a column specify what tibble the row came from.
	bind_rows(.id="sample") %>%
	# Turn samples into column and TSS names into rows.
	spread(key=sample, value=score, fill=0) %>%
	# Convert tibble to data frame, and turn the TSS name into rownames.
	as.data.frame %>%
	column_to_rownames("position") %>%
	# Convert data frame to count matrix.
	as.matrix

## Export raw counts.

counts %>%
	# Covnert back to data frame for export.
	as.data.frame %>%
	# Change the rownames back to a column.
	rownames_to_column("TSS_Identifier") %>%
	# Export raw counts.
	write.table(
		., "Raw_TSS_Counts.tsv",
		sep="\t", col.names=T, row.names=F, quote=F
	)

## TMM Normalization of Counts
## ----------

## Create EdgeR object.

edger <- DGEList(
	counts=counts,
	group=c(rep(1,3), rep(2,3))
)

## Normalize the counts.

edger <- calcNormFactors(edger)

## Get TMM normalized counts.

TMM <- counts %>%
	# After running calcNormFactors edgeR::cpm returns TMM
	cpm %>%
	# Convert to data frame.
	as.data.frame %>%
	# Convert the rownames to a column.
	rownames_to_column("TSS_position") %>%
	# Convert to tibble for convenience.
	as_tibble

## Export the TMM normalized counts.

write.table(
	TMM, "TMM_Normalized_Counts.tsv",
	sep="\t", col.names=T, row.names=F, quote=F
)

## Plotting Correlation
## ----------


