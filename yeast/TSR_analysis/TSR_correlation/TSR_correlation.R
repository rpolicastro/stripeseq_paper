#!/usr/bin/env Rscript

library("GenomicRanges")
library("tidyverse")
library("gtools")
library("edgeR")

######################################
## TSR Correlation Between Replicates
######################################

## Loading Data
## ----------

TSRs <- readRDS("../get_TSRs/TSRs.RDS")

## Merge TSRs
## ----------

merged.ranges <- TSRs %>%
	as("GRangesList") %>%
	unlist %>%
	GenomicRanges::reduce(ignore.strand=FALSE)

names(merged.ranges) <- sprintf("TSR_%s", 1:length(merged.ranges))

## Get Overlapping Regions
## ----------

return.overlaps <- function(x) {
	overlaps <- findOverlaps(query=merged.ranges, subject=TSRs[[x]]) %>%
                as.tibble %>%
                mutate(
                        nTAGs=TSRs[[x]][subjectHits]$nTAGs,
                        TSR_name=names(merged.ranges[queryHits])
                ) %>%
                dplyr::select(-queryHits, -subjectHits) %>%
                group_by(TSR_name) %>%
                summarize(nTAGs=sum(nTAGs))

	return(overlaps)
}

overlapping <- map(names(TSRs), ~return.overlaps(.))
names(overlapping) <- names(TSRs)

## TMM Normalize Reads
## ----------

## Make count matrix.

count.matrix <- overlapping %>%
	bind_rows(.id="sample") %>%
	spread(key=sample, value=nTAGs, fill=0) %>%
	column_to_rownames("TSR_name") %>%
	as.matrix

## Get TMM normalized reads.

log2.TMM <- count.matrix %>%
	DGEList %>%
	calcNormFactors %>%
	cpm %>%
	as_tibble(rownames="TSR_name") %>%
	mutate_if(is.numeric, ~log2(.+1))

## Correlation Plot
## ----------

## Comparison combinations.

# Clean up sample names.
names(log2.TMM) <- log2.TMM %>%
	names %>%
	str_replace_all("-", "_")

# Get combinations.
sample.combinations <- log2.TMM %>%
	names %>%
	combinations(6, 2, ., repeats.allowed=FALSE) %>%
	as_tibble %>%
	dplyr::rename("sample_1"=1, "sample_2"=2)

## Create output directory.

dir.create("TSR_correlation")

## Function to plot

pmap(
	sample.combinations,
	function(sample_1, sample_2) {
		p <- ggplot(log2.TMM, aes_string(x=sample_1, y=sample_2)) +
			geom_point(size=0.5, color="#34698c") +
			geom_abline(slope=1, intercept=0, lty=2) +
			theme_bw() +
			labs(
				x=paste0("log2(", sample_1, "+1)"),
				y=paste0("log2(", sample_2, "+1)")
			)

		ggsave(
			file.path("TSR_correlation", paste0(sample_1, "_vs_", sample_2, ".tiff")),
			plot=p, device="tiff", width=5, height=5
		) 
	}
)
