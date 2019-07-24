#!/usr/bin/env Rscript

library("GenomicRanges")
library("tidyverse")
library("gtools")

######################################
## TSR Correlation Between Replicates
######################################

## Loading Data
## ----------

TSRs <- readRDS("../get_TSRs/TSRs.RDS")

## Merge TSRs
## ----------



## Get Overlapping Regions
## ----------

## Comparison combinations.

sample.combinations <- TSRs %>%
	names %>%
	combinations(6, 2, ., repeats.allowed=FALSE) %>%
	as_tibble %>%
	dplyr::rename("sample_1"=1, "sample_2"=2)

## Jaccard Index
## ----------

## Bin TSRs into deciles.

TSRs.decile <- map(
	TSRs,
	~as_tibble(.) %>%
		mutate(decile=ntile(nTAGs,10)) %>%
		makeGRangesFromDataFrame(keep.extra.columns=TRUE)
)

## Calculate jaccard index per decile.

find.jaccard <- function(x, y) {
	x <- makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
	y <- makeGRangesFromDataFrame(y, keep.extra.columns=TRUE)

	intersects <- GenomicRanges::intersect(x, y, ignore.strand = FALSE)
	intersection <- intersects %>% width %>% sum

	union <- GenomicRanges::union(x, y, ignore.strand = FALSE) %>%
		width %>% sum

	ans <- DataFrame(
		intersection,
		union,
		jaccard = intersection/union,
		n_intersections = length(intersects)
	)

	return(ans)
}
