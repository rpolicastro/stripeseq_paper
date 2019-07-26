#!/bin/bash

library("tidyverse")

######################
## Get Max UTR Length
######################

## Prepare Data
## ----------

## Load annotated TSSs.

TSSs <- readRDS("../annotate_TSSs/TSSs_annotated.RDS")

## Keep only important colummns.

TSSs <- map(TSSs, ~dplyr::select(., geneId, distanceToTSS))

## Get TSS with minimum distance to start codon.

TSSs.min <- map(
	TSSs,
	~filter(., between(distanceToTSS, -1000, 1000)) %>%
		group_by(geneId) %>%
		summarize(minDistance=min(distanceToTSS))
)

## Plotting Results
## ----------

## Create output directory.

dir.create("max_UTR_length_plots")

## Function to plot results.

plot.max.utr <- function(x) {
	p <- ggplot(TSSs.min[[x]], aes(x=minDistance)) +
		geom_density(fill="#431352", color="#431352") +
		xlim(-1000,1000) +
		theme_bw() +
		labs(
			x="TSS With Minimal Distance to Start Codon",
			y="Density"
		) +
		geom_vline(xintercept=0, lty=2)

	ggsave(
		file.path("max_UTR_length_plots", paste0("Max-UTR-Length_", x, ".pdf")),
		plot=p, device=cairo_pdf, width=5, height=4
	)
}

## Plot results.

map(names(TSSs.min), ~plot.max.utr(.))
