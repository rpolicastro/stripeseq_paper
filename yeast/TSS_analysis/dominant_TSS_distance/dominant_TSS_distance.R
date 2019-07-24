#!/usr/bin/env Rscript

library("tidyverse")

########################################
## Distance of Major TSS to Start Codon
########################################

## Load and Prepare Data
## ----------

## Load annotated TSSs.

TSSs <- readRDS("../annotate_TSSs/TSSs_annotated.RDS")

## Keep only important columns.

TSSs <- map(TSSs, ~dplyr::select(., geneId, distanceToTSS, score))

## Get strongest TSS for each gene.

dominant.TSSs <- map(
	TSSs,
	~filter(., score >= 3) %>%
		group_by(geneId) %>%
		filter(score == max(score)) %>%
		ungroup
)

## Plot Data
## ----------

## Get data ready for plotting.

distance.summary <- map(
	dominant.TSSs,
	~count(., distanceToTSS) %>%
		filter(distanceToTSS >= -2000 & distanceToTSS <= 2000) %>%
		pmap(function(distanceToTSS, n) rep(distanceToTSS, n)) %>%
		unlist %>%
		enframe(value="distanceToTSS") %>%
		dplyr::select(-name)
)

## output directory.

dir.create("dominant_TSS_distance_plots")

## Function to plot dominant TSS distance to start codon.

plot.dominant.distance <- function(x) {
	p <- ggplot(distance.summary[[x]], aes(distanceToTSS)) +
		geom_density(fill="#431352", color="#431352") +
		xlim(-200,200) +
		theme_bw() +
		labs(
			x="Start Codon",
			y="Density"
		) +
		geom_vline(xintercept=0, lty=2)

	ggsave(
		file.path("dominant_TSS_distance_plots", paste0("Dominant-TSS-Distance_", x, ".pdf")),
		plot=p, device="pdf", height=4, width=5
	)
}

## Plot data.

map(names(distance.summary), ~plot.dominant.distance(.))
