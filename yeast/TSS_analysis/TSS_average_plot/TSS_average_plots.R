#!/usr/bin/env Rscript

library("tidyverse")

#####################
## TSS Average Plots
#####################

## Load and Prepare Data
## ----------

## Load data.

TSSs <- readRDS("../annotate_TSSs/TSSs_annotated.RDS")

## Grab info needed for plotting.

TSSs <- TSSs %>%
	map(
		~dplyr::select(., distanceToTSS, score) %>%
		filter(score >= 3)
	)

## Split out individual TSSs from summed TSSs.

TSSs.split <- map(
	TSSs,
	~pmap(., function(distanceToTSS, score) rep(distanceToTSS, score)) %>%
		unlist %>%
		enframe(value="distanceToTSS") %>%
		dplyr::select(-name)
)

## Plot Data
## ----------

dir.create("TSS_average_plots")

plot.TSS.averages <- function(x) {
	p <- ggplot(TSSs.split[[x]], aes(distanceToTSS)) +
		geom_density(fill="#431352", color="#431352") +
		xlim(-2000,2000) +
		labs(
			x="Start Codon",
			y="Density",
			title=x
		) +
		theme_bw()

	ggsave(
		file.path("TSS_average_plots", paste0("TSS-Average-Plot_", x, ".pdf")),
		plot=p, device=cairo_pdf, height=4, width=4
	)
}

map(names(TSSs.split), ~plot.TSS.averages(.))
