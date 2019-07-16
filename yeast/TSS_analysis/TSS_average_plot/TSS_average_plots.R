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

TSSs <- map(TSSs, ~ dplyr::select(., distanceToTSS, score))

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
		geom_density(fill="dodgerblue", color="dodgerblue") +
		xlim(-2000,2000) +
		theme_bw() +
		theme(text=element_text(size=24))

	png(file.path("TSS_average_plots", paste0("TSS-Average-Plot_", x, ".png")), height=480, width=850)
	print(p)
	dev.off()
}

map(names(TSSs.split), ~plot.TSS.averages(.))
