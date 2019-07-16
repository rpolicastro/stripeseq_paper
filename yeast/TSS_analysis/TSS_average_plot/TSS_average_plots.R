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

TSSs <- map(
	TSSs,
	~ dplyr::select(., distanceToTSS, score) %>%
	group_by(distanceToTSS) %>%
	summarize(sum=sum(score))
)

## Plot Data
## ----------

dir.create("TSS_average_plots")

plot.TSS.averages <- function(x) {
	p <- ggplot(TSSs[[x]], aes(x=distanceToTSS, y=sum)) +
		geom_line(fill="dodgerblue", color="dodgerblue", alpha=0.8) +
		xlim(-2000,2000) +
		theme_bw() +
		theme(text=element_text(size=24))

	png(file.path("TSS_average_plots", paste0("TSS-Average-Plot_", x, ".png")), height=480, width=850)
	print(p)
	dev.off()
}

map(names(TSSs), ~plot.TSS.averages(.))
