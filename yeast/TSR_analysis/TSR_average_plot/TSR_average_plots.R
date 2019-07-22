#!/usr/bin/env Rscript

library("tidyverse")

####################
## TSR Average Plot
####################

## Load and Prepare Data
## -----------

## Load data.

TSRs <- readRDS("../annotate_TSRs/Annotated_TSRs.RDS")

## Select only important columns.

TSRs <- map(TSRs,
	~dplyr::select(., distanceToTSS) %>%
	filter(distanceToTSS >= -2000 & distanceToTSS <= 2000)
)

## Plot Data.
## -----------

## Create output directory.

dir.create("TSR_average_plots")

## Plotting function.

plot.TSR.average <- function(x) {
	p <- ggplot(TSRs[[x]], aes(distanceToTSS)) +
		geom_density(color="#34698c", fill="#34698c") +
		theme_bw() +
		labs(
			x="Start Codon",
			y="Density"
		)

	ggsave(
		file.path("TSR_average_plots", paste0("TSR-Average-Plot_", x, ".pdf")),
		plot=p, device="pdf", height=5, width=5
	)
}

## Plot data.

map(names(TSRs), ~plot.TSR.average(.))
