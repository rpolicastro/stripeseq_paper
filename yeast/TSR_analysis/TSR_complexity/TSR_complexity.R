#!/usr/bin/env Rscript

library("tidyverse")

####################################################
## Measures of TSR Complexity (Unique TSSs per TSR)
####################################################

## Load and Prepare Data
## ----------

## Load TSR data.

TSRs <- readRDS("../TSR_distribution/TSR_data.RDS")

## Fix factor levels on data.

TSRs <- map(
	TSRs,
	~mutate(., cleaned.annotations=factor(
		cleaned.annotations,
		levels=c("promoter", "exon", "intron", "downstream", "intergenic")
	))
)

## Plot TSR Complexity
## ----------

## Make outpit directory.

dir.create("complexity_plots")

## Function to make TSR complexity histogram.

plot.complexity <- function(x) {
	p <- ggplot(TSRs[[x]], aes(x=cleaned.annotations, y=log2(nTSSs+1))) +
		geom_jitter(color="grey", size=0.25, width=0.25) +
		geom_boxplot(fill=NA, aes(color=cleaned.annotations), outlier.shape=NA, lwd=1, width=0.25) +
		scale_color_viridis_d() +
		theme_bw() +
		ylim(0,NA) +
		theme(
			text=element_text(size=18),
			legend.position="none"
		) +
		labs(
			y="Log2(Unique TSSs/TSR + 1)",
			x="Annotation"
		)

	ggsave(
		file.path("complexity_plots", paste0("TSR-Complexity-Plot_", x, ".tiff")),
		plot=p, device="tiff", type="cairo", width=6.5, height=4
	)
}

## Plot TSR complexity.

map(names(TSRs), ~plot.complexity(.))
