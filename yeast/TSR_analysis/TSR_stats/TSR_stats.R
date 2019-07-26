#!/usr/bin/env Rscript

library("tidyverse")

####################################
## TSR Intensity and Shape Analysis
####################################

## Load and Prepare Data
## ----------

## Loading data.

TSRs <- readRDS("../TSR_distribution/TSR_data.RDS")

## Preparing data for plotting.

TSRs.cleaned <- TSRs %>%
	map(
		~dplyr::select(., -geneId, -transcriptId, -annotation, -cleaned.annotations) %>%
		mutate(
			log2_nTAGs = log2(nTAGs),
			log2_nTSSs = log2(nTSSs+1),
			log2_tsrWdth = log2(tsrWdth+1)
		) %>%
		dplyr::select(-nTAGs, -nTSSs, -tsrWdth) %>%
		gather(key="stat", value="stat_value") %>%
		mutate("sample"="samp")
	)

## Plotting Data
## ----------

## Rename columns.



## Create output directory.

dir.create("TSR_stat_plots")

## Function to plots TSR stats.

plot.tsr.stats <- function(x) {
	p <- ggplot(TSRs.cleaned[[x]], aes(x=sample, y=stat_value, color=stat)) +
		geom_jitter(size=0.1) +
		facet_wrap(~stat, ncol=4, scales="free") +
		scale_color_viridis_d() +
		geom_boxplot(color="#1c1c1c", width=0.25, fill=NA, outlier.shape=NA) +
		theme_minimal() +
		theme(
			axis.text.x=element_blank(),
			legend.position="none",
			axis.title.x=element_blank(),
			axis.title.y=element_blank()
		)

	ggsave(
		file.path("TSR_stat_plots", paste0("TSR-Stat-Plot_", x, ".tiff")),
		plot=p, device="tiff", type="cairo", height=3.5, width=5
	)
}

## Plot TSR stats.

map(names(TSRs.cleaned), ~plot.tsr.stats(.))
