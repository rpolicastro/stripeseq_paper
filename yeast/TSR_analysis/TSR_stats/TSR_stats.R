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
			log2_nTSSs = log2(nTSSs),
			log2_tsrWdth = log2(tsrWdth)
		) %>%
		gather(key="stat", value="stat_value")
	)

## Plotting Data
## ----------


