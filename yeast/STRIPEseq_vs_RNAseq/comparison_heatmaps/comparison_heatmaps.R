#!/usr/bin/env Rscript

library("tidyverse")

######################################
## RNA-seq versus STRIPE-seq Heatmaps
######################################

## Load and Prepare Data
## ----------

## Load TMM normalized read counts.

tmm <- read.delim(
	"../TMM_normalize/TMM_normalized_reads/TMM_normalized_reads.tsv",
	sep="\t", header=TRUE, stringsAsFactors=FALSE
) %>% as_tibble

## Get sorting order of genes.

sorting.order <- tmm %>%
	mutate(rnaseq_avg=(yeast_rnaseq_1 + yeast_rnaseq_2 + yeast_rnaseq_3)/3) %>%
	arrange(desc(rnaseq_avg)) %>%
	pull(sgd_id)

## Prepare data for plotting.

tmm.cleaned <- tmm %>%
	mutate_if(is.numeric, ~log2(.+1)) %>%
	gather(key="sample", value="expression", -sgd_id) %>%
	mutate(sgd_id=factor(sgd_id, levels=sorting.order))

## Plot Data
## ----------

ggplot(tmm.cleaned, aes(x=sample, y=fct_rev(sgd_id), fill=expression)) +
	geom_tile(width=0.95) +
	theme_minimal() +
	scale_fill_viridis_c() +
	ylab("log2(TMM+1)") +
	theme(
		axis.text.y=element_blank(),
		panel.grid=element_blank(),
		axis.text.x=element_text(angle=45, hjust=1),
		axis.title.x=element_blank()
	)
