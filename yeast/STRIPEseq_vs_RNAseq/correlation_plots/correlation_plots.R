#!/usr/bin/env Rscript

library("tidyverse")
library("gtools")

###########################################
## STRIPE-seq vs RNA-seq Correlation Plots
###########################################

## Prepare Data
## ----------

## Load data.

TMM <- read.delim(
	"../TMM_normalize/TMM_normalized_reads/TMM_normalized_reads.tsv",
	sep="\t", header=T, stringsAsFactors=F
) %>% 
	as_tibble %>%
	dplyr::select(-sgd_id)

## Get log2+1 counts.

counts <- mutate_all(TMM, ~log2(.+1))

## Get combinations of comparisons.

comparisons <- counts %>%
	colnames %>%
	combinations(9, 2, .) %>%
	as_tibble %>%
	dplyr::rename("sample_1"=1, "sample_2"=2) %>%
	filter(sample_1 != sample_2)

## Plot Data
## ----------

dir.create("corr_plots")

pmap(
	comparisons,
	function(sample_1, sample_2) {
		p <- ggplot(counts, aes_string(x=sample_1, y=sample_2)) +
			geom_point() +
			theme_bw() +
			theme(text=element_text(size=24)) +
			geom_abline(intercept=0, slope=1, lty=2)
		
		ggsave(
			file.path("corr_plots", paste0(sample_1, "-vs-", sample_2, ".tiff")),
			plot=p, device="tiff", type="cairo", height=6.5, width=6.5
		)
	}
)

## Spearman's Correlation Matrix
## ----------

## Create output directory.

dir.create("spearmans_corr")

## Calculate Spearman's rho.

correlation <- TMM %>%
	as.matrix %>%
	cor(method="spearman") %>%
	as_tibble(rownames="sample_1") %>%
	gather(key="sample_2", value="spearman", -sample_1) %>%
	mutate(spearman=round(spearman, 3))

## Plot correlation matrix.

p <- ggplot(correlation, aes(x=sample_1, y=sample_2, fill=spearman, label=spearman)) +
	geom_tile(color="white", lwd=0.5) +
	geom_label(color="white", label.size=NA, fill=NA) +
	scale_fill_viridis_c(limits=c(0.7,1), name="Spearman") +
	theme_minimal() +
	theme(
		axis.text.x=element_text(angle=45, hjust=1),
		panel.grid=element_blank(),
		axis.title=element_blank()
	)

dir.create("corr_matrix_plot")

ggsave(
	file.path("corr_matrix_plot", "corr_matrix_plot.pdf"),
	plot=p, device=cairo_pdf, height=5.5, width=7
)

## Export text results.

write.table(
	correlation, file.path("spearmans_corr", "spearmans_correlation_matrix.tsv"),
	sep="\t", col.names=T, row.names=F, quote=F
)
