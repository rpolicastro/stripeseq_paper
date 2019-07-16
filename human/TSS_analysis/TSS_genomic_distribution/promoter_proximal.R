#!/usr/bin/env Rscript

library("GenomicRanges")
library("tidyverse")

####################################
## Get Promoter Proximal Percentage
####################################

## Load Annotated TSSs
## ----------

TSSs.annotated <- readRDS("../annotate_TSSs/TSSs_annotated.RDS")

## Calculate Promoter Proximal
## ----------

promoter.data <- tibble()

for (threshold in 1:50) {
	# Get info on unqiue TSSs.
	perc.unique <- TSSs.annotated %>%
		map(
			~ filter(., score >= threshold) %>%
			count(annotation, name="n.unique") %>%
			mutate("perc.unique"=n.unique/sum(n.unique)) %>%
			filter(annotation=="Promoter") %>%
			mutate("threshold"=threshold)
		) %>%
		bind_rows(.id="sample")

	# Get info on total TSS counts.
	perc.total <- TSSs.annotated %>%
		map(
			~ filter(., score >= threshold) %>%
			group_by(annotation) %>%
			summarize("n.total"=sum(score)) %>%
			mutate("perc.total"=n.total/sum(n.total)) %>%
			filter(annotation=="Promoter") %>%
			mutate("threshold"=threshold) %>%
			dplyr::select(-annotation)
		) %>%
		bind_rows(.id="sample")

	# Get gene info.
	n.genes <- TSSs.annotated %>%
		map(
			~ filter(., score >= threshold) %>%
			pull(geneId) %>%
			unique %>%
			length
		) %>%
		bind_rows %>%
		gather(key="sample", value="n.genes") %>%
		mutate("threshold"=threshold)

	# Combine the unique TSS, total TSS, and gene info.
	threshold.data <- left_join(perc.unique, n.genes, by=c("sample", "threshold"))
	threshold.data <- left_join(threshold.data, perc.total, by=c("sample", "threshold"))
	promoter.data <- bind_rows(promoter.data, threshold.data)
}

## Reorder columns.

promoter.data <- promoter.data %>% dplyr::select(
	sample, annotation, threshold, n.unique,
	perc.unique, n.total, perc.total, n.genes
)

## Export Results
## ----------

## Export to RDS.

promoter.data %>%
	group_split(sample) %>%
	saveRDS(., "Promoter_Data.RDS")

## Export to files.

dir.create("promoter_data")

promoter.data %>%
	group_split(sample) %>%
	map(
		~ write.table(
			., file.path("promoter_data", paste0("Promoter-Data_", (pull(., sample) %>% unique), ".tsv")),
			sep="\t", col.names=T, row.names=F, quote=F
		)
	)

## Plot Promoter Proximal
## ----------

## Split data.

promoter.data.split <- group_split(promoter.data, sample)

## Plot data.

dir.create("promoter_plots")

plot.promoters <- function(x) {
	p <- ggplot(x, aes(x=threshold, y=perc.unique)) +
		geom_line() +
		geom_point(aes(size=n.genes, color=n.genes)) +
		scale_color_viridis_c() +
		ggtitle(x %>% pull(sample) %>% unique) +
		theme(text=element_text(size=24))

	file.name <- paste0("Promoter-Plot_", (x %>% pull(sample) %>% unique), ".png")
	
	png(file.path("promoter_plots", file.name), width=850, height=480)
	print(p)
	dev.off()
}

map(promoter.data.split, ~plot.promoters(.))
