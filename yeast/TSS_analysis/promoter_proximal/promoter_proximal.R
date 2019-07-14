#!/usr/bin/env Rscript

library("GenomicRagnes")
library("GenomicFeatures")
library("ChIPseeker")
library("tidyverse")

####################################
## Get Promoter Proximal Percentage
####################################

## Annotate TSSs
## ----------

## Load TSS data.

TSSs <- readRDS("../get_TSSs/Yeast_TSSs.RDS")

## Load genomic annotation.

annotation <- makeTxDbFromGFF(
	"../../genome/Saccharomyces_cerevisiae.R64-1-1.97.gtf",
	format="gtf"
)

## Annotate TSSs.

TSSs.annotated <- map(
	TSSs,
	~ annotatePeak(.,
		tssRegion=c(-1000,100),
		TxDb=annotation,
		sameStrand=TRUE
	) %>%
	as_tibble
)

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

## Export results.

dir.create("promoter_data")

promoter.data %>%
	group_split(sample) %>%
	map(
		~ write.table(
			., file.path("promoter_data", paste0("Promoter-Data_", (pull(., sample) %>% unique), ".tsv")),
			sep="\t", col.names=T, row.names=F, quote=F
		)
	)
