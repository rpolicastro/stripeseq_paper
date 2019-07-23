#!/usr/bin/env

library("tidyverse")

############################
## TSR Genomic Distribution
############################

## Load and Prepare Data
## ----------

## Loading TSR data.

TSRs <- readRDS("../annotate_TSRs/Annotated_TSRs.RDS")

## Cleaning TSR data.

TSRs <- map(
	TSRs,
	~ dplyr::select(.,
		geneId, transcriptId,
		nTAGs, nTSSs, tsrPeak, tsrWdth, tsrTrq, tsrSI, tsrMSI,
		annotation, distanceToTSS
	) %>%
	mutate(cleaned.annotations=case_when(
		annotation == "Promoter" ~ "promoter",
		grepl(annotation, pattern="Exon") ~ "exon",
		grepl(annotation, pattern="Intron") ~ "intron",
		grepl(annotation, pattern="Downstream") ~ "downstream",
		annotation == "Distal Intergenic" ~ "intergenic"
	))
)

## Export cleaned TSR data as R object.

saveRDS(TSRs, "TSR_data.RDS")

## Plotting TSR Genomic Distribution
## ----------

## Function to plot distribution.

plot.distribution <- function(x) {
	plot.data <- TSRs[[x]] %>%
		mutate(
			"sample"="genomic_distribution",
			"cleaned.annotations"=factor(
				cleaned.annotations,
				levels=c("intergenic", "downstream", "intron", "exon", "promoter")
			)
		)
	
	p <- ggplot(plot.data, aes(sample, fill=cleaned.annotations)) +
		geom_bar(position="fill") +
		scale_fill_viridis_d(direction=-1, name="Annotation") +
		coord_flip() +
		ylab("Fraction") +
		theme_minimal() +
		theme(
			axis.text.y=element_blank(),
			axis.title.y=element_blank(),
			text=element_text(size=18),
			panel.grid=element_blank()
		)

	ggsave(
		file.path("distribution_plots", paste0("Genomic-Distribution_", x, ".pdf")),
		plot=p, device="pdf", width=5, height=2
	)
}

## Create output directory.

dir.create("distribution_plots")

## Plot distributions.

map(names(TSRs), ~plot.distribution(.))

## Exporting Promoter Stats
## ----------

## Create output directory.

dir.create("distribution_stats")

## Get summary stats for TSR distribution.

promoter.stats <- map(
	TSRs,
	~count(., cleaned.annotations, name="nTSRs") %>%
		arrange(desc(nTSRs)) %>%
		mutate("perc"=round(nTSRs/sum(nTSRs), 3))
)

## Export summary stats to tsv.

map(
	names(promoter.stats),
	~write.table(
		promoter.stats[[.]], file.path("distribution_stats", paste0("TSR-Genomic-Distribution_", ., ".tsv")),
		sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
	)
)
