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
		score, nTSSs,
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

## Plotting TSR Genomic Distribution
## ----------

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
			text=element_text(size=18)
		)

	png(
		file.path("distribution_plots", paste0("Genomic-Distribution_", x, ".png")),
		width=500, height=150
	)
	print(p)
	dev.off() 
}

dir.create("distribution_plots")

map(names(TSRs), ~plot.distribution(.))
