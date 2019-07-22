#!/usr/bin/env Rscript

library("tidyverse")

#####################
## Make TSR Heatmaps
#####################

## Load and Prepare Data
## ----------

## Load data.

TSRs <- readRDS("../annotate_TSRs/Annotated_TSRs.RDS")

## Select useful columns.

TSRs <- map(
	TSRs,
	~dplyr::select(.,
		strand, start, end, geneId,
		geneStart, geneEnd, nTAGs
	)
)

## Prepare data.

TSRs <- TSRs %>% map(
	~mutate(.,
		startDist = case_when(
			strand == "+" ~ start - geneStart,
			strand == "-" ~ -(end - geneEnd)
		),
		endDist = case_when(
			strand == "+" ~ end - geneStart,
			strand == "-" ~ -(start - geneEnd)
		)
	) %>%
	dplyr::select(startDist, endDist, nTAGs, geneId) %>%
	pmap(function(startDist, endDist, nTAGs, geneId) {
		tibble(
			position=seq(startDist, endDist, 1),
			score=nTAGs,
			gene=geneId
		)
	}) %>%
	bind_rows()
)

## Plotting Data

plot.data <- map(
	TSRs,
	~dplyr::filter(., position >= -250 & position <= 250) %>%
	mutate(position=factor(position, levels=-250:250))
)

plot.data %>%
	pluck(1) %>%
	ggplot(aes(x=position, y=gene, fill=score)) +
		geom_tile() +
		theme_minimal() +
		scale_fill_viridis_c() +
		theme(
			axis.text.y=element_blank(),
			axis.text.x=element_blank()
		) +
		scale_fill_continuous(na.value = 0)
