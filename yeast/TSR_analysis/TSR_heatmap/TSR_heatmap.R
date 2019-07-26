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

## Export TSRs to make average plot later.

saveRDS(TSRs, "TSR_Density.RDS")

## Plotting Data
## ----------

## Prepare data for plotting.

plot.data <- map(
	TSRs,
	~dplyr::filter(., position >= -250 & position <= 250) %>%
		mutate(position=factor(position, levels=-250:250)) %>%
		complete(gene, position, fill=list(score=0)) %>%
		mutate(
			score=log2(score+1),
			position=as.integer(position)
		)
)

## Sort genes by total TSR signal.

sorting.order <- map(
	plot.data,
	~filter(., score > 0) %>%
		count(gene, score) %>%
		group_by(gene) %>%
		summarize(total.score=sum(n)) %>%
		arrange(desc(total.score)) %>%
		rowid_to_column %>%
		dplyr::select(-total.score)
)

plot.data <- map2(
	plot.data, sorting.order,
	~mutate(.x, gene = factor(gene, levels=(arrange(.y, desc(rowid)) %>% pull(gene))))
)

## Put uper limit on score for scale.

plot.data <- map(
	plot.data,
	~mutate(., score=case_when(
		score > 8 ~ 8,
		score <= 8 ~ score
	))
)

## Create output directory.

dir.create("TSR_heatmaps")

## Plot heatmaps.

plot.TSR.heatmaps <- function(x) {
	p <- ggplot(plot.data[[x]], aes(x=position, y=gene, fill=score)) +
		geom_tile() +
		theme_minimal() +
		scale_fill_viridis_c(
			breaks=c(0,2,4,6,8),
			labels=c(0,2,4,6,">8"),
			limits=c(0,8)
		) +
		theme(
			axis.text.y=element_blank(),
			panel.grid=element_blank()
		) +
		geom_vline(xintercept=250, color="white", lty=2) +
		labs(
			x="Position",
			y="Gene"
		) +
		scale_x_continuous(
			breaks=c(0, 150, 250, 350, 500),
			labels=c(-250, -150, 0, 150, 250)
		)

		ggsave(
			file.path("TSR_heatmaps", paste0("TSR-Heatmap_", x, ".tiff")),
			plot=p, device="tiff", type="cairo", height=4, width=3
		)
}

map(names(plot.data), ~plot.TSR.heatmaps(.))
