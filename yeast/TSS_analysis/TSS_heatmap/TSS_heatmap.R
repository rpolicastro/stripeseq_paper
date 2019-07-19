#!/usr/bin/env Rscript

library("tidyverse")

###############
## TSS Heatmap
###############

## Load TSS Data
## ----------

## Load data.

TSSs <- readRDS("../annotate_TSSs/TSSs_annotated.RDS")

## Only keep necessary columns.

TSSs <- map(
	TSSs, 
	~dplyr::select(., geneId, distanceToTSS, score) %>%
		filter(
			(distanceToTSS >= -1000 & distanceToTSS <= 1000),
			score >= 3
		) %>%
		group_by(geneId, distanceToTSS) %>%
		summarize(score=sum(score)) %>%
		ungroup %>%
		arrange(geneId, distanceToTSS)
)

## Process Data
## ----------

## Generate count matrix.

count.matrix <- TSSs %>%
	map(., 
		~spread(., key="distanceToTSS", value="score") %>%
		bind_rows %>%
		replace(., is.na(.), 0)
	)

## Arrange by row average.

row.order <- count.matrix %>%
	map(.,
		~gather(., key="position", value="score", -geneId) %>%
		group_by(geneId) %>%
		summarize(average=mean(score)) %>%
		arrange(desc(average)) %>%
		pull(geneId)
	)

count.matrix <- map2(
	count.matrix, row.order,
	~.x[match(.y, .x$geneId),]
)

## Export Count Matrix
## ----------

## Log2 transform.

export.matrix <- map(count.matrix, ~ mutate_if(., is.numeric, ~log2(.+1)))

## To RDS file.

saveRDS(export.matrix, "Log2_TSS_Matrix.RDS")

## To tsv.

dir.create("Log2_TSS_matrices")

map(
	names(export.matrix),
	~write.table(
		export.matrix[[.]], file.path("Log2_TSS_matrices", paste0("Log2-TSS-Matrix_", ., ".tsv")),
		sep="\t", col.names=T, row.names=F, quote=F
	)
)

## Plotting Data
## -----------

## Formatting data for plotting.

to.plot <- map(
	count.matrix,
	~mutate_if(., is.numeric, ~log2(.+1)) %>%
		gather(key="Position", value="log2_score", -geneId) %>%
		mutate(Position=as.integer(Position)) %>%
		dplyr::rename("Gene"=geneId)
)

## Create output directory.

dir.create("TSS_heatmaps")

## Function for plotting data.

plot.TSS.heatmap <- function(x) {
	p <- to.plot[[x]] %>%
		filter(Position >= -250 & Position <= 250) %>%
		ggplot(aes(x=Position, y=Gene, fill=log2_score)) +
			geom_tile() +
			theme_minimal() +
			theme(
				axis.text.y=element_blank(),
				panel.grid=element_blank()
			) +
			scale_fill_viridis_c(limits=c(0,4)) +
			geom_vline(xintercept=0, lty=2, color="white") +
			labs(fill="log2(Score + 1)")

	ggsave(
		file.path("TSS_heatmaps", paste0("TSS-heatmap_", x, ".png")),
		plot=p, device="png", height=4, width=4
	)
}

## Plot heatmaps.

map(names(to.plot), ~plot.TSS.heatmap(.))
