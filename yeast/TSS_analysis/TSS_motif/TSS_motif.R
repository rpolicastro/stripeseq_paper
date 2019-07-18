#!/usr/bin/env Rscript

library("GenomicRanges")
library("GenomicFeatures")
library("tidyverse")
library("Rsamtools")

#####################
## TSS Motif Finding
#####################

## Load Necessary Data
## ----------

## Load TSSs.

TSSs <- readRDS("../get_TSSs/Yeast_TSSs.RDS")

## Load genome annotation GTF as TxDb.

GTF <- makeTxDbFromGFF("../../genome/Saccharomyces_cerevisiae.R64-1-1.97.gtf", "gtf")

## Load genome assembly FASTA as Biostrings.

FASTA <- FaFile("../../genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa")

## Get TSS Sequences
## ----------

## Filter out TSSs with less than 3 reads.

TSSs.filtered <- map(TSSs, ~.[score(.) >= 3])

## Add 5 bases to each side of TSS.

TSS.region <- map(
	TSSs.filtered,
	~resize(., width=5, fix="start") %>%
	resize(width=width(.)+5, fix="end")
)

## Get corresponding sequences.

TSS.sequences <- map(TSS.region, ~getSeq(FASTA, .))

## Export Sequence Fasta Files
## ----------

## Create export directory.

dir.create("TSS_sequences")

## Rename seq names.

rename.seqs <- function(x) {names(x) <- sprintf("TSS_%s",1:length(x)); return(x)}

TSS.seq.export <- map(TSS.sequences, ~rename.seqs(.))

## Export as Fasta.

export.fasta <- function(x) {
	writeXStringSet(
		TSS.seq.export[[x]],
		file.path("TSS_sequences", paste0("TSS-Sequences_", x, ".fasta")),
		format="fasta"
	)
}

map(names(TSS.seq.export), ~export.fasta(.))

## Plot Data
## ----------

## Prepare Data for Plotting

TSS.formatted <- TSS.sequences %>%
	map(
		~as.character(.) %>%
		str_split(pattern="", simplify=TRUE) %>%
		as_tibble %>%
		rowid_to_column(var="TSS") %>%
		gather(key="position", value="base", -TSS) %>%
		mutate(position=factor(position, levels=sprintf("V%s", 1:20)))
	)

## Create directory to export to.

dir.create("TSS_motif_heatmaps")

## Function to plot data.

plot.heatmap <- function(x) {
	p <- ggplot(TSS.formatted[[x]], aes(x=position, y=TSS)) +
		geom_tile(aes(fill=base)) +
		scale_fill_viridis_d() +
		theme_minimal() +
		theme(
			axis.title.y=element_blank(),
			axis.text=element_blank(),
			legend.title=element_blank(),
			text=element_text(size=16),
			panel.grid=element_blank()
		)

	png(
		file.path("TSS_motif_heatmaps", paste0("TSS-Sequence-Heatmap_", x, ".png")),
		height=450, width=450
	)
	print(p)
	dev.off()
}

## Plot TSS sequence heatmaps.

map(names(TSS.formatted), ~plot.heatmap(.))
