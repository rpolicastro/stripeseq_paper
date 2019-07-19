#!/usr/bin/env Rscript

library("GenomicRanges")
library("GenomicFeatures")
library("tidyverse")
library("Rsamtools")
library("ggseqlogo")

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

## Get TSS rank order to order sequences later.

rank.order <- TSSs.filtered %>%
	map(
		~as_tibble(.) %>%
		mutate(rank=rank(desc(score), ties="first")) %>%
		pull(rank)
	)

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

TSSs.formatted <- map2(
	TSS.sequences, rank.order,
	~as.character(.x) %>%
		.[.y] %>%
		str_split(pattern="", simplify=TRUE) %>%
		as_tibble %>%
		setNames(1:10) %>%
		rowid_to_column(var="sequence") %>%
		gather(key="Position", value="base", -sequence) %>%
		mutate(Position=as.integer(Position))
)

## Create directory to export to.

dir.create("TSS_motif_heatmaps")

## Function to plot data.

plot.heatmap <- function(x) {
	p <- ggplot(TSSs.formatted[[x]], aes(x=Position, y=sequence)) +
		geom_tile(aes(fill=base)) +
		scale_fill_viridis_d() +
		theme_minimal() +
		theme(
			axis.title.y=element_blank(),
			axis.text.y=element_blank(),
			legend.title=element_blank(),
			axis.title.x=element_text(size=16, margin=margin(t=15)),
			panel.grid=element_blank()
		) +
		scale_x_continuous(breaks=c(1,5,10))

	ggsave(
		file.path("TSS_motif_heatmaps", paste0("TSS-Sequence-Heatmap_", x, ".png")),
		plot=p, device="png", width=4, height=4
	)
}

## Plot TSS sequence heatmaps.

map(names(TSSs.formatted), ~plot.heatmap(.))

## Make Sequence Logo
## ----------

## Get consensus matrix for base.

consensus.matrix <- map(
	TSS.sequences,
	~consensusMatrix(.,as.prob=TRUE) %>%
		.[rownames(.) %in% c("A", "C", "G", "T"),]
)

## Make output directory.

dir.create("sequence_logos")

## Create viridis color scheme for bases.

viridis.bases <- make_col_scheme(
	chars=c("A", "C", "G", "T"),
	groups=c("A", "C", "G", "T"),
	cols=c("#431352", "#34698c", "#44b57b", "#fde540")
)

## Output sequence logos.

plot.seqlogos <- function(x) {
	p <- ggplot() +
		geom_logo(consensus.matrix[[x]], col_scheme=viridis.bases) +
		theme_logo()

	ggsave(
		file.path("sequence_logos", paste0("Sequence-Logo_", x, ".png")),
		plot=p, device="png", width=5, height=1.5
	)
}

map(names(consensus.matrix), ~plot.seqlogos(.))
