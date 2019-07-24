#!/usr/bin/env Rscript

library("tidyverse")
library("Rsamtools")
library("GenomicRanges")

#############################
## TSS Dinucleotide Analysis
#############################

## Loading Data
## ----------

## Loading TSSs.

TSSs <- readRDS("../get_TSSs/Yeast_TSSs.RDS")

## Loading genome assembly.

FASTA <- FaFile("../../genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa")

## Get TSS Sequences
## ----------

## Filter out TSSs with less than 3 reads.

TSSs.filtered <- map(TSSs, ~.[score(.) >= 3])

## Add base before TSS.

TSS.region <- map(TSSs.filtered, ~resize(., width=2, fix="end"))

## Get corresponding sequences.

TSS.sequences <- map(TSS.region, ~getSeq(FASTA, .))

## Find Dinucleotide Frequency
## ----------

## Get frequencies.

dinucleotide.frequencies <- map(
	TSS.sequences,
	~as.character(.) %>%
		enframe(name="chr", value="dinucleotide") %>%
		count(dinucleotide, sort=TRUE, name="occurences") %>%
		mutate("frequency"=occurences/sum(occurences))
)

## Plot Frequencies
## ----------

## Set factor order for dinucleotide.

dinucleotide.frequencies <- map(
	dinucleotide.frequencies,
	~mutate(., "dinucleotide"=factor(dinucleotide, levels=pull(., dinucleotide) %>% unique))
)

## Set export directory.

dir.create("dinucleotide_frequencies")

## Function to plot dinucleotide frequencies.

plot.dinucleotide.frequencies <- function(x) {
	p <- ggplot(dinucleotide.frequencies[[x]], aes(x=fct_rev(dinucleotide), y=frequency)) +
		geom_col(width=0.5, aes(fill=frequency)) +
		theme_classic() +
		scale_fill_viridis_c(name="Frequency") +
		coord_flip() +
		labs(
			x="Dinucleotide",
			y="Frequency"
		)

	ggsave(
		file.path("dinucleotide_frequencies", paste0("Dinucleotide-Frequencies_", x, ".pdf")),
		plot=p, device="pdf", width=5, height=3
	)
}

map(names(dinucleotide.frequencies), ~plot.dinucleotide.frequencies(.))
