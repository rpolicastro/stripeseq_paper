#!/usr/bin/env Rscript

library("GenomicRanges")
library("tidyverse")
library("edgeR")
library("gtools")

##########################################
## Get TSS Correlation Between Replicates
##########################################

## Prepare Data
## ----------

## Load GRanges object.

TSSs <- readRDS("../get_TSSs/Yeast_TSSs.RDS")

## Prepare counts sheet.

counts <- TSSs %>%
	# Make the TSS name a concatenation of the chromosome, start, end, and strand.
	map(
		~as_tibble(.) %>%
		mutate(position=paste(seqnames, start, end, strand, sep="_")) %>%
		dplyr::select(position, score)
	) %>%
	# Select the TSS name column and the score.
	# Turn the list of tibbles into one tibble with a column specify what tibble the row came from.
	bind_rows(.id="sample") %>%
	# Turn samples into column and TSS names into rows.
	spread(key=sample, value=score, fill=0) %>%
	# Convert tibble to data frame, and turn the TSS name into rownames.
	as.data.frame %>%
	column_to_rownames("position") %>%
	# Convert data frame to count matrix.
	as.matrix

## Export raw counts.

counts %>%
	# Convert back to tibble for export.
	as_tibble(rownames="TSS_Identifier") %>%
	# Export raw counts.
	write.table(
		., "Raw_TSS_Counts.tsv",
		sep="\t", col.names=T, row.names=F, quote=F
	)

## TMM Normalization of Counts
## ----------

## Only keep TSS positions that have at least 3 samples with a read

# disabled for now
#filtered.counts <- counts[rowSums(counts > 1) >= 3,]

## Create EdgeR object.

edger <- DGEList(
	counts=counts,
	group=c(rep(1,3), rep(2,3))
)

## Normalize the counts.

edger <- calcNormFactors(edger, method="TMM")

## Get TMM normalized counts.

TMM <- counts %>%
	# After running calcNormFactors edgeR::cpm returns TMM
	cpm %>%
	# Convert to tibble.
	as_tibble(rownames="TSS_position")

## Export the TMM normalized counts.

write.table(
	TMM, "TMM_Normalized_Counts.tsv",
	sep="\t", col.names=T, row.names=F, quote=F
)

## Plotting Correlation
## ----------

## Prepare data for plotting.

colnames(TMM) <- TMM %>%
        colnames %>%
        gsub(., pattern="-", replacement="_")

## Get combinations of samples to plot.

sample.combinations <- TMM %>%
	colnames %>%
	discard(. == "TSS_position") %>%
	combinations(6, 2, .) %>%
	as_tibble %>%
	dplyr::rename("sample_1"=1, "sample_2"=2) %>%
	dplyr::filter(sample_1 != sample_2)

## Plot the pairwise correlations.

dir.create("corr_plots")

pmap(
	sample.combinations,
	function(sample_1, sample_2) {
		p <- ggplot(TMM, aes_string(x=sample_1, y=sample_2)) +
			geom_point(size=0.25) +
			theme_bw() +
			scale_fill_viridis_d() +
			geom_abline(intercept=0, slope=1, lty=2)
		
		ggsave(
			file.path("corr_plots", paste0(sample_1, "-vs-", sample_2, ".png")),
			plot=p, device="png", width=5, height=5
		)
	}
)
