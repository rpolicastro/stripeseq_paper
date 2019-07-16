#!/usr/bin/env Rscript

library("tidyverse")
library("edgeR")

########################
## TMM Normalize Counts
########################

## Load RNA-seq Data
## ----------

## Load counts file.

rna.seq <- read.delim(
	"../count_files/Human_RNAseq_counts.tsv",
	sep="\t", header=T, stringsAsFactors=F
) %>% as_tibble

# Clean up column names.

colnames(rna.seq) <- str_replace_all(
	colnames(rna.seq),
	pattern="(X\\..*aligned\\.|_Aligned.*\\.bam)",
	replacement=""
)

rna.seq <- rna.seq %>%
	dplyr::rename(
		"ensembl_id"=Geneid,
		"human_rnaseq_1"=HRNASEQ_01_K562_1,
		"human_rnaseq_2"=HRNASEQ_02_K562_2,
		"human_rnaseq_3"=HRNASEQ_03_K562_3
	) %>%
	dplyr::select(
		ensembl_id,
		human_rnaseq_1,
		human_rnaseq_2,
		human_rnaseq_3
	)

## Load STRIPE-seq Data
## ----------

## Load counts file.

stripe.seq <- read.delim(
	"../count_files/Human_STRIPEseq_counts.tsv",
	sep="\t", header=T, stringsAsFactors=F
) %>% as_tibble

## Clean up column names.

colnames(stripe.seq) <- str_replace_all(
	colnames(stripe.seq),
	pattern="(X\\..*aligned\\.|_Aligned.*\\.bam)",
	replacement=""
)

## Rename columns.

stripe.seq <- stripe.seq %>%
	dplyr::rename(
		"ensembl_id"=Geneid,
		"human_stripeseq_pooled_2"=K562.pooled_WT.100ng_2,
		"human_stripeseq_pooled_1"=K562.pooled_WT.100ng_1,
		"human_stripeseq_pooled_3"=K562.pooled_WT.100ng_3,
		"human_stripeseq_unpooled_3"=K562.unpooled_WT.100ng_3,
		"human_stripeseq_unpooled_2"=K562.unpooled_WT.100ng_2,
		"human_stripeseq_unpooled_1"=K562.unpooled_WT.100ng_1
	) %>%
	dplyr::select(
		ensembl_id,
		human_stripeseq_unpooled_1,
		human_stripeseq_unpooled_2,
		human_stripeseq_unpooled_3,
		human_stripeseq_pooled_1,
		human_stripeseq_pooled_2,
		human_stripeseq_pooled_3
	)

## Merge and Prepare Data
## ----------

## Merge data.

merged <- left_join(rna.seq, stripe.seq, by="ensembl_id")

## Export raw reads.

dir.create("TMM_normalized_reads")

write.table(
	merged, file.path("TMM_normalized_reads", "raw_counts.tsv"),
	sep="\t", col.names=T, row.names=F, quote=F
)

## Format reads into counts matrix.

merged <- merged %>%
	column_to_rownames("ensembl_id") %>%
	as.matrix

## TMM Normalize Reads
## ----------

## Create edgeR object.

groups <- c(rep(1,3), rep(2,3), rep(3,3))

edger.obj <- DGEList(counts=merged, group=groups)

# TMM normalize read counts.

edger.obj <- calcNormFactors(edger.obj, method="TMM")

## Get TMM normalized reads.

TMM <- edger.obj %>%
	cpm %>%
	as_tibble(rownames="ensembl_id")

## Export TMM normalized reads.

write.table(
	TMM, file.path("TMM_normalized_reads", "TMM_normalized_reads.tsv"),
	sep="\t", col.names=T, row.names=F, quote=F
)
