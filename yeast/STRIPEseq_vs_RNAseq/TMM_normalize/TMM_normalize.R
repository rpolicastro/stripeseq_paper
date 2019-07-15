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
	"../count_files/Yeast_RNAseq_counts.tsv",
	sep="\t", header=T, stringsAsFactors=F
) %>% as_tibble

## Clean up column names.

colnames(rna.seq) <- str_replace_all(
	colnames(rna.seq),
	pattern="(X\\..*aligned\\.|_Aligned.*\\.bam)",
	replacement=""
)

## Rename columns.

rna.seq <- rna.seq %>%
	dplyr::rename(
		"sgd_id"=Geneid,
		"yeast_rnaseq_1"=YRNASEQ_01_S288C_1,
		"yeast_rnaseq_2"=YRNASEQ_02_S288C_2,
		"yeast_rnaseq_3"=YRNASEQ_03_S288C_3
	) %>%
	dplyr::select(
		sgd_id,
		yeast_rnaseq_1,
		yeast_rnaseq_2,
		yeast_rnaseq_3
	)

## Load STRIPE-seq Data
## ----------

## Load counts file.

stripe.seq <- read.delim(
	"../count_files/Yeast_STRIPEseq_counts.tsv",
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
		"sgd_id"=Geneid,
		"yeast_stripeseq_pooled_2"=S288C.pooled_WT.100ng_2,
		"yeast_stripeseq_pooled_1"=S288C.pooled_WT.100ng_1,
		"yeast_stripeseq_pooled_3"=S288C.pooled_WT.100ng_3,
		"yeast_stripeseq_unpooled_3"=S288C.unpooled_WT.100ng_3,
		"yeast_stripeseq_unpooled_2"=S288C.unpooled_WT.100ng_2,
		"yeast_stripeseq_unpooled_1"=S288C.unpooled_WT.100ng_1
	) %>%
	dplyr::select(
		sgd_id,
		yeast_stripeseq_unpooled_1,
		yeast_stripeseq_unpooled_2,
		yeast_stripeseq_unpooled_3,
		yeast_stripeseq_pooled_1,
		yeast_stripeseq_pooled_2,
		yeast_stripeseq_pooled_3
	)

## Merge and Prepare Data
## ----------

## Merge data.

merged <- left_join(rna.seq, stripe.seq, by="sgd_id")

## Export raw reads.

dir.create("TMM_normalized_reads")

write.table(
	merged, file.path("TMM_normalized_reads", "raw_counts.tsv"),
	sep="\t", col.names=T, row.names=F, quote=F
)

## Format reads into counts matrix.

merged <- merged %>%
	column_to_rownames("sgd_id") %>%
	as.matrix

## TMM Normalize Reads
## ----------

## Create edgeR object.

groups <- c(rep(1,3), rep(2,3), rep(3,3))

edger.obj <- DGEList(counts=merged, group=groups)

## TMM normalize read counts.

edger.obj <- calcNormFactors(edger.obj, method="TMM")

## Get TMM normalized reads.

TMM <- edger.obj %>%
	cpm %>%
	as_tibble(rownames="sgd_id")

## Export TMM normalized reads.

write.table(
	TMM, file.path("TMM_normalized_reads", "TMM_normalized_reads.tsv"),
	sep="\t", col.names=T, row.names=F, quote=F
)
