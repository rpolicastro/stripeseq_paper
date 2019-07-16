# Analysis of *S. cerevisiae* STRIPE-seq TSS Data

This directory contains the TSS analysis for yeast STRIPE-seq data.

## Contents

- **get_TSSs**: Retrieve TSSs from BAM files using TSRchitect.
- **TSS_correlation**: TMM normalized correlation of TSSs.
- **make_bedgraphs**: Create bedgraphs using TMM normalized counts.
- **annotate_TSSs**: Annotate TSSs to nearest transcript.
- **TSS_genomic_distribution**: Get fraction of TSSs that are promoter proximal while iterating through global threshold values.
- **TSS_average_plot**: Average plots of TSSs cenetered on start codon (for yeast data).

## Run Analysis.

To redo the execution of the analysis enter `bash RUN.sh`.
