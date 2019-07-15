# Analysis of *S. cerevisiae* STRIPE-seq TSS Data

This directory contains the TSS analysis for yeast STRIPE-seq data.

## Contents

- **get_TSSs**: Retrieve TSSs from BAM files using TSRchitect.
- **TSS_correlation**: TMM normalized correlation of TSSs.
- **make_bedgraphs**: Create bedgraphs using TMM normalized counts.
- **annotate_TSSs**: Annotate TSSs to nearest transcript.
- **promoter_proximal**: Get fraction of TSSs that are promoter proximal while iterating through global threshold values.
