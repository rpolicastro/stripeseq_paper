#!/usr/env/bin Rscript

library(Gviz)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

#create genomic axis track
axis.track <- GenomeAxisTrack(col="black",
                              scale=0.1,
                              col.range="black")

options(ucscChromosomeNames=FALSE)

#create gene annotation track
genome.track <- GeneRegionTrack(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,
                                genome="sacCer3",
                                shape="arrow",
                                names="Genes",
                                col="black",
                                showId=TRUE,
                                fill="black",
                                thinBoxFeature=c("utr",
                                                 "ncRNA", 
                                                 "utr3", 
                                                 "utr5", 
                                                 "3UTR", 
                                                 "5UTR", 
                                                 "miRNA", 
                                                 "lincRNA", 
                                                 "three_prime_UTR", 
                                                 "five_prime_UTR"),
                                trancriptAnnotation="gene_symbol",
                                collapseTranscripts="meta")

#load data tracks and plot at HHT1 (narrow TSR)
S288C_unpooled.1.pos.track <- DataTrack(range="pos_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkblue",
                                        fill.histogram="darkblue",
                                        ylim=c(0,90),
                                        baseline=0,
                                        col.baseline="darkblue")

S288C_unpooled.1.min.track <- DataTrack(range="min_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkcyan",
                                        fill.histogram="darkcyan",
                                        ylim=c(0,1),
                                        baseline=0,
                                        col.baseline="darkcyan")

pdf(file="tracks/HHT1.pdf")
plotTracks(list(axis.track,
                S288C_unpooled.1.pos.track,
                S288C_unpooled.1.min.track,
                genome.track),
           chromosome="chrII",
           from=256250,
           to=256800,
           background.title="white",
           col.title="black",
           col.axis="black")
dev.off()

#load data tracks and plot at RPS8B (multiple TSRs)
S288C_unpooled.1.pos.track <- DataTrack(range="pos_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkblue",
                                        fill.histogram="darkblue",
                                        ylim=c(0,400),
                                        baseline=0,
                                        col.baseline="darkblue")

S288C_unpooled.1.min.track <- DataTrack(range="min_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkcyan",
                                        fill.histogram="darkcyan",
                                        ylim=c(0,1),
                                        baseline=0,
                                        col.baseline="darkcyan")

pdf(file="tracks/RPS8B.pdf")
plotTracks(list(axis.track,
                S288C_unpooled.1.pos.track,
                S288C_unpooled.1.min.track,
                genome.track),
           chromosome="chrV",
           from=362500,
           to=363900,
           background.title="white",
           col.title="black",
           col.axis="black")
dev.off()

#load data tracks and plot at PEX22 (TSR includes start codon)
S288C_unpooled.1.pos.track <- DataTrack(range="pos_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkblue",
                                        fill.histogram="darkblue",
                                        ylim=c(0,55),
                                        baseline=0,
                                        col.baseline="darkblue")

S288C_unpooled.1.min.track <- DataTrack(range="min_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkcyan",
                                        fill.histogram="darkcyan",
                                        ylim=c(0,1),
                                        baseline=0,
                                        col.baseline="darkcyan")

pdf(file="tracks/PEX22.pdf")
plotTracks(list(axis.track,
                S288C_unpooled.1.pos.track,
                S288C_unpooled.1.min.track,
                genome.track),
           chromosome="chrI",
           from=42000,
           to=43000,
           background.title="white",
           col.title="black",
           col.axis="black")
dev.off()

#load data tracks and plot at AIM39 (TSR downstream of annotated start codon)
S288C_unpooled.1.pos.track <- DataTrack(range="pos_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkblue",
                                        fill.histogram="darkblue",
                                        ylim=c(0,60),
                                        baseline=0,
                                        col.baseline="darkblue")

S288C_unpooled.1.min.track <- DataTrack(range="min_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkcyan",
                                        fill.histogram="darkcyan",
                                        ylim=c(0,2.5),
                                        baseline=0,
                                        col.baseline="darkcyan")

pdf(file="tracks/AIM39.pdf")
plotTracks(list(axis.track,
                S288C_unpooled.1.pos.track,
                S288C_unpooled.1.min.track,
                genome.track),
           chromosome="chrXV",
           from=230000,
           to=231400,
           background.title="white",
           col.title="black",
           col.axis="black")
dev.off()

#load data tracks and plot at GCN4 (upstream TSR)
S288C_unpooled.1.pos.track <- DataTrack(range="pos_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkblue",
                                        fill.histogram="darkblue",
                                        ylim=c(0,2.5),
                                        baseline=0,
                                        col.baseline="darkblue")

S288C_unpooled.1.min.track <- DataTrack(range="min_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkcyan",
                                        fill.histogram="darkcyan",
                                        ylim=c(0,100),
                                        baseline=0,
                                        col.baseline="darkcyan")

pdf(file="tracks/GCN4.pdf")
plotTracks(list(axis.track,
                S288C_unpooled.1.pos.track,
                S288C_unpooled.1.min.track,
                genome.track),
           chromosome="chrV",
           from=138800,
           to=140500,
           background.title="white",
           col.title="black",
           col.axis="black")
dev.off()

#load data tracks and plot at CCT8 (sense/antisense TSRs)
S288C_unpooled.1.pos.track <- DataTrack(range="pos_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkblue",
                                        fill.histogram="darkblue",
                                        ylim=c(0,25),
                                        baseline=0,
                                        col.baseline="darkblue")

S288C_unpooled.1.min.track <- DataTrack(range="min_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkcyan",
                                        fill.histogram="darkcyan",
                                        ylim=c(0,65),
                                        baseline=0,
                                        col.baseline="darkcyan")

pdf(file="tracks/CCT8.pdf")
plotTracks(list(axis.track,
                S288C_unpooled.1.pos.track,
                S288C_unpooled.1.min.track,
                genome.track),
           chromosome="chrX",
           from=419750,
           to=423000,
           background.title="white",
           col.title="black",
           col.axis="black")
dev.off()

#load data tracks and plot at RPS15 (broad TSR)
S288C_unpooled.1.pos.track <- DataTrack(range="pos_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkblue",
                                        fill.histogram="darkblue",
                                        ylim=c(0,2.5),
                                        baseline=0,
                                        col.baseline="darkblue")

S288C_unpooled.1.min.track <- DataTrack(range="min_S288C-unpooled_WT-100ng_1.bedgraph", 
                                        genome="sacCer3", 
                                        name="S288C pooled rep 1",
                                        type="hist",
                                        col.histogram="darkcyan",
                                        fill.histogram="darkcyan",
                                        ylim=c(0,300),
                                        baseline=0,
                                        col.baseline="darkcyan")

pdf(file="tracks/RPS15.pdf")
plotTracks(list(axis.track,
                S288C_unpooled.1.pos.track,
                S288C_unpooled.1.min.track,
                genome.track),
           chromosome="chrXV",
           from=253100,
           to=253700,
           background.title="white",
           col.title="black",
           col.axis="black")
dev.off()