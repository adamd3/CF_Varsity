#!/usr/bin/env Rscript

##------------------------------------------------------------------------------
## Load libraries
##------------------------------------------------------------------------------
library(ComplexHeatmap)
library(DESeq2)
library(ggplot2)
library(NbClust)
library(org.Hs.eg.db)
library(plyr)
library(RColorBrewer)
library(readr)
library(stringr)
library(topGO)


##------------------------------------------------------------------------------
## Parameters
##------------------------------------------------------------------------------
workingDir <- "./"
set1Pal <- brewer.pal(9, "Set1") ## colour palette for plots


##------------------------------------------------------------------------------
## Functions
##------------------------------------------------------------------------------
## GO term enrichment using topGO
topGO_enrich <- function(clu_number){
    clu_number <- as.character(clu_number)
    cluster_ids <- rownames(FPKM_NE_zscores)[
        median_zscores_NE$Gene_cluster==clu_number]
    sig_idx <- match(cluster_ids, rownames(res_NE_CF0))
    backG <- setdiff(rownames(res_NE_CF0), cluster_ids)
    inSelection <- rownames(res_NE_CF0) %in% cluster_ids
    alg <- factor(as.integer(rownames(res_NE_CF0) %in% cluster_ids))
    names(alg) <- rownames(res_NE_CF0)
    tgd <- new(
        "topGOdata", ontology = "BP", allGenes = alg, nodeSize = 10,
        annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl"
    )
    resultTopGO_classic <- runTest(
        tgd, algorithm = "classic", statistic = "Fisher")
    topGOres <- GenTable(tgd, Fisher.classic = resultTopGO_classic,
        orderBy = "Fisher.classic", topNodes = length(usedGO(tgd))
    )
    topGOres <- rbind.fill(topGOres)
    topGOres$Fisher.classic <- as.numeric(topGOres$Fisher.classic)
    topGOres$Fisher_padj <- p.adjust(topGOres$"Fisher.classic", method = "BH")
    topGOres <- topGOres[order(as.numeric(topGOres$Fisher.classic)),]
    topGOres
}




##------------------------------------------------------------------------------
## Read in data
##------------------------------------------------------------------------------
## sample tables
samplesDf_NE <- read_tsv(file.path("sample_tables","samplesDf_NE.tsv"))
samplesDf_Mono <- read_tsv(file.path("sample_tables","samplesDf_Mono.tsv"))

## read counts per gene
gene_counts_NE <- read.table(file.path(
    "gene_counts","gene_counts_NE.tsv"), header=TRUE, fill=TRUE)
gene_counts_Mono <- read.table(file.path(
    "gene_counts","gene_counts_Mono.tsv"), header=TRUE, fill=TRUE)

## FPKM values per gene
FPKM_NE <- read.table(file.path(
    "gene_counts","FPKM_NE.tsv"), header=TRUE, fill=TRUE)
FPKM_Mono <- read.table(file.path(
    "gene_counts","FPKM_Mono.tsv"), header=TRUE, fill=TRUE)

## deseq2 results
res_NE_CF0 <- readRDS(file.path("DESeq2_res","res_NE_CF0.rds"))
res_NE_CF14 <- readRDS(file.path("DESeq2_res","res_NE_CF14.rds"))
res_NE_CF30 <- readRDS(file.path("DESeq2_res","res_NE_CF30.rds"))
res_Mono_CF0 <- readRDS(file.path("DESeq2_res","res_Mono_CF0.rds"))
res_Mono_CF14 <- readRDS(file.path("DESeq2_res","res_Mono_CF14.rds"))
res_Mono_CF30 <- readRDS(file.path("DESeq2_res","res_Mono_CF30.rds"))


## match columns in gene counts with rows in sample table
gene_counts_NE <- gene_counts_NE[
    , match(samplesDf_NE$Sample_ID, colnames(gene_counts_NE))
]
gene_counts_Mono <- gene_counts_Mono[
    , match(samplesDf_Mono$Sample_ID, colnames(gene_counts_Mono))
]

FPKM_NE <- FPKM_NE[
    , match(samplesDf_NE$Sample_ID, colnames(FPKM_NE))
]
FPKM_Mono <- FPKM_Mono[
    , match(samplesDf_Mono$Sample_ID, colnames(FPKM_Mono))
]



##------------------------------------------------------------------------------
## Clustering + heatmap of DEGs in Neutrophils
##------------------------------------------------------------------------------
DE_NE_RNA_CF0 <- unique(na.omit(rownames(
    subset(res_NE_CF0, padj < 0.05 & abs(log2FoldChange) > 0.5))))
DE_NE_RNA_CF14 <- unique(na.omit(rownames(
    subset(res_NE_CF14, padj < 0.05 & abs(log2FoldChange) > 0.5))))
DE_NE_RNA_CF30 <- unique(na.omit(rownames(
    subset(res_NE_CF30, padj < 0.05 & abs(log2FoldChange) > 0.5))))

DE_NE_RNA_all <- unique(c(
    DE_NE_RNA_CF0, DE_NE_RNA_CF14, DE_NE_RNA_CF30
))

DE_NE_RNA_all[grep("ENSG", DE_NE_RNA_all)] <- gsub(
    "\\..*","",DE_NE_RNA_all[grep("ENSG", DE_NE_RNA_all)])

FPKM_NE_sub <- subset(FPKM_NE, rownames(FPKM_NE) %in% DE_NE_RNA_all)

## calculate z-scores per gene
FPKM_NE_zscores <- t(scale(t(as.matrix(FPKM_NE_sub))))

## median z-score by group
median_zscores_NE <- data.frame(
    Day_0 = matrixStats::rowMedians(
        FPKM_NE_zscores[,(
            samplesDf_NE$Class=="CF" & samplesDf_NE$Timepoint=="0")]
    ),
    Day_14 = matrixStats::rowMedians(
        FPKM_NE_zscores[,(
            samplesDf_NE$Class=="CF" & samplesDf_NE$Timepoint=="14")]
    ),
    Day_30_plus = matrixStats::rowMedians(
        FPKM_NE_zscores[,(
            samplesDf_NE$Class=="CF" & samplesDf_NE$Timepoint=="30")]
    ),
    HV = matrixStats::rowMedians(
        FPKM_NE_zscores[,(
            samplesDf_NE$Class=="HV")])
)


## hierarchical clustering
hclust_genes <- hclust(
    dist(median_zscores_NE, method = 'manhattan'),
    method='ward.D2'
)


## choose number of clusters using C-index method
res_nbclust <- NbClust(
    median_zscores_NE, distance = "manhattan",
    min.nc = 2, max.nc = 10, method = "ward.D2", index = "cindex"
)
nclust <- res_nbclust$Best.nc[[1]]


## cut the tree:
geneclusters <- cutree(hclust_genes, k = nclust)


median_zscores_NE$Gene_cluster <- as.factor(geneclusters)


## Plot heatmap
ann_colors_NE = list(
    Gene_cluster = c(
        "1" = set1Pal[1], "2" = set1Pal[2],
        "3" = set1Pal[3], "4" = set1Pal[4],
        "5" = set1Pal[5], "6" = set1Pal[6],
        "7" = set1Pal[7], "8" = set1Pal[8]
    )
)



png("NE_pheatmap.png",
    width = 1100, height = 2700,
    res = 300
)
ComplexHeatmap::pheatmap(
    name = "z-score",
    as.matrix(median_zscores_NE[,1:(ncol(median_zscores_NE)-1)]),
    annotation_row = subset(median_zscores_NE, select=Gene_cluster),
    fontsize = 16, annotation_legend = c(TRUE, TRUE),
    legend = c(TRUE), cluster_rows = FALSE, cluster_cols = FALSE,
    border_color = NA,
    show_rownames = TRUE, show_colnames = TRUE,
    annotation_colors = ann_colors_NE,
    annotation_names_row = FALSE, annotation_names_col = FALSE,
    row_split = median_zscores_NE$Gene_cluster
)
dev.off()



##------------------------------------------------------------------------------
## Functional enrichment of gene clusters in Neutrophils
##------------------------------------------------------------------------------
topGO_res <- lapply(unique(median_zscores_NE$Gene_cluster), topGO_enrich)

## plot enriched terms
nkeep <- 5

combinedRes <- rbind(
    topGO_res[[1]][1:nkeep,], topGO_res[[2]][1:nkeep,],
    topGO_res[[3]][1:nkeep,], topGO_res[[4]][1:nkeep,],
    topGO_res[[5]][1:nkeep,], topGO_res[[6]][1:nkeep,],
    topGO_res[[7]][1:nkeep,], topGO_res[[8]][1:nkeep,]
)


combinedRes$Cluster <- c(
    rep("1", nkeep), rep("2", nkeep), rep("3", nkeep), rep("4", nkeep),
    rep("5", nkeep), rep("6", nkeep), rep("7", nkeep), rep("8", nkeep)
)
combinedRes$Cluster <- as.factor(combinedRes$Cluster)



## barplot of top n terms by cluster
topGOdf <- combinedRes
topGOdf$Term <- str_trunc(topGOdf$Term, 40)
topGOdf$Term <- factor(topGOdf$Term, levels = (topGOdf$Term))
topGOdf$Fisher_padj <- -log10(topGOdf$Fisher_padj)

topGOdf$Cluster <- factor(
    topGOdf$Cluster, levels = (sort(unique(topGOdf$Cluster)))
)


p1 <- ggplot(topGOdf,
    aes(x = Term, y = Fisher_padj, fill = Cluster)) +
    geom_bar(
        position = "dodge", stat = "identity",
        width = 0.7,  colour = "black"
    ) +
    geom_hline(yintercept=1, linetype="dashed") +
    theme_bw() +
    ylab(expression("−log"[10]~"("*italic(P)~"adj."*")")) +
    coord_flip() +
    scale_fill_manual("", values = set1Pal) +
    theme(
        legend.position="top",
        legend.title.align = 0,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(colour = "black", face = "bold", size = 18),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.text = element_text(colour = "black", face="bold", size = 18),
        legend.title = element_text(colour = "black", face="bold",size = 18)
    )

ggsave(
    p1, file = "topGO_bar.png",
    device = "png",
    width = 7, height = 12,
    dpi = 400
)



##------------------------------------------------------------------------------
## Clustering + heatmap of DEGs in Monocytes
##------------------------------------------------------------------------------
## subset to DEGs:
DE_Mono_RNA_CF0 <- unique(na.omit(rownames(
    subset(res_Mono_CF0, padj < 0.05 & abs(log2FoldChange) > 0.5))))
DE_Mono_RNA_CF14 <- unique(na.omit(rownames(
    subset(res_Mono_CF14, padj < 0.05 & abs(log2FoldChange) > 0.5))))
DE_Mono_RNA_CF30 <- unique(na.omit(rownames(
    subset(res_Mono_CF30, padj < 0.05 & abs(log2FoldChange) > 0.5))))

DE_Mono_RNA_all <- unique(c(
    DE_Mono_RNA_CF0, DE_Mono_RNA_CF14, DE_Mono_RNA_CF30
))


FPKM_Mono_sub <- subset(FPKM_Mono, rownames(FPKM_Mono) %in% DE_Mono_RNA_all)

## calculate z-scores per gene
FPKM_Mono_zscores <- t(scale(t(as.matrix(FPKM_Mono_sub))))

## median z-score by group
median_zscores_Mono <- data.frame(
    Day_0 = matrixStats::rowMedians(
        FPKM_Mono_zscores[,(
            samplesDf_Mono$Class=="CF" & samplesDf_Mono$Timepoint=="0")]
    ),
    Day_14 = matrixStats::rowMedians(
        FPKM_Mono_zscores[,(
            samplesDf_Mono$Class=="CF" & samplesDf_Mono$Timepoint=="14")]
    ),
    Day_30_plus = matrixStats::rowMedians(
        FPKM_Mono_zscores[,(
            samplesDf_Mono$Class=="CF" & samplesDf_Mono$Timepoint=="30")]
    ),
    HV = matrixStats::rowMedians(
        FPKM_Mono_zscores[,(samplesDf_Mono$Class=="HV")])
)
dim(median_zscores_Mono)
# [1] 33  4


groups_df <- data.frame(
    Sample = factor(
        c("Day_0", "Day_14", "Day_30_plus", "HV")
    )
)
rownames(groups_df) <- colnames(median_zscores_Mono)

ann_colors_Mono = list(
    Sample = c(
        "Day_0" = "#F8766D", "Day_14" = "#00BA38",
        "Day_30_plus" = "#619CFF", "HV" = "black"
    )
)

colnames(median_zscores_Mono) <- gsub("_"," ", colnames(median_zscores_Mono))

png("Mono_pheatmap.png",width = 600, height = 1100, res = 200)
ComplexHeatmap::pheatmap(
    name = "z-score",
    as.matrix(median_zscores_Mono),
    fontsize = 13,
    clustering_distance_rows = "manhattan",
    clustering_distance_cols = "manhattan",
    clustering_method = "ward.D2",
    annotation_legend = c(TRUE, TRUE),
    legend = c(TRUE),
    cluster_rows = TRUE, cluster_cols = TRUE,
    border_color = "grey60",
    show_rownames = FALSE, show_colnames = TRUE,
    labels_row = NA,
    annotation_colors = ann_colors_Mono,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE#
)
dev.off()
