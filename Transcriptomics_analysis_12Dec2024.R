#################################
## source function
VolcanoScatter = function(DEtable, title, lFC=0, adjp = 0.1, nlabel=10){
  DEtable$Lipid = rownames(DEtable)
  DEtable = DEtable[order(DEtable$FDR),]
  p = ggplot(DEtable, aes(y=-log10(FDR), x=logFC, label=Lipid))+
    geom_point(color="darkgray") +
    geom_point(data=DEtable %>% filter(logFC >= lFC) %>% filter(FDR < adjp), col="#E64B35FF") +
    geom_point(data=DEtable %>% filter(logFC <= -lFC) %>% filter(FDR < adjp), col="#2E9FDF") +
    theme_bw() +
    geom_hline(yintercept = -log10(adjp), linetype="dotted", color="black", linewidth=1) +
    geom_vline(xintercept = lFC, linetype="dotted", color="black", linewidth=1) +
    geom_vline(xintercept = -lFC, linetype="dotted", color="black", linewidth=1) +
    geom_hline(yintercept = 0, linetype="dotted", color="grey", linewidth=1) +
    geom_vline(xintercept = 0, linetype="dotted", color="grey", linewidth=1) +
    geom_text_repel(data=DEtable %>% filter(logFC >= lFC) %>% slice_head(n=nlabel),max.overlaps=Inf,color="#E64B35FF") +
    geom_text_repel(data=DEtable %>% filter(logFC <= -lFC) %>% slice_head(n=nlabel),max.overlaps=Inf,color="#2E9FDF") +
    ggtitle(title)+
    scale_x_continuous(
      expand = c(0, 0), 
      limits = c(-4, 4)
    )
  return(p)
}

ScalebyRow <- function(x, na.rm=TRUE){
  rm <- rowMeans(x, na.rm = na.rm)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd, na.rm = na.rm)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

#################################
#install the packages
setwd("D:/研究生的科研/数据分析/IP_promgram/PE&GDM/manuscript/Transcriptomics/")
getwd()

library(edgeR)
library(limma)
library(ggplot2)
library(patchwork)
library(factoextra)
library(AnnotationDbi)
library(sva)
library(Matrix)
library(readxl)
library(dplyr)
library(tibble)
library(ggrepel)
library(scales)
library(ggsci)
library(tidyverse)
library(ggfortify)
library(gridExtra)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

#################################
### preparing the data
pal = c("#4E9040","#3D59A5","#F5931F","#ED3523")
CondCol = pal[c(1,2,3,4)]
names(CondCol) = c("CTL", "GDM", "PE","PG")

ExpCol<- c("#E64B35FF","#2E9FDF")
names(ExpCol) = c("up", "down")

ClusterCol = brewer.pal(6,"Set2")
ClusterCol

# Load the data
data_raw <- read.table("Transcriptomics_fcounts.txt", header=TRUE, row.names=1)
data <- data_raw[, 6:ncol(data_raw)]
colnames(data) <- gsub(".*STARout.(.*)A_.*", "\\1", colnames(data))
data <- data[, c("CTL1","CTL2","CTL3","CTL4","CTL5","GDM1","GDM2","GDM3","GDM4","GDM5",
                 "PE1","PE2","PE3","PG1","PG2","PG3","PG4")]

fcnt <- read.table("Transcriptomics_fcounts.txt.summary",row.names = 1, header=T, stringsAsFactors = FALSE, check.names = F)
colnames(fcnt)
fcnt <- fcnt[, c("CTL1","CTL2","CTL3","CTL4","CTL5","GDM1","GDM2","GDM3","GDM4","GDM5",
                 "PE1","PE2","PE3","PG1","PG2","PG3","PG4")]

meta_data <- read_xlsx("D:/研究生的科研/数据分析/IP_promgram/PE&GDM/manuscript/Sample_Information.xlsx")
meta_data <- meta_data[!is.na(meta_data$RNA_name),]
meta_data <- meta_data[match(colnames(fcnt), meta_data$Lipid_name), ]

# match the gene names
gene_anno <- read.table("gene_anno.txt", header=FALSE, sep="\t")
gene_anno$unique_gene_name <- make.unique(gene_anno[, 3], sep = "_")
rownames(data) <- gene_anno$unique_gene_name[match(rownames(data), gene_anno$V1)]
colnames(data)
colnames(data) <- meta_data$Lipid_name

#################################
### QC Analysis
meta_data$n_count = colSums(data)
meta_data$n_genes = apply(data, 2, function(x) sum(x>0))
meta_data$assigned_ratio = as.numeric(fcnt[1,meta_data$RNA_name]/colSums(fcnt))
meta_data$log_ng = log2(meta_data$n_genes)
meta_data$log_nc = log2(meta_data$n_count)

ggplot(meta_data, aes(x=n_count, y=n_genes)) + 
  geom_point(aes(color=Disease_Type), size=3) +
  geom_text_repel(aes(label=Lipid_name)) +
  theme_bw()+
  scale_color_manual(values = CondCol)+
  scale_fill_manual(values = CondCol)
ggsave("RNA_count_gene.png", width=6, height=5)

# delete the outlier GDM4 and CTL3
meta_data <- meta_data[-9,]
meta_data <- meta_data[-3,]
meta_data$Lipid_name
data <- data[,-9]
data <- data[,-3]
colnames(data)

#################################
### Sample normalization
y = DGEList(counts=data, group=meta_data$Disease_Type)
keep <- filterByExpr(y,group = meta_data$Disease_Type,min.count = 50, min.total.count = 200) #36601
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y) 
dim(y) #12859
logcpm = cpm(y, normalized.lib.sizes = TRUE, log = TRUE)

pheatmap(cor(logcpm),clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",filename = "RNA_correlation.png",width = 6, height = 5)

# PCA
pca1 = logcpm %>% 
  t() %>%
  prcomp(scale. = TRUE, center = TRUE)
autoplot(
  pca1,
  data = meta_data,
  x=1,
  y=2,
  fill = "Disease_Type",
  frame.colour = "Disease_Type",
  shape = 21,
  size = 3,
  alpha = 0.8,
  frame = TRUE,
)+
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  scale_fill_manual(values=CondCol) +
  scale_color_manual(values=CondCol) +
  theme(panel.grid.minor = element_blank()) +
  geom_text_repel(aes(label=Lipid_name))
ggsave("RNA_PCA.png", width=5, height=4)

#################################
### Differential expression analysis
nGenes = nrow(y) 
nrow(y) #12859

design <- model.matrix(~Disease_Type, data = meta_data)
y1 <- estimateGLMCommonDisp(y, design)
y1 <- estimateGLMTrendedDisp(y1, design)
y1 <- estimateGLMTagwiseDisp(y1, design)
fit <- glmQLFit(y1, design)

GDM_G <- topTags(glmQLFTest(fit, coef = 2), adjust.method = "fdr", n = nGenes)$table
dim(GDM_G)
GDM_G <- GDM_G[order(GDM_G$FDR),]
GDM_DEGs <- GDM_G[GDM_G$FDR < 0.25&(GDM_G$logFC>1|GDM_G$logFC<(-1)),] #1
dim(GDM_DEGs) #825
V1 <- VolcanoScatter(GDM_G, title="GDM vs CTL", lFC = log2(2), adjp = 0.25, nlabel=0)
V1
ggsave("RNA_Volcano_GDMvsCTL.png", width=4, height=4)

PE_G <- topTags(glmQLFTest(fit, coef = 3), adjust.method = "fdr", n = nGenes)$table
dim(PE_G)
PE_G <- PE_G[order(PE_G$FDR),]
PE_DEGs <- PE_G[PE_G$FDR < 0.25&(PE_G$logFC>1|PE_G$logFC<(-1)),] #2
dim(PE_DEGs) #60
V2 <- VolcanoScatter(PE_G, title="PE vs CTL", lFC = log2(2), adjp = 0.25, nlabel=0)
V2
ggsave("RNA_Volcano_PEvsCTL.png", width=4, height=4)

PG_G <- topTags(glmQLFTest(fit, coef = 4), adjust.method = "fdr", n = nGenes)$table
dim(PG_G)
PG_G <- PG_G[order(PG_G$FDR),]
PG_DEGs <- PG_G[PG_G$FDR < 0.25&(PG_G$logFC>1|PG_G$logFC<(-1)),] #1120
dim(PG_DEGs) #482
V3 <- VolcanoScatter(PG_G, title="PG vs CTL", lFC = log2(2), adjp = 0.25, nlabel=0)
V3
ggsave("RNA_Volcano_PGvsCTL.png", width=4, height=4)

write.csv(GDM_G, "GDMvsCTL_Genes.csv",row.names = T)
write.csv(PE_G, "PEvsCTL_Genes.csv",row.names = T)
write.csv(PG_G, "PGvsCTL_Genes.csv",row.names = T)

###############################
## visualization of plot
library(venn)
png("venn_overlap.png", width = 5, height = 5, units = "in", res = 300)
venn(list(
  GDM = rownames(GDM_DEGs),
  PG = rownames(PG_DEGs),
  PE = rownames(PE_DEGs)
),
ilab = c("GDM", "PG", "PE"),
col = "transparent",
ilabels = "counts",
zcolor = c("#4E9040","#3D59A5","#F5931F"),
ilcs = 1.5, sncs = 1.8)
dev.off()

## heatmap
ALL_EXP = as.data.frame(logcpm[rownames(logcpm) %in% unique(c(rownames(GDM_DEGs),rownames(PE_DEGs),rownames(PG_DEGs))),])
rownames(gene_anno) <- gene_anno$unique_gene_name
gene_protein <- gene_anno[gene_anno$V2 == "protein_coding",]
ALL_EXP <- rownames_to_column(ALL_EXP, "GeneID")
ALL_EXP <- left_join(ALL_EXP,gene_protein,by=c("GeneID"="unique_gene_name"))

dim(ALL_EXP[!is.na(ALL_EXP$V3),])
ALL_EXP$CTL = apply(ALL_EXP[,c("CTL1", "CTL2", "CTL4", "CTL5")], 1, median)
ALL_EXP$GDM = apply(ALL_EXP[,c("GDM1", "GDM2", "GDM3", "GDM5")], 1, median)
ALL_EXP$PE = apply(ALL_EXP[,c("PE1", "PE2", "PE3")], 1, median)
ALL_EXP$PG = apply(ALL_EXP[,c("PG1", "PG2", "PG3", "PG4")], 1, median)

GDM_DEGs<- rownames_to_column(GDM_DEGs, "GeneID")
PE_DEGs<- rownames_to_column(PE_DEGs, "GeneID")
PG_DEGs<- rownames_to_column(PG_DEGs, "GeneID")

GDM_DEGs$GDM_DEGs = ifelse(GDM_DEGs$logFC>1,"up","down")
PE_DEGs$PE_DEGs = ifelse(PE_DEGs$logFC>1,"up","down")
PG_DEGs$PG_DEGs = ifelse(PG_DEGs$logFC>1,"up","down")

ALL_EXP <- left_join(ALL_EXP, GDM_DEGs[,c("GeneID","GDM_DEGs")], by = "GeneID")
ALL_EXP <- left_join(ALL_EXP, PE_DEGs[,c("GeneID","PE_DEGs")], by = "GeneID")
ALL_EXP <- left_join(ALL_EXP, PG_DEGs[,c("GeneID","PG_DEGs")], by = "GeneID")

ALL_ES = ScalebyRow(ALL_EXP[!is.na(ALL_EXP$V3),c("CTL","GDM","PE","PG")])
rownames(ALL_ES) = ALL_EXP[!is.na(ALL_EXP$V3),"GeneID"]
ALL_ES[ALL_ES > 3] = 3
ALL_ES[ALL_ES < -3] = -3

pheatmap_res <- pheatmap(ALL_ES,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = T,
         cutree_rows = 6,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2"
         )
row_clusters <- cutree(pheatmap_res$tree_row, k = 6)
# write.csv(row_clusters, "RNA_DEGs_clusters.csv")

names(ClusterCol) <- unique(row_clusters)
annotation_row <- data.frame(row.names = rownames(ALL_ES),
                             Cluster = as.factor(row_clusters),
                             GDM_DEGs = factor(ALL_EXP[!is.na(ALL_EXP$V3),"GDM_DEGs"], levels=c("up","down")),
                             PE_DEGs = factor(ALL_EXP[!is.na(ALL_EXP$V3),"PE_DEGs"], levels=c("up","down")),
                             PG_DEGs = factor(ALL_EXP[!is.na(ALL_EXP$V3),"PG_DEGs"], levels=c("up","down")))
annotation_col = data.frame(row.names = colnames(ALL_ES),Disease_Type = factor(c("CTL","GDM","PE","PG"), levels=c("CTL","GDM","PE","PG")))
pheatmap(ALL_ES,
         annotation_col = annotation_col,
         annotation_colors = list(Disease_Type = CondCol,GDM_DEGs=ExpCol,PE_DEGs=ExpCol,PG_DEGs=ExpCol,Cluster=ClusterCol),
         annotation_row = annotation_row, 
         cluster_cols = FALSE, 
         show_rownames = FALSE,
         show_colnames = TRUE,
         cutree_rows = 6,  
         clustering_distance_rows = "correlation", 
         clustering_method = "ward.D2",  
         filename = "heatmap_median_cluster.png",
         width = 6, height = 8
         )

ALL_pro_gene <- ALL_EXP[!is.na(ALL_EXP$V3),]
ALL_pro_gene$Cluster <- row_clusters[match(ALL_EXP[!is.na(ALL_EXP$V3),]$GeneID, names(row_clusters))]
write.csv(ALL_pro_gene, "RNA_DEGs_ALL_Cluster.csv")

##################################
## cluster of up and down trend - cluster4 and cluster5
# cluster4
dim(ALL_pro_gene[ALL_pro_gene$Cluster %in% c("4"),]) #22
rownames(ALL_pro_gene) <- ALL_pro_gene$GeneID

Cluster4 <- ScalebyRow(ALL_pro_gene[ALL_pro_gene$Cluster %in% c("4"),2:16])
rownames(Cluster4) <- ALL_pro_gene[ALL_pro_gene$Cluster %in% c("4"),]$GeneID
annotation_col <- data.frame(row.names = colnames(Cluster4),Disease_Type = factor(meta_data$Disease_Type,levels=c("CTL","GDM","PE","PG")))
annotation_row4 <- annotation_row[rownames(Cluster4),]
pheatmap(Cluster4,
         annotation_col = annotation_col,
         annotation_colors = list(Disease_Type = CondCol,GDM_DEGs=ExpCol,PE_DEGs=ExpCol,PG_DEGs=ExpCol,Cluster=ClusterCol[4]),
         annotation_row = annotation_row4, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend = FALSE,
         annotation_legend = FALSE,
         filename = "RNA_heatmap_cluster4.png",
         width = 4, height = 3.8
         )

# cluster5
dim(ALL_pro_gene[ALL_pro_gene$Cluster %in% c("5"),]) #60
Cluster5 <- ScalebyRow(ALL_pro_gene[ALL_pro_gene$Cluster %in% c("5"),2:16])
rownames(Cluster5) <- ALL_pro_gene[ALL_pro_gene$Cluster %in% c("5"),]$GeneID
annotation_row5 <- annotation_row[rownames(Cluster5),]
pheatmap(Cluster5,
         annotation_col = annotation_col,
         annotation_colors = list(Disease_Type = CondCol,GDM_DEGs=ExpCol,PE_DEGs=ExpCol,PG_DEGs=ExpCol,Cluster=ClusterCol[5]),
         annotation_row = annotation_row5, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend = FALSE,
         annotation_legend = FALSE,
         filename = "RNA_heatmap_cluster5.png",
         width = 4, height = 8.5
)

##################################
## Cluster of Enrichment Pathways
Enrich_cluster <- read_xlsx("RNA_pathway_cluster.xlsx")
Enrich_cluster$Cluster <- factor(Enrich_cluster$Cluster,levels = c("4","6","2","5","3","1"))
ggplot(Enrich_cluster, aes(x = `Odds Ratio`, y = reorder(Term, `Odds Ratio`), fill = Cluster)) +
  geom_bar(stat = "identity") + 
  labs(title = NULL, y = NULL, x = "Odds Ratio", fill = "Cluster") +  
  scale_fill_manual(  
    values = ClusterCol  
  ) +
  facet_grid(Cluster ~ ., scales = "free_y", space = "free_y", switch = "x") +
  #scale_x_continuous(expand = c(0.05, 0.05)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        legend.position = "none") +
  scale_y_discrete(position = "right")
ggsave("RNA_Enrichpathway_cluster.png", width = 7, height = 4.5)

############################
## Gene expression
list3 <- c("GPER1","PDK4","AHR","KANK1","S100A10")
## bar plot of exp
v_plot <- as.data.frame(t(logcpm[list3,]))
v_plot$Disease_Type <- meta_data$Disease_Type
v_plot_long <- pivot_longer(v_plot, cols = -Disease_Type,
                            names_to = "Gene", values_to = "Expression")

B3 <- ggplot(v_plot_long, aes(x = Gene, y = Expression, fill = Disease_Type)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = CondCol) + 
  theme_bw() +
  labs(y = "Expression Level (logCPM)", x = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
B3
ggsave("Barplot.png", plot = B3, width = 6, height = 3)

#######################
# exp of biomarker genes other articles
logcpm_t <- as.data.frame(t(logcpm))
logcpm_t$Disease_Type <- meta_data$Disease_Type
write.csv(logcpm_t, "logcpm_t.csv")

# selected relative genes with literature reported
plots_violin <- list()
gene1 <- c("CAMK1D","IL1B", #GDM
           "RGL3","PREX1","MME","GPX3", #PE
           "TTN","PTGS2") #Outcome 
for (i in seq_along(gene1)) {
  plot <- ggplot(logcpm_t, aes(x = Disease_Type, y = .data[[gene1[i]]], fill = Disease_Type)) +
    geom_violin(trim = FALSE, width = 0.8, scale = "width") +  
    geom_boxplot(width = 0.3, position = position_dodge(0.9),fill="white") +  
    scale_fill_manual(values = CondCol) +  
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +
    labs(y = "Relative Expression",x = NULL,title = gene1[i]) +
    stat_compare_means(comparisons = list(c("GDM", "CTL"), c("PE", "CTL"), c("PG", "CTL")),method = "t.test", label = "p.signif") 
  plots_violin[[i]] <- plot
  ggsave(paste0("RNA_Reported_Violin_", gene1[i], ".png"), plot, width = 3, height = 3.5, dpi = 300)
}

# exp of violin the cluster of 2 4 6 significant of disease
gene2 <- c("ZNF683","PPDPF", #PG related gene
           "B3GNT7","KLRC3","ENOSF1", #GDM related gene
           "MERTK","CD163","MXRA7") #PE related gene
plots_violin <- list()
for (i in seq_along(gene2)) {
  plot <- ggplot(logcpm_t, aes(x = Disease_Type, y = .data[[gene2[i]]], fill = Disease_Type)) +
    geom_violin(trim = FALSE, width = 0.8, scale = "width") +  
    geom_boxplot(width = 0.3, position = position_dodge(0.9),fill="white") +  
    scale_fill_manual(values = CondCol) +  
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +
    labs(y = "Relative Expression",x = NULL,title = gene2[i]) +
    stat_compare_means(comparisons = list(c("GDM", "CTL"), c("PE", "CTL"), c("PG", "CTL")),method = "t.test", label = "p.signif") 
  plots_violin[[i]] <- plot
  ggsave(paste0("RNA_Disease_Violin_", gene2[i], ".png"), plot, width = 3, height = 3.5, dpi = 300)
}

## viloin plot of the selecte genes
gene3 <- c("LRRC56", "SPSB2", "MPO", "ELANE", "TNNT1", "SIK1B",
           "EPHA4", "GFRA2", "PDZD4")
## without "TNN"
plots_violin <- list()
for (i in seq_along(gene3)) {
  plot <- ggplot(logcpm_t, aes(x = Disease_Type, y = .data[[gene3[i]]], fill = Disease_Type)) +
    geom_violin(trim = FALSE, width = 0.8, scale = "width") +  
    geom_boxplot(width = 0.3, position = position_dodge(0.9),fill="white") +  
    scale_fill_manual(values = CondCol) +  
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +
    labs(y = "Relative Expression",x = NULL,title = gene3[i]) +
    stat_compare_means(comparisons = list(c("GDM", "CTL"), c("PE", "CTL"), c("PG", "CTL")),method = "t.test", label = "p.signif") 
  plots_violin[[i]] <- plot
  ggsave(paste0("RNA_Disease_Violin_", gene3[i], ".png"), plot, width = 3, height = 3, dpi = 300)
}
