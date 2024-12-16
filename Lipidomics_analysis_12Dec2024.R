setwd("D:/研究生的科研/数据分析/IP_promgram/PE&GDM/manuscript/Lipidomics/")
getwd()

#############################
## source function
ScalebyRow <- function(x, na.rm=TRUE){
  rm <- rowMeans(x, na.rm = na.rm)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd, na.rm = na.rm)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

Volcano_plot <- function(data, title = "Volcano Plot", pval_cutoff = 0.1, 
                         fc_cutoff = 0, nlabel = 10) {
  fold_change_col <- "logFC"      
  pval_col <- "adj.P.Val"          
  
  data$logPval <- -log10(data[[pval_col]])
  data$threshold <- ifelse(data[[fold_change_col]] > fc_cutoff & data[[pval_col]] < pval_cutoff, "up", 
                           ifelse(data[[fold_change_col]] < -fc_cutoff & data[[pval_col]] < pval_cutoff, "down", "no"))
  top_genes <- head(order(data$logPval, decreasing = TRUE), nlabel)
  highlight_genes <- rownames(data)[top_genes]
  max <- max(abs(min(data$logFC, na.rm = TRUE)), abs(max(data$logFC, na.rm = TRUE)))
  
  p <- ggplot(data, aes(x = logFC, y = logPval, color = threshold)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#2E9FDF", "darkgrey", "#E64B35FF"),
                       labels = c("down", "no", "up")) +
    geom_vline(xintercept = c(fc_cutoff,-fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(FDR)") +
    geom_label_repel(aes(label = ifelse(rownames(data) %in% highlight_genes, rownames(data), "")),
                     show.legend = FALSE) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 9)) +
    scale_x_continuous(limits = c(-max, max))
  return(p)
}

Dotplot_class = function(DEtable, title){
  max <- max(abs(min(DEtable$logFC, na.rm = TRUE)), abs(max(DEtable$logFC, na.rm = TRUE)))
  p = ggplot(DEtable, aes(x=`Sub.Class`, y=logFC, color=`Sub.Class`, size=-log10(adj.P.Val))) + 
    geom_point() +
    coord_flip() +
    guides(color=FALSE) +
    scale_color_manual(values=SubclassCol) +
    theme_bw() + 
    theme(axis.text.y = element_text(face = "bold")) +
    ggtitle(title) +
    geom_hline(yintercept=0, linetype="dotted", color="grey") +
    geom_hline(yintercept=c(-1,1), linetype="dotted", color="orange")+
    scale_y_continuous(limits = c(-max, max))
  p
  return(p)
}

ScalebyColumn <- function(x, na.rm=TRUE){
  rm <- colMeans(x, na.rm = na.rm)
  x <- sweep(x, 2, rm)
  sx <- apply(x, 2, sd, na.rm = na.rm)
  x <- sweep(x, 2, sx, "/")
  return(x)
}

VolcanoScatter_exact = function(DEtable, title, GL1, GL2, lFC=0, FDR=0.1, nlabel1=10, nlabel2=10, color_fill1="#6C61AF", color_fill2="#E7B800") {
  DEtable$Lipid = rownames(DEtable)
  DEtable = DEtable[order(DEtable$adj.P.Val), ]
  max <- max(abs(min(DEtable$logFC, na.rm = TRUE)), abs(max(DEtable$logFC, na.rm = TRUE)))
  p = ggplot(DEtable, aes(y = -log10(adj.P.Val), x = logFC, label = Lipid)) +
    geom_point(color = "gray",size=1) + 
    geom_point(data = DEtable %>% filter(Lipid %in% GL1), col = color_fill1, stroke = 0.5, size=2.5,alpha=0.5) + 
    geom_point(data = DEtable %>% filter(Lipid %in% GL2), col = color_fill2, stroke = 0.5, size=2.5,alpha=0.5) +  
    theme_bw() +
    geom_hline(yintercept = -log10(FDR), linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = lFC, linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = -lFC, linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.5) +
    geom_text_repel(data = DEtable[DEtable$Lipid %in% GL1, ] %>% slice_head(n=nlabel1), 
                    max.overlaps = Inf,nudge_x = 0.5, nudge_y = 0.2, color = color_fill1,fontface = "bold") +
    geom_text_repel(data = DEtable[DEtable$Lipid %in% GL2, ] %>% slice_head(n=nlabel2), 
                    max.overlaps = Inf,nudge_x = -0.5, nudge_y = 0.2, color = color_fill2,fontface = "bold") +
    ggtitle(title)+
    scale_x_continuous(limits = c(-max, max))
  
  return(p)
}

#############################
library(scales)
library(ggsci)
library(dplyr)
library(gridExtra)
library(RColorBrewer)
library(corrplot)
library(ggplot2)
library(readxl)
library(limma)
library(reshape2)
library(pheatmap)
library(ggrepel)
library(tidyverse)
library(ggfortify)
library(ggsci)

pal = c("#6C61AF","#4E9040","#3D59A5","#F5931F","#ED3523")
CondCol = pal[c(1,2,3,4,5)]
names(CondCol) = c("QC", "CTL", "GDM", "PE", "PG")

TrendCol = c("#2E9FDF","#E64B35FF")
names(TrendCol) = c("down", "up") 

CategoryCol = brewer.pal(n = 7, "Paired")
names(CategoryCol) = c("Glycerophospholipids", "Glycerolipids", "Saccharolipids", 
                       "Sphingolipids", "PrenolLipids", "SterolLipids", "FattyAcyls")

FunCol <- brewer.pal(3,"Set2")
names(FunCol) = c("Acid-Base Balance", "Oxygen Transport & Gas Exchange", "Osmotic Balance & Circulatory")

#############################
# read data
data = read.csv("Lipidomics_Preprocessed_data.csv", row.names = 1)
lipid_anno = data[,1:23]
Exp = data[,24:85]
rownames(Exp) = data$Name

meta = read_xlsx("Sample_Informationxlsx.xlsx")
id = c("CTL1","CTL2","CTL3","CTL4","CTL5","CTL6","CTL7","CTL8","CTL9","CTL10",
     "GDM1","GDM2","GDM3","GDM4","GDM5","GDM6","GDM7","GDM8","GDM9","GDM10","GDM11","GDM12",
     "PE1","PE2","PE3","PE4","PE5","PE6","PE7","PE8","PE9","PE10",
     "PG1","PG2","PG3","PG4","PG5","PG6","PG7","PG8","PG9","PG10")
meta_order <- meta[match(id, meta$Lipid_name),]
meta_order$Disease_Type = factor(meta_order$Disease_Type, levels=c("QC", "CTL", "GDM", "PE", "PG"))

SubclassCol = colorRampPalette(CategoryCol)(36)
names(SubclassCol) = unique(lipid_anno$Sub.Class)

#############################
## baseline description
# normalize
logExp = log2(Exp)
dim(logExp) #763  62

png("Density_intensity.png", width=5, height=5, units = "in", res = 300)
limma::plotDensities(logExp, legend = F)
dev.off()

longFormatExp = melt(t(logExp))
colnames(longFormatExp) = c("ID", "Lipids","logIntensity")
longFormatExp$Group = gsub("(\\D*)(\\d+)", "\\1", longFormatExp$ID)
longFormatExp$Group = factor(longFormatExp$Group, levels=names(CondCol))
ggplot(longFormatExp,aes(x=ID, y=logIntensity, fill=Group))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust= 1)) +
  scale_fill_manual(values=CondCol)

# after QC
logExp_noQC = logExp[,1:42]
pheatmap(cor(logExp_noQC),clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",clustering_method = "ward.D2",file="Correlation_noQC.png",width=7, height=6)

longFormatExp = melt(t(logExp_noQC))
colnames(longFormatExp) = c("ID", "Lipids","logIntensity")
longFormatExp$Group = gsub("(\\D*)(\\d+)", "\\1", longFormatExp$ID)
longFormatExp$Group = factor(longFormatExp$Group, levels=names(CondCol))
ggplot(longFormatExp,aes(x=ID, y=logIntensity, fill=Group))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust= 1)) +
  scale_fill_manual(values=CondCol)
ggsave("Boxplot_Intensity_noQC.png", width=8, height=5)

# PCA
pca1 = logExp_noQC %>%
  t() %>%
  prcomp(scale. = TRUE, center = TRUE) %>%
  autoplot(
    data = meta_order,
    x=1,
    y=2,
    fill = "Disease_Type",
    frame.colour = "Disease_Type",
    shape = 21,
    size = 3,
    alpha = 0.8,
    frame = TRUE,
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_fill_manual(values=CondCol) +
  scale_color_manual(values=CondCol) +
  geom_text_repel(aes(label=Lipid_name))
pca1
ggsave("PCA_noQC_norm.png", width=6, height=5)

## overall category
lipid_anno_plot <- lipid_anno %>%
  dplyr::count(Categories) %>%
  mutate(Percentage = n / sum(n) * 100)
lipid_anno_plot$CategoryLabel <- paste0(round(lipid_anno_plot$Percentage, 2), "% ", lipid_anno_plot$Categories)
match(names(CategoryCol),lipid_anno_plot$Categories)
lipid_anno_plot <- lipid_anno_plot[c(3,2,5,6,4,7,1),]

color_list = brewer.pal(n = 7, "Paired")
names(color_list) = lipid_anno_plot$CategoryLabel

P1 <- ggplot(lipid_anno_plot, aes(x = 2, y = Percentage, fill = CategoryLabel)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  
  theme_void() +  
  theme(legend.position = "right") +  
  theme(plot.title = element_text(hjust = 0.5)) +  
  scale_fill_manual(values = color_list) +  
  labs(title = NULL) + 
  annotate("text", x = 1.5, y = 0, label = paste("Total =", sum(lipid_anno_plot$n)), size =5 , hjust = 0.5, vjust = 4.5)
P1
ggsave("Pie_Category_PE.png", plot = P1, width=5, height=3)

## overall subclass
data_percentage <- lipid_anno %>%
  group_by(Categories, Sub.Class) %>%
  summarise(Count = n(), .groups = "drop") %>% 
  group_by(Categories) %>%
  mutate(Total = sum(Count), Percentage = (Count / Total) * 100) %>%    
  select(Categories, Sub.Class, Percentage)          

ggplot(data_percentage, aes(x = Categories, y = Percentage, fill = Sub.Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = SubclassCol) +
  labs(
    title = NULL,
    x = NULL,
    y = "Percentage (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    legend.position = "none" 
  )+
  geom_text(aes(label = ifelse(Percentage > 5, Sub.Class, "")), 
            position = position_stack(vjust = 0.5), 
            size = 3, 
            color = "black")
ggsave("Subclass_all_in.png", width=5, height=4.5)

#############################
#### differential analysis
library(ggpubr)
meta_order$Disease_Type = factor(meta_order$Disease_Type, levels=c("CTL", "GDM", "PE", "PG"))
design = model.matrix(~Disease_Type,data = meta_order)
fit1 <- lmFit(logExp_noQC, design)
fit1 <- eBayes(fit1)
nLipid = nrow(logExp_noQC)

GDMvsCTL = topTable(fit1, coef = 2, number = nLipid)
print(sum(GDMvsCTL$adj.P.Val < 0.25)) #20/12
top_GDM = GDMvsCTL[GDMvsCTL$adj.P.Val < 0.25,]
PEvsCTL = topTable(fit1, coef = 3, number = nLipid)
print(sum(PEvsCTL$adj.P.Val < 0.25)) #128/26
top_PE = PEvsCTL[PEvsCTL$adj.P.Val < 0.25,]
PGvsCTL = topTable(fit1, coef = 4, number = nLipid)
print(sum(PGvsCTL$adj.P.Val < 0.25)) #89/19
top_PG = PGvsCTL[PGvsCTL$adj.P.Val < 0.25,]

write.csv(GDMvsCTL,"GDMvsCTL_All_Lipids.csv")
write.csv(PEvsCTL,"PEvsCTL_All_Lipids.csv")
write.csv(PGvsCTL,"PGvsCTL_All_Lipids.csv")

V1 <- Volcano_plot(GDMvsCTL, "GDM vs CTL",pval_cutoff = 0.25, fc_cutoff = log2(1),nlabel = 0)
V1
ggsave("Volcano_GDMvsCTL.png", plot = V1, width=4, height=4)
V2 <- Volcano_plot(PEvsCTL, "PE vs CTL",pval_cutoff = 0.25, fc_cutoff = log2(1),nlabel = 0)
V2
ggsave("Volcano_PEvsCTL.png", plot = V2, width=4, height=4)
V3 <- Volcano_plot(PGvsCTL, "PG vs CTL",pval_cutoff = 0.25, fc_cutoff = log2(1),nlabel = 0)
V3
ggsave("Volcano_PGvsCTL.png", plot = V3, width=4, height=4)

top_GDM = rownames_to_column(top_GDM,var="Name")
GDMDL = top_GDM %>% left_join(lipid_anno[,c("Main.Class", "Sub.Class", "Name")], by="Name")
D1 <- Dotplot_class(GDMDL, title = "GDM vs CTL")
D1
ggsave("Dotplot_GDMvsCTL.png",D1, width=5, height=2)

top_PE = rownames_to_column(top_PE,var="Name")
PEDL = top_PE %>% left_join(lipid_anno[,c("Main.Class", "Sub.Class", "Name")], by="Name")
D2 <- Dotplot_class(PEDL, title = "PE vs CTL")
D2
ggsave("Dotplot_PEvsCTL.png",D2, width=5, height=5)

top_PG = rownames_to_column(top_PG,var="Name")
PGDL = top_PG %>% left_join(lipid_anno[,c("Main.Class", "Sub.Class", "Name")], by="Name")
D3 <- Dotplot_class(PGDL, title = "PG vs CTL")
D3
ggsave("Dotplot_PGvsCTL.png",D3, width=5, height=4)

## heatmap
top_GDM <- top_GDM %>% filter(adj.P.Val < 0.25&abs(logFC) > log2(2))
rownames(top_GDM) <- top_GDM$Name
dim(top_GDM)
top_PE <- top_PE %>% filter(adj.P.Val < 0.25&abs(logFC) > log2(2))
rownames(top_PE) <- top_PE$Name
dim(top_PE)
top_PG <- top_PG %>% filter(adj.P.Val < 0.25&abs(logFC) > log2(2))
rownames(top_PG) <- top_PG$Name
dim(top_PG)

top_PE$Trend = ifelse(top_PE$logFC > 0, "up", "down")
top_GDM$Trend = ifelse(top_GDM$logFC > 0, "up", "down")
top_PG$Trend = ifelse(top_PG$logFC > 0, "up", "down")

dim(LP_DE)
LP_DE$CTL_median <- apply(LP_DE[,meta_order$Disease_Type == "CTL"], 1, median)
LP_DE$GDM_median <- apply(LP_DE[,meta_order$Disease_Type == "GDM"], 1, median)
LP_DE$PE_median <- apply(LP_DE[,meta_order$Disease_Type == "PE"], 1, median)
LP_DE$PG_median <- apply(LP_DE[,meta_order$Disease_Type == "PG"], 1, median)

LP_SCALE = ScalebyRow(LP_DE[,c(43:46)])
colnames(LP_SCALE) <- c("CTL", "GDM", "PE", "PG")
LP_SCALE[LP_SCALE > 3] = 3
LP_SCALE[LP_SCALE < -3] = -3

annotation_row <- data.frame(row.names = rownames(LP_DE),Category = factor(LP_DEINFO$Categories, levels=unique(LP_DEINFO$Categories)),
                             #Subclass = factor(LP_DEINFO$Sub.Class, levels=unique(LP_DEINFO$Sub.Class)),
                             GDM = factor(ifelse(rownames(LP_DE) %in% rownames(top_GDM),top_GDM[rownames(LP_DE),]$Trend,NA)),
                             PE = factor(ifelse(rownames(LP_DE) %in% rownames(top_PE),top_PE[rownames(LP_DE),]$Trend,NA)),
                             PG = factor(ifelse(rownames(LP_DE) %in% rownames(top_PG),top_PG[rownames(LP_DE),]$Trend,NA)))
annotation_col <- data.frame(row.names = colnames(LP_SCALE),Disease_Type = factor(c("CTL", "GDM", "PE", "PG"), levels=c("CTL", "PE", "GDM", "PG")))
pheatmap(LP_SCALE,  
         cluster_cols = F, 
         cluster_rows = T,
         annotation_row=annotation_row, 
         annotation_col = annotation_col,
         annotation_colors = list(Disease_Type = CondCol[2:5],Category = CategoryCol[unique(LP_DEINFO$Categories)],
                                  #Subclass = SubclassCol[unique(LP_DEINFO$Sub.Class)],
                                  GDM = TrendCol, PE = TrendCol, PG = TrendCol),
         clustering_distance_rows = "correlation", 
         filename = "Heatmap_DELipid_median.png", width=6.5, height=8,
         clustering_method = "ward.D2",
         cutree_rows = 6)

################################
## PE predictor
library(venn)
png("Venn_PE_PG.png", width=4, height=4, units = "in", res = 300)
venn(list(PE=top_PE$Name,PG=top_PG$Name),counts = T,
     col = "transparent",
     ilabels="counts",
     zcolor = c("#F5931F","#ED3523"),
     ilcs = 2, sncs = 1.5)
dev.off()

library(glmnet)
PE_meta = meta_order %>% filter(Disease_Type == "PE"|Disease_Type == "CTL"|Disease_Type == "PG")
PE_meta$Disease_Type <- factor(PE_meta$Disease_Type, levels = c("CTL","PE","PG"))
PE_exp = logExp_noQC[,PE_meta$Lipid_name]
PE_DE <- PE_exp[intersect(top_PE$Name,top_PG$Name),]
PE_DE = as.data.frame(PE_DE)
PE_DE_t = as.data.frame(t(PE_DE))
PE_DE_t <- rownames_to_column(PE_DE_t, var = "Lipid_name")
PE_DE_t$Disease_Type = PE_meta$Disease_Type[match(PE_DE_t$Lipid_name, PE_meta$Lipid_name)]
PE_DE_t$Disease_Type <- factor(PE_DE_t$Disease_Type, levels = c("CTL","PE","PG"))

cor_matrix <- cor(t(PE_DE[intersect(top_PE$Name,top_PG$Name), ]), method = "pearson")
hc <- hclust(as.dist(1 - abs(cor_matrix)), method = "ward.D2")
png(filename = "Lipids_cor(PEPG).png",width = 7, height = 5, units = "in", res = 300)
corrplot(cor_matrix[hc$order, hc$order],type="upper",tl.pos = "tp", tl.col = "black")
corrplot(cor_matrix[hc$order, hc$order],add=TRUE, type="lower", method="number",diag=FALSE,tl.pos="n", cl.pos="n")
dev.off()

## predict with overlap lipids -- 8
library(pROC)
colors <- brewer.pal(n = 8, "Paired")
auc_table <- data.frame(Lipid = character(),AUC = numeric(),CI_lower = numeric(),CI_upper = numeric(),stringsAsFactors = FALSE)
png("ROC_lipids.png", width=7, height=7, units = "in", res = 300)
plot(1, type = "n", xlab = "Specificity", ylab = "Sensitivity", 
     xlim = c(1, 0), ylim = c(0, 1), main = NULL)
abline(a = 1, b =-1, col = "gray", lwd=3, lty = 1)

lipids <- intersect(top_PE$Name,top_PG$Name)
for (i in seq_along(lipids)) {
  roc_result <- roc(PE_DE_t$Disease_Type, PE_DE_t[[lipids[i]]])
  auc_value <- auc(roc_result)
  ci_values <- ci.auc(roc_result)
  lines(roc_result, col = colors[i], lwd = 3, lty = 1)
  
  auc_table <- rbind(auc_table, data.frame(Lipid = lipids[i],
                                                       AUC = round(auc_value, 3),
                                                       CI_lower = round(ci_values[1], 3),
                                                       CI_upper = round(ci_values[3], 3)))
}
legend("bottomright", legend = lipids, col = colors[seq_along(lipids)], lwd = 3, lty = 1)
dev.off()
auc_table
write.csv(auc_table, "ROC_lipids.csv", row.names = FALSE)

## calculate PLRS -- Ridge regression score
PE_lipids_scale <- as.matrix(ScalebyColumn(PE_DE_t[,lipids]))
rownames(PE_lipids_scale) <- rownames(PE_DE_t)

Y <- factor(PE_DE_t$Disease_Type, levels = c("CTL","PE","PG"))
ridge_model <- cv.glmnet(PE_lipids_scale, Y, alpha = 0, family = "multinomial",nfolds = 10)
best_lambda_ridge <- ridge_model$lambda.min
ridge_best <- glmnet(PE_lipids_scale, Y, alpha = 0, family = "multinomial", lambda = best_lambda_ridge)
ridge_coefficients <- coef(ridge_best)

PE_lipids <- as.data.frame(t(logExp_noQC[rownames(logExp_noQC) %in% lipids,]))
PE_lipids <- as.data.frame(ScalebyColumn(PE_lipids))
PE_lipids <- tibble::rownames_to_column(PE_lipids, var = "Lipid_name")
PE_lipids$PE <- c(rep("0",22),rep("1",20))
PE_lipids$Disease_Type <- meta_order$Disease_Type[match(PE_lipids$Lipid_name, meta_order$Lipid_name)]
dim(PE_lipids)

weights <- (ridge_coefficients$PE[-1]+ridge_coefficients$PG[-1])/2
# [1] -0.26254522  0.15614200  0.16307295  0.07829566  0.35224435  0.05596146 -0.13137539
# [8] -0.11137356
lipids <- c("LPE(18:2)","SM(d36:0)(rep)","SM(d36:0)","SM(d38:3)(rep)",
            "PC(18:0/22:6)(rep)","SM(d45:5)", "LPC(14:0)(rep)","LPC(18:4)")

PLRS <- as.matrix(PE_lipids[,lipids]) %*% as.numeric(weights)
PLRS <- as.vector(PLRS)
PE_lipids$PLRS <- PLRS
# c("LPE(18:2)","SM(d36:0)(rep)","SM(d36:0)","SM(d38:3)(rep)",
#   "PC(18:0/22:6)(rep)","SM(d45:5)", "LPC(14:0)(rep)","LPC(18:4)")
write.csv(PE_lipids[,c("Lipid_name","PLRS")], "PLRS_PE.csv", row.names = FALSE)

B2<-ggboxplot(PE_lipids,x="Disease_Type",y="PLRS",
              color = "Disease_Type",palette = "jco",add = "jitter",
              short.panel.labs = FALSE)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank())+
  ylab("PLRS")+
  scale_color_manual(values=CondCol[2:5])+
  theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),legend.position = "none")+xlab(NULL)
B2
my_comparisons<-list(c("CTL","GDM"),c("CTL","PE"),c("CTL","PG"),c("GDM","PG"),c("PE","PG"))
B2<-B2+stat_compare_means(label = "p.signif",method="wilcox.test",comparisons = my_comparisons)
B2
ggsave("PLRS_PE_allgroup.png",plot = B2, width=4, height=4)

## correlation analysis with the PLRS-Baby_weight
library(ggpubr)
PE_meta$PLRS <- PE_lipids[PE_lipids$Lipid_name %in% PE_meta$Lipid_name,"PLRS"]
colnames(PE_meta)
CR2 <- ggplot(PE_meta,aes(x=PLRS, y=SBP)) +
  geom_point(aes(color = Disease_Type),size=2) + 
  scale_color_manual(values=CondCol[c(2,4,5)])+
  geom_smooth(method = "lm", se = TRUE,  color = "skyblue", fill = "gray") + 
  stat_cor(method = "pearson") + 
  theme_bw() + 
  labs(title = NULL,
       x = "PLRS", 
       y = "SBP(mmHg)")
ggsave("Correlation_PLRS_SBP.png", plot = CR2, width=4.5, height=3.5)

CR3 <- ggplot(PE_meta,aes(x=PLRS, y=DBP)) +
  geom_point(aes(color = Disease_Type),size=2) + 
  scale_color_manual(values=CondCol[c(2,4,5)])+
  geom_smooth(method = "lm", se = TRUE,  color = "skyblue", fill = "gray") + 
  stat_cor(method = "pearson") + 
  theme_bw() + 
  labs(title = NULL,
       x = "PLRS", 
       y = "DBP(mmHg)")
ggsave("Correlation_PLRS_DBP.png", plot = CR3, width=4.5, height=3.5)

#############################
### Correlation of Lipids and blood-gas
colnames(meta_order)
Lipids <- rownames(logExp_noQC)
BLGAS <- c("HCO3(std)","BE(ecf)","pH","pO2","sO2(est)","mOsm","RI","tHb(est)")

bg <- meta_order[!is.na(meta_order$pH),c(2,61:82)]
logExp_noQC_t <- t(logExp_noQC)
logExp_noQC_t <- data.frame(Lipid_name = rownames(logExp_noQC_t), logExp_noQC_t)
colnames(logExp_noQC_t) <- c("Lipid_name", rownames(logExp_noQC))
logExp_noQC_t$Disease_Type <- meta_order$Disease_Type
bg <- left_join(bg, logExp_noQC_t, by = "Lipid_name")
dim(bg) #20 787

correlation_LB_r <- matrix(NA, nrow = length(Lipids), ncol = length(BLGAS))
correlation_LB_p <- matrix(NA, nrow = length(Lipids), ncol = length(BLGAS))
colnames(correlation_LB_r) <- BLGAS
rownames(correlation_LB_r) <- Lipids
colnames(correlation_LB_p) <- BLGAS
rownames(correlation_LB_p) <- Lipids

for (i in seq_along(Lipids)) {
  for (j in seq_along(BLGAS)) {
    test <- cor.test(bg[[Lipids[i]]], bg[[BLGAS[j]]])
    if (test$p.value < 0.05) {
      correlation_LB_r[i, j] <- test$estimate  
    } else {
      correlation_LB_r[i, j] <- NA  
    }
    correlation_LB_p[i, j] <- test$p.value  
  }
}
correlation_LB_r <- as.data.frame(correlation_LB_r)
correlation_LB_p <- as.data.frame(correlation_LB_p)

Oxygen <- correlation_LB_r[!is.na(correlation_LB_r$`sO2(est)`)&
                             !is.na(correlation_LB_r$pO2)&
                             !is.na(correlation_LB_r$`tHb(est)`),
                           c("sO2(est)","pO2","tHb(est)")]
Oxygen
Oxygen <- Oxygen[order(Oxygen$`sO2(est)`),]
Oxygen$lipid <- factor(rownames(Oxygen), levels = rownames(Oxygen))
Oxygen_log <- melt(Oxygen)
Oxygen_log$pvalue <- melt(rownames_to_column(correlation_LB_p[rownames(Oxygen), c("sO2(est)","pO2","tHb(est)")], var = "lipid"))$value
ggplot(Oxygen_log, aes(x = variable, y = lipid, color = value, size = pvalue)) +
  geom_point(alpha = 1) +
  scale_size(range = c(3, 8)) +
  scale_color_gradient(low="#2E9FDF",high="#E64B35FF",name = "correlation") +
  theme_bw() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(#legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10))
ggsave("Oxygen Transport & Gas Exchange.png", width=4, height=4)

AcidBase <- correlation_LB_r[!is.na(correlation_LB_r$pH)&
                               !is.na(correlation_LB_r$`HCO3(std)`)&
                               !is.na(correlation_LB_r$`BE(ecf)`),
                             c("pH","HCO3(std)","BE(ecf)")]
AcidBase
AcidBase <- AcidBase[order(AcidBase$pH),]
AcidBase$lipid <- factor(rownames(AcidBase), levels = rownames(AcidBase))
AcidBase_log <- melt(AcidBase)
AcidBase_log$pvalue <- melt(rownames_to_column(correlation_LB_p[rownames(AcidBase), c("pH","HCO3(std)","BE(ecf)")], var = "lipid"))$value
ggplot(AcidBase_log, aes(x = variable, y = lipid, color = value, size = pvalue)) +
  geom_point(alpha = 1) +
  scale_size(range = c(3, 8)) +
  scale_color_gradient(low="#2E9FDF",high="#E64B35FF",name = "correlation") +
  theme_bw() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
ggsave("Acid-Base Balance.png", width=3, height=2.2)
Osmotic <- correlation_LB_r[!is.na(correlation_LB_r$mOsm)&
                              !is.na(correlation_LB_r$RI),
                            c("mOsm", "RI"), drop = F]
Osmotic
Osmotic <- Osmotic[order(Osmotic$mOsm),]
Osmotic$lipid <- factor(rownames(Osmotic), levels = rownames(Osmotic))
Osmotic_log <- melt(Osmotic)
Osmotic_log$pvalue <- melt(rownames_to_column(correlation_LB_p[rownames(Osmotic), c("mOsm","RI")], var = "lipid"))$value
ggplot(Osmotic_log, aes(x = variable, y = lipid, color = value, size = pvalue)) +
  geom_point(alpha = 1) +
  scale_size(range = c(3, 8)) +
  scale_color_gradient(low="#2E9FDF",high="#E64B35FF",name = "correlation") +
  theme_bw() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))
ggsave("Osmotic Balance & Circulatory.png", width=2.5, height=1.5)

list1 <- unique(c(Osmotic$lipid, AcidBase$lipid, Oxygen$lipid))
# list1 <- c("PC(16:0e/22:6)", "PC(40:6e)(rep)",
#   "SM(d43:1)(rep)", "PC(32:0e)", "PC(38:3)(rep)", "PC(20:0p/20:4)",
#   "PC(18:1/24:0)", "LPI(18:0)", "PC(30:0)(rep)", "PC(16:0/16:1)",
#   "PE(16:0/18:1)", "PI(16:0/16:1)", "PC(32:1)(rep)", "PC(40:4)",
#   "PC(28:0)", "dMePE(16:0/16:1)", "dMePE(18:2/18:2)")

for (i in seq_along(Lipids)) {
  for (j in seq_along(BLGAS)) {
    test <- cor.test(bg[[Lipids[i]]], bg[[BLGAS[j]]])
    correlation_LB_r[i, j] <- test$estimate
    correlation_LB_p[i, j] <- test$p.value
  }
}
heatmap_data <- correlation_LB_r[rownames(correlation_LB_r) %in% list1,c("HCO3(std)", "BE(ecf)", "pH", "pO2", "sO2(est)","tHb(est)","mOsm", "RI")]
pheatmap_res <- pheatmap(heatmap_data, 
         cluster_cols = T, 
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T,
         cutree_rows = 2,
         cutree_cols = 2)
row_clusters <- cutree(pheatmap_res$tree_row, k = 2)
annotation_row <- data.frame(row.names = rownames(heatmap_data),
                             Category = lipid_anno[lipid_anno$Name %in% rownames(heatmap_data), "Categories"],
                             Subclass = lipid_anno[lipid_anno$Name %in% rownames(heatmap_data), "Sub.Class"],
                             Cluster = as.factor(row_clusters))
annotation_col <- data.frame(
  row.names = c("pH","HCO3(std)","BE(ecf)","sO2(est)","pO2","tHb(est)","mOsm","RI"),
  Outcome = c(rep("down-worse", 5), rep("up-worse", 3)),
  Function = c("Acid-Base Balance",               # pH
               "Acid-Base Balance",               # HCO3(std)
               "Acid-Base Balance",               # BE(ecf)
               "Oxygen Transport & Gas Exchange",  # sO2(est)
               "Oxygen Transport & Gas Exchange",  # pO2
               "Oxygen Transport & Gas Exchange", # tHb(est)
               "Osmotic Balance & Circulatory",       # mOsm
               "Osmotic Balance & Circulatory")       # RI
)
annotation_col <- annotation_col[match(colnames(heatmap_data), rownames(annotation_col)), ]
pheatmap(heatmap_data, 
         cluster_cols = T, 
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T,
         display_numbers = T,
         cutree_rows = 2,
         cutree_cols = 2,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = list(Category=CategoryCol[unique(annotation_row$Category)], 
                                  Subclass=SubclassCol[unique(annotation_row$Subclass)],
                                  Outcome=c("up-worse"="#E64B35FF","down-worse"="#2E9FDF"),
                                  Function=FunCol,
                                  Cluster=c("1"="#E7B800","2"="#6C61AF")),
         filename = "Heatmap_outcome_lipids.png", 
         width=10, height=6)
write.csv(heatmap_data,"Outcome_lipid.csv",row.names = T)

#########################################
## volcano plot of Outcome-related Lipids
GL1 <- c("PI(16:0/16:1)","PC(16:0/16:1)","PE(16:0/18:1)","dMePE(16:0/16:1)","PC(20:0p/20:4)",
         "PC(28:0)","PC(38:3)(rep)","PC(30:0)(rep)","PC(32:1)(rep)","PC(32:0e)","PC(40:4)","SM(d43:1)(rep)")
GL2 <- c("LPI(18:0)","PC(18:1/24:0)","dMePE(18:2/18:2)","PC(16:0e/22:6)","PC(40:6e)(rep)")
VC1 <- VolcanoScatter_exact(GDMvsCTL,title = "GDM vs CTL",GL1,GL2,lFC = log2(2),FDR=0.25,nlabel1=0,nlabel2=2)
VC2 <- VolcanoScatter_exact(PEvsCTL,title = "PE vs CTL",GL1,GL2,lFC = log2(2),FDR=0.25,nlabel1=5,nlabel2=2)
VC3 <- VolcanoScatter_exact(PGvsCTL,title = "PG vs CTL",GL1,GL2,lFC = log2(2),FDR=0.25,nlabel1=1,nlabel2=2)
ggsave("Volcano_GDMvsCTL_outcome.png",VC1,width=3,height=3)
ggsave("Volcano_PEvsCTL_outcome.png",VC2,width=3,height=3)
ggsave("Volcano_PGvsCTL._outcome.png",VC3,width=3,height=3)

## violin plot of Outcome-related Lipids
GL1 <- c("PC(40:4)","PC(16:0/16:1)","PI(16:0/16:1)","PC(32:1)(rep)","dMePE(16:0/16:1)")
GL2 <- c("PC(18:1/24:0)","dMePE(18:2/18:2)")
plots_violin <- list()
for (i in seq_along(GL1)) {
  plot <- ggplot(logExp_noQC_t, aes(x = Disease_Type, y = .data[[GL1[i]]], fill = Disease_Type)) +
    geom_violin(trim = FALSE, width = 0.8, scale="width") +  
    geom_boxplot(width = 0.3, position = position_dodge(0.9),fill="white") +  
    scale_fill_manual(values = CondCol) +  
    theme_bw() +  
    labs(y = "Relatvie Expression", x = "", title = GL1[i]) +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +
    stat_compare_means(comparisons = list(c("GDM", "CTL"), c("PE","CTL"),c("PG","CTL")),method = "t.test", label = "p.signif", size = 3)
  plots_violin[[i]] <- plot
  ggsave(paste0("violin_up_(",i,").png"),plot, width=3, height=3)
}
for (i in seq_along(GL2)) {
  plot <- ggplot(logExp_noQC_t, aes(x = Disease_Type, y = .data[[GL2[i]]], fill = Disease_Type)) +
    geom_violin(trim = FALSE, width = 0.8, scale="width") +  
    geom_boxplot(width = 0.3, position = position_dodge(0.9),fill="white") +  
    scale_fill_manual(values = CondCol) +  
    theme_bw() +  
    labs(y = "Relatvie Expression", x = "", title = GL2[i]) +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +
    stat_compare_means(comparisons = list(c("GDM", "CTL"), c("PE","CTL"),c("PG","CTL")),method = "t.test", label = "p.signif", size = 3)
  plots_violin[[i]] <- plot
  ggsave(paste0("violin_down_(",i,").png"),plot, width=3, height=3)
}

###########################
## venn and exp
library(venn)
png("venn.png", width=5, height=5, units = 'in', res = 300)
venn(list(Disease_related = rownames(LP_DE), Outcome_related = list1),
     ilab = c("Disease_related", "Outcome_related"),
     col = "transparent",
     ilabels = "counts",
     zcolor = c("#E7B800","#2E9FDF"),
     ilcs = 2, sncs = 1.5)
dev.off()
intersect(list1,rownames(LP_DE)) #"PC(18:1/24:0)"

v_plot <- logExp_noQC_t[,c("Lipid_name","Disease_Type","PC(18:1/24:0)")]
colnames(v_plot) <- c("Lipid_name","Disease_Type","Value")
VL1 <- ggplot(v_plot, aes(x = Disease_Type, y = Value, fill = Disease_Type)) +
  geom_violin(trim = FALSE, width = 0.8,scale = "width") +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = CondCol[2:5]) +
  labs(title = "PC(18:1/24:0)", y = "Relative Expression", x = "") +
  theme_bw()+
  theme(legend.position = "none")
VL1
ggsave("PC.png", width=3.5, height=3.5)

lasso_features <- c("PC(18:2/18:2)", "PC(36:4)(rep)(rep)", "SM(d42:4)", "SM(d44:3)(rep)", 
                    "LPC(18:2)(rep)", "DGDG(30:3)", "SM(d42:4)", "SM(d38:4)", "DGDG(30:3)", 
                    "PC(35:2p)", "LPC(18:2)(rep)", "SM(d42:6)", "SM(d38:5)", "PC(18:0/22:6)(rep)", 
                    "SM(d42:1)(rep)", "SM(d45:5)")
list1 <- c("PC(16:0e/22:6)", "PC(40:6e)(rep)","SM(d43:1)(rep)", "PC(32:0e)", "PC(38:3)(rep)", "PC(20:0p/20:4)", 
          "PC(18:1/24:0)", "LPI(18:0)", "PC(30:0)(rep)", "PC(16:0/16:1)","PE(16:0/18:1)", "PI(16:0/16:1)", 
          "PC(32:1)(rep)", "PC(40:4)","PC(28:0)", "dMePE(16:0/16:1)", "dMePE(18:2/18:2)")
plots_violin <- list()
for (i in seq_along(list1)) {
  plot <- ggplot(logExp_noQC_t, aes(x = Disease_Type, y = .data[[list1[i]]], fill = Disease_Type)) +
    geom_violin(trim = FALSE, width = 0.8) +  
    geom_boxplot(width = 0.2, position = position_dodge(0.9)) +  
    scale_fill_manual(values = CondCol) +  
    theme_bw() +  
    labs(y = list1[i], x = "")+
    stat_compare_means(comparisons = list(c("GDM", "CTL"), c("PE","CTL"),c("PG","CTL")),method = "t.test", label = "p.signif", size = 3)
  plots_violin[[i]] <- plot
}
VL0 <- do.call(ggarrange, c(plots_violin, list(ncol = 5, nrow = 4, common.legend = TRUE, legend = "right")))
VL0
ggsave("violin_of_outcome_lipids_t.png",VL0, width=15, height=12)
