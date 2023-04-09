rm(list = ls())
load("./data analysis/mouse olfactory bulb/MOB.RData")
load("./data analysis/mouse olfactory bulb/MOB_results.RData")

##----------------------------------------------------------
## Moran's I
##----------------------------------------------------------
source("./funcs/ValidationFuncs.R")
library(spdep)

genes.overlap <- intersect(genes_rep_bh, genes_rep_maxp)
genes.bh.only  <- genes_rep_bh[!genes_rep_bh%in%genes.overlap]
genes.jump.only <- genes_rep_jump[!genes_rep_jump%in%genes.overlap]

## MI based on ST Rep9 data
stat.res1 <- corrValue(counts1, location1)
MI1 <- stat.res1
mi_all1 <- MI1[overlap]
mi_jump_only1 <- MI1[genes.jump.only]

# MI of additional replicable SVGs identified by JUMP
mean(mi_all1); mean(mi_jump_only1); median(mi_all1); median(mi_jump_only1)
# 0.0503; 0.1147; 0.0392; 0.1085
factor1 <- factor(rep(c("JUMP only", "All (ST)"), 
                   times = c(length(genes.jump.only), length(overlap))))
dataset1 <- data.frame(value = c(mi_jump_only1, mi_all1), group = factor1)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset1, col = c("antiquewhite1", "lightSalmon"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6) # 4 * 6 inches

## MI based on 10X Visium data
stat.res2 <- corrValue(counts2, location2)
MI2 <- stat.res2
mi_all2 <- MI2[overlap]
mi_jump_only2 <- MI2[genes.jump.only]

mean(mi_all2); mean(mi_jump_only2); median(mi_all2); median(mi_jump_only2)
# 0.1934; 2672; 0.1478; 0.2341
factor2 <- factor(rep(c("JUMP only", "All (10X)"), 
                      times = c(length(genes.jump.only), length(overlap))))
dataset2 <- data.frame(value = c(mi_jump_only2, mi_all2), group = factor2)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset2, col = c("antiquewhite1", "lightSalmon"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)

##--------------------------------------------------------------------
## Summarized spatial expression patterns
##--------------------------------------------------------------------
library(amap)
source("./funcs/PlotFuncs.R")

# Summarized spatial expression pattern based on ST data
vst_count1 <- var_stabilize(counts1) # R function in funcs.R
sig_vst_count1 <- vst_count1[genes_rep_jump, ]
sig_vst_res1 <- t(apply(sig_vst_count1, 1, LMReg, T = log(apply(counts1, 2, sum))))
hc1 <- hcluster(sig_vst_res1, method = "euc", link = "ward", nbproc = 1, 
                doubleprecision = TRUE)
numC <- 3
memb1 <- cutree(hc1, k = numC)
# The mean residuals of the three patterns for each location
cent1 <- NULL
for (k in 1:numC) {
  cent1 <- cbind(cent1, colMeans(sig_vst_res1[memb1 == k, , drop = FALSE]))
}
position_cord1 <- location1
rownames(position_cord1) <- rownames(cent1)
# The relative residuals
rel_cent1 <- t(apply(cent1, 1, relative_func))
# rel_cent1 <- rel_cent1[,c(3,1,2)]
pd1 <- setNames(cbind.data.frame(position_cord1, rel_cent1), c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP1 <- lapply(1:numC, function(x) {
  pattern_plot2(pd1, x, xy = T, main = T, titlesize = 1.7)
})
grid.arrange(grobs = MBP1[numC:1], nrow = 1) # 3*8 inches

# Summarized spatial expression pattern based on 10X Visium data
vst_count2 <- var_stabilize(counts2) # R function in funcs.R
sig_vst_count2 <- vst_count2[genes_rep_jump, ]
sig_vst_res2 <- t(apply(sig_vst_count2, 1, LMReg, T = log(apply(counts2, 2, sum))))
hc2 <- hcluster(sig_vst_res2, method = "euc", link = "ward", nbproc = 1, 
                doubleprecision = TRUE)
numC <- 3
memb2 <- cutree(hc2, k = numC)
# The mean residuals of the three patterns for each location
cent2 <- NULL
for (k in 1:numC) {
  cent2 <- cbind(cent2, colMeans(sig_vst_res2[memb2 == k, , drop = FALSE]))
}
position_cord2 <- location2
rownames(position_cord2) <- rownames(cent2)
# The relative residuals
rel_cent2 <- t(apply(cent2, 1, relative_func))
# rel_cent2 <- rel_cent2[,c(3,1,2)]
pd2 <- setNames(cbind.data.frame(position_cord2, rel_cent2), c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP2 <- lapply(1:numC, function(x) {
  pattern_plot2(pd2, x, xy = T, main = T, titlesize = 3.3)
})
grid.arrange(grobs = MBP2[numC:1], nrow = 1) # 4*15 inches

##--------------------------------------------------------------------
## UMAP clusters
##--------------------------------------------------------------------
library(umap)
library(amap)
source("./funcs/PlotFuncs.R")

# ST data
sig_counts <- counts1[genes.jump.only,]
scaled_counts <- t(scale(t(sig_counts)))
umap_results <- umap::umap(scaled_counts)
labels <- memb1[genes.jump.only]
labels[which(labels == 1)] <- "III"
labels[which(labels == 2)] <- "II"
labels[which(labels == 3)] <- "I"
df <- data.frame(umap_results$layout[,1], umap_results$layout[,2], labels)
colnames(df) <- c("UMAP1","UMAP2", "Cluster")
# 4 * 5 inches
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster), color = labels) + 
  geom_point() + 
  scale_colour_manual(name="Cluster",  
                      values = c("I"="#3cc08f", "II"="#F4A460", "III"="Steelblue"))

# 10X Visium data
sig_counts <- counts2[genes_rep_jump,]
scaled_counts <- t(scale(t(as.matrix(sig_counts))))
umap_results <- umap::umap(scaled_counts)
labels <- memb2
labels[which(labels == 1)] <- "III"
labels[which(labels == 2)] <- "II"
labels[which(labels == 3)] <- "I"
df <- data.frame(umap_results$layout[,1], umap_results$layout[,2], labels)
colnames(df) <- c("UMAP1","UMAP2", "Cluster")
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster), color = labels) + 
  geom_point() + 
  scale_colour_manual(name="Cluster",  
                      values = c("I"="#3cc08f", "II"="#F4A460", "III"="Steelblue"))

##--------------------------------------------------------------------
## Spatial expression patterns of randomly selected SVGs
##--------------------------------------------------------------------
source("./funcs/PlotFuncs.R")

# three representative genes corresponding to the main patterns
gene_plot <- c("Kctd12", "Plcxd2", "Ctxn1") 

randi <- sample(length(genes.jump.only), 30, replace = FALSE)
gene_plot <- genes.jump.only[randi]

# Spatial expression patterns based on ST data
vst_count1 <- var_stabilize(counts1) # R function in funcs.R
sig_vst_ct1 <- vst_count1[gene_plot, ]
rel_vst_ct1 <- apply(sig_vst_ct1, 1, relative_func)
pltdat1 <- cbind.data.frame(location1[,1:2],rel_vst_ct1)
genetitle <- gene_plot
pp1 <- lapply(1:(ncol(pltdat1)-2),
              function(x){pattern_plot2(pltdat1,x,main=T,titlesize=1.7,title=genetitle[x])})
grid.arrange(grobs=pp1, nrow=3) # 7.5 * 25 inches

# Spatial expression patterns based on 10X Visium data
vst_count2 <- var_stabilize(counts2) # R function in funcs.R
sig_vst_ct2 <- vst_count2[gene_plot, ]
rel_vst_ct2 <- apply(sig_vst_ct2, 1, relative_func)
# rel_vst_ct2 <- relative_func(sig_vst_ct2)
pltdat2 <- cbind.data.frame(location2[,1:2],rel_vst_ct2)
genetitle <- gene_plot
pp2 <- lapply(1:(ncol(pltdat2)-2),
              function(x){pattern_plot2(pltdat2,x,main=T,titlesize=3.6,title=genetitle[x])})
grid.arrange(grobs=pp2, nrow=3) # 12 * 50 inches

##--------------------------------------------------------------------
## Validation gene sets
##--------------------------------------------------------------------
## The highlighted marker genes in the original study 
spatial.genes <- c("Doc2g", "Slc17a7", "Reln", "Cdhr1", "Sv2b", "Shisa3", "Plcxd2", "Nmb", "Uchl1", "Rcan2")
length(intersect(genes_rep_maxp, spatial.genes)) # 6/618
length(intersect(genes_rep_bh, spatial.genes)) # 6/638
length(intersect(genes_rep_jump, spatial.genes)) # 7/810

length(intersect(genes.bh.only, spatial.genes)) # 0/20
length(intersect(genes.jump.only, spatial.genes)) # 1/189

## Genes related to the main olfactory bulb in Harmonizome database: Allen Brain Atlas
library(rjson)
library(stringr)
Harmonizome_Glomerular <- fromJSON(paste(readLines("./data analysis/mouse olfactory bulb/validation sets/Harmonizome_Glomerular.json")))
Harmonizome_Granule <- fromJSON(paste(readLines("./data analysis/mouse olfactory bulb/validation sets/Harmonizome_Granule.json")))
Harmonizome_Mitral <- fromJSON(paste(readLines("./data analysis/mouse olfactory bulb/validation sets/Harmonizome_Mitral.json")))

Harmonizome_Glomerular <- Harmonizome_Glomerular[["associations"]]
Harmonizome_Granule <- Harmonizome_Granule[["associations"]]
Harmonizome_Mitral <- Harmonizome_Mitral[["associations"]]

Allen.genes <- NULL
for (i in 1:length(Harmonizome_Glomerular)){
  Allen.genes   <- c(Allen.genes, Harmonizome_Glomerular[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_Granule)){
  Allen.genes   <- c(Allen.genes, Harmonizome_Granule[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_Mitral)){
  Allen.genes   <- c(Allen.genes, Harmonizome_Mitral[[i]][["gene"]][["symbol"]])
}
Allen.genes <- Allen.genes[!duplicated(Allen.genes)]
Allen.genes <- tolower(Allen.genes)
Allen.genes <- str_to_title(Allen.genes) # 3485 genes

length(intersect(genes.overlap, Allen.genes)) # 202/618
length(intersect(genes.bh.only, Allen.genes)) # 4/20
length(intersect(genes.jump.only, Allen.genes)) # 50/189


# Genes relative to olfactorybulb: BioGPS
Harmonizome_BioGPS <- fromJSON(paste(readLines("./data analysis/mouse olfactory bulb/validation sets/Harmonizome_BioGPS.json")))
Harmonizome_BioGPS <- Harmonizome_BioGPS[["associations"]]
BioGPS.genes <- NULL
for (i in 1:length(Harmonizome_BioGPS)){
  BioGPS.genes <- c(BioGPS.genes, Harmonizome_BioGPS[[i]][["gene"]][["symbol"]])
}
BioGPS.genes <- BioGPS.genes[!duplicated(BioGPS.genes)]
BioGPS.genes <- tolower(BioGPS.genes)
BioGPS.genes <- str_to_title(BioGPS.genes) # 2031 genes

length(intersect(genes.overlap, BioGPS.genes)) # 210/618
length(intersect(genes.bh.only, BioGPS.genes)) # 4/20
length(intersect(genes.jump.only, BioGPS.genes)) # 38/189

##--------------------------------------------------------------------
## GO enrichment analysis
##--------------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
# options(connectionObserver = NULL) # run if library(org.Mm.eg.db) failed

genes.all <- bitr(overlap, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
genes.jump <- bitr(genes_rep_jump, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.jump <- enrichGO(
  gene          = genes.jump$ENTREZID,
  universe      = genes.all$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "All" ,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE,
  pool = TRUE
)
sum(go.jump$p.adjust<0.01) # 846

genes.bh <- bitr(genes_rep_bh, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.bh <- enrichGO(
  gene          = genes.bh$ENTREZID,
  universe      = genes.all$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "All" ,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE,
  pool = TRUE
)
sum(go.bh$p.adjust<0.01) # 764

genes.maxp <- bitr(genes_rep_maxp, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.maxp <- enrichGO(
  gene          = genes.maxp$ENTREZID,
  universe      = genes.all$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "All" ,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE,
  pool = TRUE
)
sum(go.maxp$p.adjust<0.01) # 737

sig.go.jump <- filter(go.jump, p.adjust<.01)
sig.go.bh <- filter(go.bh, p.adjust<.01)
sig.go.maxp <- filter(go.maxp, p.adjust<.01)

library(VennDiagram)
venn.diagram(
  x = list(sig.go.jump@result$ID, sig.go.bh@result$ID, sig.go.maxp@result$ID),
  category.names = c("JUMP", "BH", "MaxP"),
  fill = c("#F4B183", "#56B4E9", "#8FBC8F"),
  filename = './data analysis/mouse olfactory bulb/venn_go_MOB.tiff',
  output = TRUE,
  margin = 0.02,
  cex = 1.5,
  cat.cex = 1.5,
  cat.dist = c(0.07, 0.05, 0.05),
)

sig.go.overlap <- intersect(sig.go.bh@result$ID, sig.go.jump@result$ID)
go.jump.only <- sig.go.jump[!sig.go.jump@result$ID%in%sig.go.overlap]

write.csv(go.jump,file = "./data analysis/mouse olfactory bulb/MOB_go_jump.csv", quote = FALSE)
write.csv(go.jump.only,file = "./data analysis/mouse olfactory bulb/MOB_go_jump_only.csv", quote = FALSE)

results <- go.jump@result
results <- results[sample(nrow(results)),]
results <- results[,c("ID","ONTOLOGY", "pvalue", "Count", "Description")]
results[,"Category"] = NA
results$Category[which(results$ONTOLOGY == "BP")] = 1
results$Category[which(results$ONTOLOGY == "CC")] = 2
results$Category[which(results$ONTOLOGY == "MF")] = 3
results <- results[order(results$Category),]

results[,"BP"] <- NA
results$BP[which(results$Category == 1)] = 1:sum(results$Category == 1)
results$BP[which(results$Category == 2)] = 1:sum(results$Category == 2)
results$BP[which(results$Category == 3)] = 1:sum(results$Category == 3)

don <- results %>% 
  
  # Compute group size
  group_by(ONTOLOGY) %>% 
  summarise(len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(len)-len) %>%
  dplyr::select(-len) %>%
  
  # Add this info to the initial dataset
  left_join(results, ., by=c("ONTOLOGY"="ONTOLOGY")) %>%
  
  # Add a cumulative position of group
  arrange(ONTOLOGY, BP) %>%
  mutate(BPcum=BP+tot)

axisdf = don %>% group_by(ONTOLOGY) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
annotated = c("G protein-coupled receptor signaling pathway", "trans-synaptic signaling",
              "neuron to neuron synapse", "regulation of nervous system development",
              "extracellular matrix", "dendritic tree", "neuronal cell body", 
              "regulation of neurotransmitter levels", "brain development")

ggplot(don, aes(x = BPcum, y=-log10(pvalue))) +
  geom_point(aes(color = as.factor(ONTOLOGY), size = Count), alpha=0.8) +
  scale_colour_manual(name="",  
                      values = c("BP"="Tan", "CC"="DarkSeaGreen", "MF"="#6496D2")) +
  scale_x_continuous(label = axisdf$ONTOLOGY, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +     # remove space between plot area and x axis
  geom_hline(yintercept = -log10(0.0025), color = '#545454', size= 1.2, linetype = "dashed") + 
  guides(color = "none") + 
  theme_bw() +
  theme( 
    legend.position = c(0.02,0.98),
    legend.justification = c(0.02,0.98),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, linetype = "solid")
  ) + 
  geom_text_repel(
    data = don[don$Description %in% annotated,],
    aes(label = Description),
    size = 4,
    segment.color = "black", show.legend = FALSE) # 5 * 7.5 inches/ 4.5*6.75
