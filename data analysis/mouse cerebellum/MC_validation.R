rm(list = ls())
load("./data analysis/mouse cerebellum/MC.RData")
load("./data analysis/mouse cerebellum/MC_Results.RData")

##----------------------------------------------------------
## Moran's I
##----------------------------------------------------------
source("./funcs/ValidationFuncs.R")
library(spdep)

genes.overlap <- intersect(intersect(genes_rep_bh, genes_rep_maxp), genes_rep_jump)
genes.bh.only  <- genes_rep_bh[!genes_rep_bh%in%genes.overlap]
genes.jump.only <- genes_rep_jump[!genes_rep_jump%in%genes.overlap]

## MI based on Slide-seq data
stat.res1 <- corrValue(counts1, location1)
MI1 <- stat.res1
mi_all1 <- MI1[overlap]
mi_jump_only1 <- MI1[genes.jump.only]

# Moran's I for SVGs additionally identified by JUMP
mean(mi_all1); mean(mi_jump_only1); median(mi_all1); median(mi_jump_only1)
# 0.0011; 0.0079; 1.77e-4; 0.0061
factor1 <- factor(rep(c("JUMP only", "All (Slide-seq)"), 
                      times = c(length(genes.jump.only), length(overlap))))
dataset1 <- data.frame(value = c(mi_jump_only1, mi_all1), group = factor1)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset1, col = c("antiquewhite1", "lightSalmon"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6) # 5 * 7.5 inches

## MI based on Slide-seqV2 data
stat.res2 <- corrValue(counts2, location2)
MI2 <- stat.res2
mi_all2 <- MI2[overlap]
mi_jump_only2 <- MI2[genes.jump.only]

# Moran's I for SVGs additionally identified by JUMP
mean(mi_all2); mean(mi_jump_only2); median(mi_all2); median(mi_jump_only2)
# 0.0031; 0.0176; 0.0011; 0.0143
factor2 <- factor(rep(c("JUMP only", "All (Slide-seqV2)"), 
                      times = c(length(genes.jump.only), length(overlap))))
dataset2 <- data.frame(value = c(mi_jump_only2, mi_all2), group = factor2)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset2, col = c("antiquewhite1", "lightSalmon"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.7) # 5 * 7.5

save(stat.res1, stat.res2, file = "./data analysis/mouse cerebellum/MC_MI.RData")

##--------------------------------------------------------------------
## Summarized spatial expression patterns
##--------------------------------------------------------------------
library(amap)
source("R/PlotFuncs.R")

## Summarized patterns in Slide-seq data
vst_count1 <- var_stabilize(counts1) # R function in funcs.R
sig_vst_count1 <- vst_count1[genes_rep_jump, ]
sig_vst_res1 <- t(apply(sig_vst_count1, 1, LMReg, T = log(apply(counts1, 2, sum))))

hc1 <- hcluster(sig_vst_res1, method = "euc", link = "ward", nbproc = 1, 
                doubleprecision = TRUE)
numC <- 3
memb1 <- cutree(hc1, k = numC)
cent1 <- NULL
for (k in 1:numC) {
  cent1 <- cbind(cent1, colMeans(sig_vst_res1[memb1 == k, , drop = FALSE]))
}
position_cord1 <- location1
rownames(position_cord1) <- rownames(cent1)

rel_cent1 <- t(apply(cent1, 1, relative_func))
# rel_cent1 <- rel_cent1[,c(2,1,3)]
pd1 <- setNames(cbind.data.frame(position_cord1, rel_cent1), 
                c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP1 <- lapply(1:numC, function(x) {
  pattern_plot3(pd1, x, xy = T, main = T, titlesize = 9.5)
})
grid.arrange(grobs = MBP1[numC:1], nrow = 3) # 12.5 * 37.5

## Summarized patterns in Slide-seqV2 data
vst_count2 <- var_stabilize(counts2) 
sig_vst_count2 <- vst_count2[genes_rep_jump, ]
sig_vst_res2 <- t(apply(sig_vst_count2, 1, LMReg, T = log(apply(counts2, 2, sum))))

hc2 <- hcluster(sig_vst_res2, method = "euc", link = "ward", nbproc = 1, 
                doubleprecision = TRUE)
numC <- 3
memb2 <- cutree(hc2, k = numC)
cent2 <- NULL
for (k in 1:numC) {
  cent2 <- cbind(cent2, colMeans(sig_vst_res2[memb2 == k, , drop = FALSE]))
}
position_cord2 <- location2
rownames(position_cord2) <- rownames(cent2)

rel_cent2 <- t(apply(cent2, 1, relative_func))
rel_cent2 <- rel_cent2[,c(3,2,1)]
pd2 <- setNames(cbind.data.frame(position_cord2, rel_cent2), 
                c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP2 <- lapply(1:numC, function(x) {
  pattern_plot3(pd2, x, xy = T, main = T, titlesize = 9.5)
})
grid.arrange(grobs = MBP2[numC:1], nrow = 3) # 12.5 * 37.5

##--------------------------------------------------------------------
## Spatial expression patterns of randomly selected SVGs
##--------------------------------------------------------------------
source("R/PlotFuncs.R")

# three representative genes corresponding to the main patterns
gene_plot <- c("Pcp2", "Mbp", "Snap25")

# 24 genes randomly selected from 169 SVGs additionally identified by JUMP compared to the overlaps
randi <- sample(length(genes.jump.only), 24, replace = FALSE)
gene_plot <- genes.jump.only[randi]

# Spatial expression patterns based on Slide-seq data
vst_count1 <- var_stabilize(counts1) # R function in funcs.R
sig_vst_ct1 <- vst_count1[gene_plot, ]
rel_vst_ct1 <- apply(sig_vst_ct1, 1, relative_func)
pltdat1 <- cbind.data.frame(location1,rel_vst_ct1)
genetitle <- gene_plot
pp1 <- lapply(1:(ncol(pltdat1)-2),
              function(x){pattern_plot3(pltdat1,x,main=T,titlesize=9.5,title=genetitle[x])})
grid.arrange(grobs=pp1, nrow=3) # 37.5 * 100 inches

# Spatial expression patterns based on Slide-seqV2 data
vst_count2 <- var_stabilize(counts2) # R function in funcs.R
sig_vst_ct2 <- vst_count2[gene_plot, ]
rel_vst_ct2 <- apply(sig_vst_ct2, 1, relative_func)
pltdat2 <- cbind.data.frame(location2,rel_vst_ct2)
genetitle <- gene_plot
pp2 <- lapply(1:(ncol(pltdat2)-2),
              function(x){pattern_plot3(pltdat2,x,main=T,titlesize=9.5,title=genetitle[x])})
grid.arrange(grobs=pp2, nrow=3) # 37.5 * 100 inches

##--------------------------------------------------------------------
## Mouse cerebellum related studies validation
##--------------------------------------------------------------------
# Genes related to the cerebellum in Harmonizome database (3000 genes)
library(rjson)
library(stringr)
Harmonizome_Allen <- fromJSON(paste(readLines("./data analysis/mouse cerebellum/validation sets/Harmonizome_Allen Atlas.json")))
Harmonizome_Allen_cortex <- fromJSON(paste(readLines("./data analysis/mouse cerebellum/validation sets/Harmonizome_Allen_cortex.json")))
Harmonizome_Allen_hemisphere <- fromJSON(paste(readLines("./data analysis/mouse cerebellum/validation sets/Harmonizome_Allen_hemisphere.json")))

Harmonizome_Allen <- Harmonizome_Allen[["associations"]]
Harmonizome_Allen_cortex <- Harmonizome_Allen_cortex[["associations"]]
Harmonizome_Allen_hemisphere <- Harmonizome_Allen_hemisphere[["associations"]]

Harmonizome.genes <- NULL
for (i in 1:length(Harmonizome_Allen)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_Allen_cortex)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen_cortex[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_Allen_hemisphere)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen_hemisphere[[i]][["gene"]][["symbol"]])
}

Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)]
Harmonizome.genes <- tolower(Harmonizome.genes)
Harmonizome.genes <- str_to_title(Harmonizome.genes) # 3000 genes

length(intersect(genes_rep_maxp, Harmonizome.genes)) # 84/279
length(intersect(genes.bh.only, Harmonizome.genes)) # 28/115
length(intersect(genes.jump.only, Harmonizome.genes)) # 41/169

# Kozareva et al. (3976 genes)
library(readxl)
cluster_genes <- read_xlsx("./data analysis/mouse cerebellum/validation sets/Kozareva et al.xlsx",
                              col_names = TRUE)
cluster_genes = cluster_genes[which(abs(cluster_genes$logFC)>=0.5),] 
cluster_genes = as.vector(unique(cluster_genes$gene)) # 3976

length(intersect(genes_rep_maxp, cluster_genes)) # 214/279
length(intersect(genes.bh.only, cluster_genes)) # 73/115
length(intersect(genes.jump.only, cluster_genes)) # 107/169

##--------------------------------------------------------------------
## GO enrichment
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
sum(go.jump$p.adjust<0.01) # 452

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
sum(go.bh$p.adjust<0.01) # 418

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
sum(go.maxp$p.adjust<0.01) # 282

sig.go.jump <- filter(go.jump, p.adjust<.01)
sig.go.bh <- filter(go.bh, p.adjust<.01)
sig.go.maxp <- filter(go.maxp, p.adjust<.01)

library(VennDiagram)
venn.diagram(
  x = list(sig.go.jump@result$ID, sig.go.bh@result$ID, sig.go.maxp@result$ID),
  category.names = c("JUMP", "BH", "MaxP"),
  fill = c("#F4B183", "#56B4E9", "#8FBC8F"),
  filename = './data analysis/mouse cerebellum/validation sets/venn_go_MC.tiff',
  output = TRUE,
  margin = 0.02,
  cex = 1.5,
  cat.cex = 1.5,
  cat.dist = c(0.07, 0.05, 0.05),
)

sig.go.overlap <- intersect(sig.go.bh@result$ID, sig.go.jump@result$ID) # 394
go.jump.only <- sig.go.jump[!sig.go.jump@result$ID%in%sig.go.overlap] # 58

write.csv(go.jump,file = "./data analysis/mouse cerebellum/validation sets/MC_go_jump.csv", quote = FALSE)
write.csv(go.jump.only,file = "./data analysis/mouse cerebellum/validation sets/MC_go_jump_only.csv", quote = FALSE)

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
# MaxP_S
annotated = c("glutamatergic synapse", "neurotransmitter transport", "synapse organization",
              "exocytic process", "parallel fiber to Purkinje cell synapse",
              "GABA-ergic synapse", "regulation of neuron projection development")

ggplot(don, aes(x = BPcum, y=-log10(pvalue))) +
  geom_point(aes(color = as.factor(ONTOLOGY), size = Count), alpha=0.8) +
  scale_colour_manual(name="",  
                      values = c("BP"="Tan", "CC"="DarkSeaGreen", "MF"="#6496D2")) +
  scale_x_continuous(label = axisdf$ONTOLOGY, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +     # remove space between plot area and x axis
  geom_hline(yintercept = -log10(0.0012), color = '#545454', size= 1.2, linetype = "dashed") + 
  guides(color = "none") + 
  theme_bw() +
  theme( 
    legend.position = c(0.02,0.98),
    legend.justification = c(0.02,0.98),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12.5),
    axis.text.y = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, linetype = "solid")
  ) + 
  geom_text_repel(
    data = don[don$Description %in% annotated,],
    aes(label = Description),
    size = 4,
    segment.color = "black", show.legend = FALSE) # 5 * 7.5 inches
