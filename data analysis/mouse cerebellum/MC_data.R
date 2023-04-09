#############################################################################
#          Mouse cerebellum data analysis (Slide-seq + Slide-seqV2)         #
#############################################################################
library(SPARK)
##---------------------------------------------------------------------------
## Slide-seq mouse cerebellum data
##---------------------------------------------------------------------------
counts1 <- read.csv("data analysis/mouse cerebellum/Slide-seq/Slide-seq_counts.csv", header = FALSE)
counts1 <- as.matrix(counts1)
genes1 <- as.vector(read.csv("data analysis/mouse cerebellum/Slide-seq/GeneNames.csv", header = FALSE))
location1 <- read.csv("data analysis/Mouse cerebellum/Slide-seq/Coordinates.csv", header = TRUE)
rownames(counts1) <- t(genes1)
colnames(location1) <- c("barcodes", "x", "y")
colnames(counts1) <- paste(location1$x, location1$y, sep = "x")
rownames(location1) <- colnames(counts1) # 18671 * 25551

mt_idx1  <- grep("mt-",rownames(counts1)) # 25 mitochondrial genes
counts1 <- counts1[-mt_idx1,] # 18646 * 25551

counts1 <- counts1[,which(colSums(counts1)>=50)]
counts1 <- counts1[which(rowSums(counts1)>0),] # 17550 * 16436
location1 <- location1[colnames(counts1),2:3]

# SPARK-X analysis
sparkx1 <- sparkx(counts1, location1, numCores=1, option="mixture") # 17481 * 14667

rm(mt_idx1, genes1)

##---------------------------------------------------------------------------
## Slide-seqV2 mouse cerebellum data
##---------------------------------------------------------------------------
load("./data analysis/mouse cerebellum/Slide-seqV2/SlideseqV2_ROI.rds")
# 20141 genes on 11626 beads
# remove mt genes
mt_idx2   <- grep("mt-",rownames(sp_count)) # 24 mitochondrial genes
counts2 <- as.matrix(sp_count[-mt_idx2,]) # 20117 genes on 11626 beads
location2 = as.data.frame(location)
rownames(location2) <- colnames(counts2)
rm(mt_idx2, sp_count, location)

##---------------------------------------------------------------------------
# SPARK-X analysis
sparkx2 <- sparkx(counts2, location2, numCores=1, option="mixture") # 20117 * 11626

##---------------------------------------------------------------------------
save(counts1, counts2, location1, location2, sparkx1, sparkx2, 
     file = "./data analysis/mouse cerebellum/MC.RData")

