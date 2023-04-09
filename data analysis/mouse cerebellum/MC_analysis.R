load("./data analysis/mouse cerebellum/MC.RData")

# Slide-seq
p1 <- sparkx1$res_mtest$combinedPval
names(p1) <- rownames(sparkx1$res_mtest)
# Slide-seqV2
p2 <- sparkx2$res_mtest$combinedPval
names(p2) <- rownames(sparkx2$res_mtest)
alpha = 0.05

overlap <- intersect(names(p1),names(p2))
pvals1 = p1[overlap]
pvals2 = p2[overlap]

# BH method
p1.bh <- p.adjust(pvals1, method = 'BH')
genes1.bh <- overlap[which(p1.bh<=alpha)]
p2.bh <- p.adjust(pvals2, method = 'BH')
genes2.bh <- overlap[which(p2.bh<=alpha)]
genes_rep_bh <- intersect(genes1.bh, genes2.bh)

# MaxP method
max_pvals <- apply(cbind(pvals1, pvals2), 1, max)
maxp.bh <- p.adjust(max_pvals, method = 'BH')
genes_rep_maxp <- overlap[which(maxp.bh<=alpha)]

# JUMP method
library(JUMPrcpp)
jump.obj <- JUMP(pvals1, pvals2, alpha)
jump.thr <- jump.obj$jump.thr
p.max <- jump.obj$p.max
genes_rep_jump <- overlap[which(p.max<=jump.thr)]

# Venn diagram
library(VennDiagram)
venn.diagram(
  x = list(genes_rep_jump, genes_rep_bh, genes_rep_maxp),
  category.names = c("JUMP", "BH", "MaxP"),
  fill = c("#F4B183", "#56B4E9", "#8FBC8F"),
  filename = './data analysis/mouse cerebellum/venn_MC.tiff',
  output = TRUE,
  margin = 0.02,
  scale = TRUE,
  cex = 1.5,
  cat.cex = 1.5,
  cat.dist = c(0.02, -0.01, -0.1),
)

save(overlap, pvals1, pvals2, genes_rep_bh, genes_rep_maxp, genes_rep_jump, 
     file = "./data analysis/mouse cerebellum/MC_results.RData")
