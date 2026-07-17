if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")}

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)

setwd("C:/Users/Mrinalini/Desktop/MultipleSclerosis")
deg <- read.csv(
  "SigDegs_new.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)
colnames(deg)
View(deg)

deg_filtered <- subset(
  deg,
  P.Value < 0.05 &
    logFC > 0
)



genes <- unique(deg_filtered$Gene.Symbol)

genes <- genes[genes != ""]

genes <- na.omit(genes)

length(genes)

gene_ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


head(gene_ids)

dim(gene_ids)

go_bp <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

head(as.data.frame(go_bp))

write.csv(
  as.data.frame(go_bp),
  "GO_BP_results.csv",
  row.names = FALSE
)

dotplot(
  go_bp,
  showCategory = 20
)
kegg <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)


kegg <- setReadable(
  kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

write.csv(
  as.data.frame(kegg),
  "KEGG_results.csv",
  row.names = FALSE
)


dotplot(
  kegg,
  showCategory = 20
)

deg_down <- subset(
  deg,
  P.Value < 0.05 &
    logFC < 0
)

genes_down <- unique(deg_down$Gene.Symbol)

genes <- genes_down[genes_down != ""]

genes <- na.omit(genes)

length(genes)

gene_ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_ids)

dim(gene_ids)

go_bp <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)
head(as.data.frame(go_bp))

write.csv(
  as.data.frame(go_bp),
  "GO_Down_BP_results.csv",
  row.names = FALSE
)

dotplot(
  go_bp,
  showCategory = 20
)


kegg <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)


kegg <- setReadable(
  kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

write.csv(
  as.data.frame(kegg),
  "KEGG_down_results.csv",
  row.names = FALSE
)


dotplot(
  kegg,
  showCategory = 20
)



