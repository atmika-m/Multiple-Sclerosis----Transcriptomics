############################################################
# GSE21942 RAW MICROARRAY PREPROCESSING
#
# Dataset:
# GSE21942
#
# Tissue:
# Peripheral blood mononuclear cells (PBMC)
#
# Platform:
# GPL570 - Affymetrix Human Genome U133 Plus 2.0 Array
#
# Output:
# 1. Raw AffyBatch object
# 2. RMA-normalized ExpressionSet
# 3. Expression matrix from exprs()
# 4. Sample metadata
############################################################


#############################
# 1. Install packages
#############################

# Install BiocManager if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
bioc_packages <- c(
  "GEOquery",
  "affy",
  "Biobase",
  "limma",
  "simpleaffy",
  "affyPLM",
  "hgu133plus2.db"
)

# Install missing Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# CRAN packages
cran_packages <- c(
  "R.utils"
)

# Install missing CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


#############################
# 2. Load libraries
#############################

library(GEOquery)
library(affy)
library(Biobase)
library(limma)
library(simpleaffy)
library(affyPLM)
library(hgu133plus2.db)
library(R.utils)


#############################
# 3. Create project folders
#############################

# Change this path to your preferred project location
project_dir <- "C:/Users/Mrinalini/Desktop/MultipleSclerosis"

raw_dir <- file.path(project_dir, "raw_data")
cel_dir <- file.path(project_dir, "CEL_files")
results_dir <- file.path(project_dir, "results")
plots_dir <- file.path(project_dir, "QC_plots")

dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

setwd(project_dir)

#############################
# 4. Download GEO metadata
#############################

gse_list <- getGEO(
  "GSE21942",
  GSEMatrix = TRUE,
  AnnotGPL = FALSE,
  getGPL = FALSE
)

# Check how many ExpressionSet objects were returned
length(gse_list)

# GSE21942 contains one platform, GPL570
gse <- gse_list[[1]]

# View basic dataset information
gse

# Extract sample metadata
metadata <- pData(gse)

# Inspect metadata
dim(metadata)
colnames(metadata)
head(metadata[, 1:min(10, ncol(metadata))])

# Save complete GEO metadata
write.csv(
  metadata,
  file = file.path(results_dir, "GSE21942_GEO_metadata.csv"),
  row.names = TRUE
)


#############################
# 5. Download raw CEL files
#############################

# This downloads the GEO supplementary files.
# For GSE21942, GEO provides GSE21942_RAW.tar.

getGEOSuppFiles(
  GEO = "GSE21942",
  makeDirectory = FALSE,
  baseDir = raw_dir
)

# List downloaded files
list.files(raw_dir, full.names = TRUE)


#############################
# 6. Extract the TAR archive
#############################

tar_file <- file.path(raw_dir, "GSE21942_RAW.tar")

if (!file.exists(tar_file)) {
  stop(
    "GSE21942_RAW.tar was not found. Check the contents of: ",
    raw_dir
  )
}

untar(
  tarfile = tar_file,
  exdir = cel_dir
)

# Check extracted files
list.files(cel_dir)[1:10]


#############################
# 7. Decompress CEL.gz files
#############################

# GEO Affymetrix CEL files are often stored as .CEL.gz files
gz_files <- list.files(
  cel_dir,
  pattern = "\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

length(gz_files)

if (length(gz_files) > 0) {
  for (gz_file in gz_files) {
    
    message("Decompressing: ", basename(gz_file))
    
    R.utils::gunzip(
      gz_file,
      remove = FALSE,
      overwrite = TRUE
    )
  }
}

# Identify decompressed CEL files
cel_files <- list.files(
  cel_dir,
  pattern = "\\.CEL$",
  full.names = TRUE,
  ignore.case = TRUE
)

length(cel_files)

# GSE21942 should contain 29 CEL files
if (length(cel_files) == 0) {
  stop("No CEL files were detected in: ", cel_dir)
}

print(basename(cel_files))


#############################
# 8. Read raw CEL files
#############################

raw_affy <- ReadAffy(
  filenames = cel_files
)

# Inspect raw AffyBatch object
raw_affy

# Number of probes and arrays
dim(exprs(raw_affy))

# Sample names
sampleNames(raw_affy)


#############################
# 9. Standardize sample names
#############################

# Extract GSM IDs from CEL filenames
gsm_ids <- sub(
  pattern = "^(GSM[0-9]+).*",
  replacement = "\\1",
  x = basename(cel_files)
)

# Confirm extraction
data.frame(
  original_file = basename(cel_files),
  extracted_GSM = gsm_ids
)

# Assign GSM IDs as sample names
sampleNames(raw_affy) <- gsm_ids

# Confirm sample names
sampleNames(raw_affy)


#############################
# 10. Match raw CEL files with GEO metadata
#############################

# Metadata row names should be GSM accessions
head(rownames(metadata))

# Reorder metadata so that it follows the CEL-file order
metadata_raw <- metadata[
  match(sampleNames(raw_affy), rownames(metadata)),
  ,
  drop = FALSE
]

# Check that every CEL file matched a metadata row
if (any(is.na(rownames(metadata_raw)))) {
  warning("Some CEL files did not match GEO metadata.")
}

# Stronger matching check
stopifnot(
  all(sampleNames(raw_affy) == rownames(metadata_raw))
)

# Add metadata to the AffyBatch object
pData(raw_affy) <- metadata_raw

# Save reordered metadata
write.csv(
  metadata_raw,
  file = file.path(results_dir, "GSE21942_metadata_CEL_order.csv"),
  row.names = TRUE
)


#############################
# 11. Inspect experimental groups
#############################

# GEO sample titles clearly identify controls and MS patients
table(metadata_raw$title)

# Create a broad Control/MS group from sample titles
group <- ifelse(
  grepl(
    pattern = "control",
    x = metadata_raw$title,
    ignore.case = TRUE
  ),
  "Control",
  ifelse(
    grepl(
      pattern = "MS patient",
      x = metadata_raw$title,
      ignore.case = TRUE
    ),
    "MS",
    NA
  )
)

group <- factor(
  group,
  levels = c("Control", "MS")
)

table(group, useNA = "ifany")

# Add group information to metadata
pData(raw_affy)$group <- group


#############################
# 12. Raw-data QC: boxplot
#############################

pdf(
  file.path(plots_dir, "01_raw_boxplot.pdf"),
  width = 12,
  height = 7
)

boxplot(
  raw_affy,
  main = "GSE21942: Raw CEL intensity distributions",
  las = 2,
  cex.axis = 0.6,
  outline = FALSE
)

dev.off()


#############################
# 13. Raw-data QC: density plot
#############################

pdf(
  file.path(plots_dir, "02_raw_density_plot.pdf"),
  width = 10,
  height = 7
)

hist(
  raw_affy,
  main = "GSE21942: Raw probe-intensity density"
)

dev.off()


#############################
# 14. RNA degradation assessment
#############################

rna_deg <- AffyRNAdeg(raw_affy)

pdf(
  file.path(plots_dir, "03_RNA_degradation_plot.pdf"),
  width = 10,
  height = 7
)

plotAffyRNAdeg(rna_deg)

dev.off()

# RNA degradation summary
rna_deg_summary <- summaryAffyRNAdeg(rna_deg)

write.csv(
  rna_deg_summary,
  file = file.path(results_dir, "GSE21942_RNA_degradation_summary.csv")
)


#############################
# 15. Optional QC: MA plots
#############################

pdf(
  file.path(plots_dir, "04_raw_MA_plots.pdf"),
  width = 10,
  height = 8
)

MAplot(
  raw_affy,
  pairs = TRUE,
  plot.method = "smoothScatter"
)

dev.off()


#############################
# 16. RMA preprocessing
#############################

# RMA performs:
# 1. Background correction
# 2. Quantile normalization
# 3. Probe-set summarization
# 4. Log2 transformation

eset_rma <- rma(raw_affy)

# Inspect normalized ExpressionSet
eset_rma


#############################
# 17. Extract expression matrix using exprs()
#############################

expression_matrix <- exprs(eset_rma)

# Dimensions:
# rows    = Affymetrix probe-set IDs
# columns = samples
dim(expression_matrix)

# View first few values
expression_matrix[1:6, 1:6]

# Confirm range of normalized values
summary(as.vector(expression_matrix))


#############################
# 18. Add metadata to normalized ExpressionSet
#############################

# Make sure metadata and normalized matrix have the same order
stopifnot(
  all(colnames(expression_matrix) == rownames(metadata_raw))
)

pData(eset_rma) <- metadata_raw
pData(eset_rma)$group <- group

# Verify
head(pData(eset_rma)[, c("title", "group")])


#############################
# 19. Normalized-data QC: boxplot
#############################

pdf(
  file.path(plots_dir, "05_RMA_normalized_boxplot.pdf"),
  width = 12,
  height = 7
)

boxplot(
  expression_matrix,
  main = "GSE21942: RMA-normalized expression",
  las = 2,
  cex.axis = 0.6,
  outline = FALSE
)

dev.off()


#############################
# 20. Normalized-data QC: density plot
#############################

pdf(
  file.path(plots_dir, "06_RMA_normalized_density.pdf"),
  width = 10,
  height = 7
)

plotDensities(
  expression_matrix,
  legend = FALSE,
  main = "GSE21942: RMA-normalized expression density"
)

dev.off()


#############################
# 21. PCA after normalization
#############################

pca <- prcomp(
  t(expression_matrix),
  center = TRUE,
  scale. = FALSE
)

pca_variance <- round(
  100 * (pca$sdev^2 / sum(pca$sdev^2)),
  2
)

pca_data <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = group
)

pdf(
  file.path(plots_dir, "07_RMA_PCA_plot.pdf"),
  width = 8,
  height = 7
)

plot(
  pca_data$PC1,
  pca_data$PC2,
  pch = 19,
  xlab = paste0("PC1: ", pca_variance[1], "% variance"),
  ylab = paste0("PC2: ", pca_variance[2], "% variance"),
  main = "GSE21942 PCA after RMA normalization"
)

text(
  pca_data$PC1,
  pca_data$PC2,
  labels = pca_data$Sample,
  pos = 3,
  cex = 0.55
)

legend(
  "topright",
  legend = levels(group),
  pch = 19
)

dev.off()


#############################
# 22. Sample correlation heatmap
#############################

sample_correlation <- cor(
  expression_matrix,
  method = "pearson"
)

pdf(
  file.path(plots_dir, "08_sample_correlation_heatmap.pdf"),
  width = 11,
  height = 10
)

heatmap(
  sample_correlation,
  symm = TRUE,
  margins = c(8, 8),
  main = "GSE21942 sample correlation"
)

dev.off()


#############################
# 23. Save normalized outputs
#############################

# Save complete R objects
saveRDS(
  raw_affy,
  file = file.path(results_dir, "GSE21942_raw_AffyBatch.rds")
)

saveRDS(
  eset_rma,
  file = file.path(results_dir, "GSE21942_RMA_ExpressionSet.rds")
)

saveRDS(
  expression_matrix,
  file = file.path(results_dir, "GSE21942_RMA_expression_matrix.rds")
)

# Save expression matrix as CSV
write.csv(
  expression_matrix,
  file = file.path(results_dir, "GSE21942_RMA_expression_matrix.csv"),
  row.names = TRUE
)

# Save PCA data
write.csv(
  pca_data,
  file = file.path(results_dir, "GSE21942_PCA_coordinates.csv"),
  row.names = FALSE
)


#############################
# 24. Final checks
#############################

cat("\nPreprocessing completed successfully.\n")

cat(
  "\nExpression matrix dimensions:",
  nrow(expression_matrix),
  "probe sets x",
  ncol(expression_matrix),
  "arrays\n"
)

cat("\nGroup counts:\n")
print(table(group))

cat("\nExpression-value summary:\n")
print(summary(as.vector(expression_matrix)))

cat(
  "\nNormalized expression matrix saved at:\n",
  file.path(results_dir, "GSE21942_RMA_expression_matrix.csv"),
  "\n"
)

dim(expression_matrix)
dim(metadata_raw)

colnames(metadata_raw)

length(metadata_raw$group)
length(group)

metadata_raw$title
group <- ifelse(
  grepl("control", metadata_raw$title, ignore.case = TRUE),
  "Control",
  ifelse(
    grepl("MS patient", metadata_raw$title, ignore.case = TRUE),
    "MS",
    NA
  )
)

metadata_raw$group <- group
table(metadata_raw$group, useNA = "ifany")
length(metadata_raw$group)
colnames(expression_matrix)
rownames(metadata_raw)
all(colnames(expression_matrix) == rownames(metadata_raw))

library(ggplot2)

# Confirm dimensions before PCA
stopifnot(ncol(expression_matrix) == nrow(metadata_raw))

# PCA requires samples as rows and genes/probes as columns
pca <- prcomp(
  t(expression_matrix),
  center = TRUE,
  scale. = FALSE
)

# Percentage variance explained
variance_explained <- round(
  100 * pca$sdev^2 / sum(pca$sdev^2),
  2
)

# Create plotting dataframe
plot_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = metadata_raw$group,
  stringsAsFactors = FALSE
)

# Check plot dataframe
dim(plot_df)
head(plot_df)

# PCA plot
p <- ggplot(
  plot_df,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    label = Sample
  )
) +
  geom_point(size = 3) +
  labs(
    title = "PCA of GSE21942 after RMA normalization",
    x = paste0("PC1: ", variance_explained[1], "% variance"),
    y = paste0("PC2: ", variance_explained[2], "% variance")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p)

library(ggrepel)

p_labelled <- ggplot(
  plot_df,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    label = Sample
  )
) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "PCA of GSE21942 after RMA normalization",
    x = paste0("PC1: ", variance_explained[1], "% variance"),
    y = paste0("PC2: ", variance_explained[2], "% variance")
  ) +
  theme_bw()

print(p_labelled)
ggsave(
  filename = file.path(plots_dir, "GSE21942_PCA.png"),
  plot = p_labelled,
  width = 9,
  height = 7,
  dpi = 300
)
library(ggplot2)

boxplot_df <- data.frame(
  Expression = as.vector(expression_matrix),
  Sample = rep(
    colnames(expression_matrix),
    each = nrow(expression_matrix)
  )
)

p_box <- ggplot(
  boxplot_df,
  aes(
    x = Sample,
    y = Expression
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.3
  ) +
  labs(
    title = "GSE21942 RMA-normalized expression distributions",
    x = "Sample",
    y = "Log2 expression"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    )
  )

print(p_box)
ggsave(
  filename = file.path(plots_dir, "GSE21942_normalized_boxplot.png"),
  plot = p_box,
  width = 12,
  height = 7,
  dpi = 300
)
cat("Expression matrix dimensions:\n")
print(dim(expression_matrix))

cat("\nMetadata dimensions:\n")
print(dim(metadata_raw))

cat("\nGroup column exists:\n")
print("group" %in% colnames(metadata_raw))

cat("\nNumber of group labels:\n")
print(length(metadata_raw$group))

cat("\nGroup counts:\n")
print(table(metadata_raw$group, useNA = "ifany"))

cat("\nSample order identical:\n")
print(all(colnames(expression_matrix) == rownames(metadata_raw)))

cat("\nMissing group labels:\n")
print(sum(is.na(metadata_raw$group)))




#Convert Probe IDs into Gene Symbols

library(hgu133plus2.db)
library(AnnotationDbi)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = rownames(expression_matrix),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

expr <- data.frame(
  ProbeID = rownames(expression_matrix),
  Gene = gene_symbols,
  expression_matrix,
  check.names = FALSE
)
expr <- expr[!is.na(expr$Gene), ]
expr <- expr[expr$Gene != "", ]
dim(expr)
library(limma)

expr_gene <- avereps(
  as.matrix(expr[, -(1:2)]),
  ID = expr$Gene
)
keep <- rowMeans(expr_gene) > 5

expr_gene <- expr_gene[keep, ]
dim(expr_gene)
write.csv(
  expr_gene,
  "GSE21942_GeneExpression.csv"
)
colnames(metadata_raw)
group <- factor(metadata_raw$group)

table(group)
design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

design
library(limma)

fit <- lmFit(expr_gene, design)
contrast.matrix <- makeContrasts(
  
  MSvsControl = MS - Control,
  
  levels = design
  
)
fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)
deg <- topTable(
  fit2,
  coef = "MSvsControl",
  number = Inf,
  adjust.method = "BH"
)
head(deg)
sig_deg <- subset(
  deg,
  adj.P.Val < 0.05 &
    abs(logFC) >0.5
)
upregulated <- sig_deg[sig_deg$logFC > 0.5, ]

downregulated <- sig_deg[sig_deg$logFC < -0.5, ]
nrow(upregulated)

nrow(downregulated)
write.csv(sig_deg, "Significant_DEGs.csv")
nrow(sig_deg)
write.csv(upregulated, "Upregulated.csv")

write.csv(downregulated, "Downregulated.csv")
