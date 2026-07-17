# Install BiocManager if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_packages <- c(
  "GEOquery",
  "limma",
  "Biobase",
  "lumi",
  "illuminaHumanv3.db",
  "AnnotationDbi"
)
for(pkg in bioc_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
library(GEOquery)
library(limma)
library(Biobase)
BiocManager::install(
  "illuminaHumanv3.db",
  update = FALSE,
  ask = FALSE
)
library(illuminaHumanv3.db)
library(AnnotationDbi)
BiocManager::install(
  "lumi",
  update = FALSE,
  ask = FALSE
)
library(lumi)
BiocManager::version()
setwd("C:/Users/Mrinalini/Desktop/MultipleSclerosis")

# Create a new folder for GSE17048
dir.create(
  "GSE17048",
  showWarnings = FALSE
)

dir.create("raw_data", showWarnings = FALSE)

dir.create("normalized_data", showWarnings = FALSE)

dir.create("metadata", showWarnings = FALSE)

dir.create("QC_plots", showWarnings = FALSE)

dir.create("R_objects", showWarnings = FALSE)

dir.create("Results", showWarnings = FALSE)

dir.create(
  "Results/RRMS_vs_Control",
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  "Results/PPMS_vs_Control",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/SPMS_vs_Control",
  recursive = TRUE,
  showWarnings = FALSE
)

list.files()

series_matrix_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE17nnn/GSE17048/matrix/",
  "GSE17048_series_matrix.txt.gz"
)

download.file(
  series_matrix_url,
  "raw_data/GSE17048_series_matrix.txt.gz",
  mode = "wb",
  method = "libcurl"
)

file.exists(
  "raw_data/GSE17048_series_matrix.txt.gz"
)

dir.create(
  "Results/PPMS_vs_Control",
  recursive = TRUE,
  showWarnings = FALSE

  

  gse <- getGEO(
    filename = "raw_data/GSE17048_series_matrix.txt.gz",
    GSEMatrix = TRUE,
    AnnotGPL = FALSE,
    getGPL = FALSE
  )
  if (is.list(gse)) {
    gse <- gse[[1]]
  }
  metadata <- pData(gse)
  
  dim(metadata)
  
  head(metadata$title)
  
  colnames(metadata)


  write.csv(
    metadata,
    "metadata/GSE17048_original_metadata.csv",
    row.names = TRUE
  )
  
  saveRDS(
    metadata,
    "R_objects/GSE17048_original_metadata.rds"
  )


  expression_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    "GSE17nnn/GSE17048/suppl/",
    "GSE17048_non-normalized_data.txt.gz"
  )
  
  download.file(
    expression_url,
    "raw_data/GSE17048_non-normalized_data.txt.gz",
    mode = "wb",
    method = "libcurl"
  )

  file.exists(
    "raw_data/GSE17048_non-normalized_data.txt.gz"
  )
  

  raw_table <- read.delim(
    gzfile(
      "raw_data/GSE17048_non-normalized_data.txt.gz"
    ),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  dim(raw_table)
  
  head(raw_table[, 1:6])
  
  head(raw_table$ID_REF)
  
  probe_ids <- raw_table$ID_REF
  
  expression_raw <- as.matrix(
    raw_table[, -1]
  )
  
  storage.mode(expression_raw) <- "numeric"
  
  rownames(expression_raw) <- probe_ids
  
  dim(expression_raw)
  
  sum(is.na(expression_raw))
  
  sum(expression_raw <= 0)
  
  summary(
    as.vector(expression_raw)
  )
  
  
  write.csv(
    expression_raw,
    "raw_data/GSE17048_raw_expression_matrix.csv",
    row.names = TRUE
  )
  
  saveRDS(
    expression_raw,
    "R_objects/GSE17048_raw_expression_matrix.rds"
  )
  
  expression_log2 <- log2(
    expression_raw
  )
  summary(
    as.vector(expression_log2)
  )
  
  quantile(
    expression_log2,
    probs = c(
      0,
      0.01,
      0.25,
      0.50,
      0.75,
      0.99,
      1
    ),
    na.rm = TRUE
  )
  
  expression_norm <- normalizeBetweenArrays(
    expression_log2,
    method = "quantile"
  )
  
  dim(expression_norm)
  
  summary(
    as.vector(expression_norm)
  )
  write.csv(
    expression_norm,
    "normalized_data/GSE17048_quantile_normalized_expression.csv",
    row.names = TRUE
  )
  
  saveRDS(
    expression_norm,
    "R_objects/GSE17048_quantile_normalized_expression.rds"
  )  
  
  pdf(
    "QC_plots/GSE17048_before_normalization_boxplot.pdf",
    width = 14,
    height = 7
  )
  
  boxplot(
    expression_log2,
    outline = FALSE,
    xaxt = "n",
    main = "GSE17048 Before Quantile Normalization",
    ylab = "Log2 expression"
  )
  
  dev.off()
  
  
  
  pdf(
    "QC_plots/GSE17048_after_normalization_boxplot.pdf",
    width = 14,
    height = 7
  )
  
  boxplot(
    expression_norm,
    outline = FALSE,
    xaxt = "n",
    main = "GSE17048 After Quantile Normalization",
    ylab = "Normalized log2 expression"
  )
  
  dev.off()
  
  
  
  pdf(
    "QC_plots/GSE17048_before_normalization_density.pdf",
    width = 10,
    height = 7
  )
  
  plotDensities(
    expression_log2,
    legend = FALSE,
    main = "GSE17048 Before Quantile Normalization"
  )
  
  dev.off()
  pdf(
    "QC_plots/GSE17048_after_normalization_density.pdf",
    width = 10,
    height = 7
  )
  
  plotDensities(
    expression_norm,
    legend = FALSE,
    main = "GSE17048 After Quantile Normalization"
  )
  
  dev.off()
  
  while (!is.null(dev.list())) {
    dev.off()
  }  
  
  metadata_aligned <- metadata[
    match(
      colnames(expression_norm),
      metadata$description
    ),
    ,
    drop = FALSE
  ]
  
  
  all(
    colnames(expression_norm) ==
      metadata_aligned$description
  )
  
  metadata_aligned$group <- NA_character_
  
  metadata_aligned$group[
    grepl(
      "^healthy control",
      metadata_aligned$title,
      ignore.case = TRUE
    )
  ] <- "Control"
  
  metadata_aligned$group[
    grepl(
      "^RRMS",
      metadata_aligned$title,
      ignore.case = TRUE
    )
  ] <- "RRMS"
  
  metadata_aligned$group[
    grepl(
      "^SPMS",
      metadata_aligned$title,
      ignore.case = TRUE
    )
  ] <- "SPMS"
  
  metadata_aligned$group[
    grepl(
      "^PPMS",
      metadata_aligned$title,
      ignore.case = TRUE
    )
  ] <- "PPMS"
  
  table(
    metadata_aligned$group,
    useNA = "ifany"
  )
  
  write.csv(
    metadata_aligned,
    "metadata/GSE17048_aligned_metadata.csv",
    row.names = TRUE
  )
  
  saveRDS(
    metadata_aligned,
    "R_objects/GSE17048_aligned_metadata.rds"
  )
  
  pca <- prcomp(
    t(expression_norm),
    center = TRUE,
    scale. = FALSE
  )
  
  variance_explained <- round(
    100 * pca$sdev^2 /
      sum(pca$sdev^2),
    2
  )
  
  pca_df <- data.frame(
    Sample = colnames(expression_norm),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = metadata_aligned$group
  )
  
  plot(
    pca_df$PC1,
    pca_df$PC2,
    pch = 19,
    xlab = paste0(
      "PC1 (",
      variance_explained[1],
      "%)"
    ),
    ylab = paste0(
      "PC2 (",
      variance_explained[2],
      "%)"
    ),
    main = "PCA - GSE17048"
  )
  
  legend(
    "topright",
    legend = unique(pca_df$Group),
    pch = 19
  )
  
  write.csv(
    pca_df,
    "QC_plots/GSE17048_PCA_coordinates.csv",
    row.names = FALSE
  )
  
  pdf(
    "QC_plots/GSE17048_PCA_all_samples.pdf",
    width = 9,
    height = 7
  )
  
  plot(
    pca_df$PC1,
    pca_df$PC2,
    pch = 19,
    xlab = paste0(
      "PC1 (",
      variance_explained[1],
      "%)"
    ),
    ylab = paste0(
      "PC2 (",
      variance_explained[2],
      "%)"
    ),
    main = "PCA - GSE17048"
  )
  
  dev.off()
  
  
  probe_symbols <- mapIds(
    illuminaHumanv3.db,
    keys = rownames(expression_norm),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  probe_entrez <- mapIds(
    illuminaHumanv3.db,
    keys = rownames(expression_norm),
    column = "ENTREZID",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  probe_genenames <- mapIds(
    illuminaHumanv3.db,
    keys = rownames(expression_norm),
    column = "GENENAME",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  rrms_index <- metadata_aligned$group %in% c(
    "Control",
    "RRMS"
  )
  
  expr_rrms <- expression_norm[
    ,
    rrms_index,
    drop = FALSE
  ]
  
  meta_rrms <- metadata_aligned[
    rrms_index,
    ,
    drop = FALSE
  ]
  
  
  table(meta_rrms$group)
  
  all(
    colnames(expr_rrms) ==
      meta_rrms$description
  )
  
  group_rrms <- factor(
    meta_rrms$group,
    levels = c(
      "Control",
      "RRMS"
    )
  )
  
  design_rrms <- model.matrix(
    ~0 + group_rrms
  )
  
  colnames(design_rrms) <- levels(
    group_rrms
  )
  
  design_rrms
  
  
  fit_rrms <- lmFit(
    expr_rrms,
    design_rrms
  )
  
  contrast_rrms <- makeContrasts(
    RRMS_vs_Control = RRMS - Control,
    levels = design_rrms
  )
  
  fit_rrms2 <- contrasts.fit(
    fit_rrms,
    contrast_rrms
  )
  
  fit_rrms2 <- eBayes(
    fit_rrms2
  )
  
  deg_rrms <- topTable(
    fit_rrms2,
    coef = "RRMS_vs_Control",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  dim(deg_rrms)
  
  head(deg_rrms)
  
  deg_rrms$PROBEID <- rownames(
    deg_rrms
  )
  
  deg_rrms$SYMBOL <- probe_symbols[
    deg_rrms$PROBEID
  ]
  
  deg_rrms$GENENAME <- probe_genenames[
    deg_rrms$PROBEID
  ]
  
  deg_rrms$ENTREZID <- probe_entrez[
    deg_rrms$PROBEID
  ]
  
  
  deg_rrms <- deg_rrms[
    ,
    c(
      "PROBEID",
      "SYMBOL",
      "GENENAME",
      "ENTREZID",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    )
  ]
  
  write.csv(
    deg_rrms,
    "Results/RRMS_vs_Control/GSE17048_RRMS_vs_Control_all_results.csv",
    row.names = FALSE
  )
  
  saveRDS(
    deg_rrms,
    "R_objects/GSE17048_RRMS_vs_Control_all_results.rds"
  )
  
  rrms_fdr05 <- subset(
    deg_rrms,
    adj.P.Val < 0.05
  )
  
  nrow(rrms_fdr05)
  
  write.csv(
    rrms_fdr05,
    "Results/RRMS_vs_Control/GSE17048_RRMS_FDR05_DEGs.csv",
    row.names = FALSE
  )
  
  rrms_strict <- subset(
    deg_rrms,
    adj.P.Val < 0.05 &
      abs(logFC) > 0.5
  )
  
  nrow(rrms_strict)
  
  
  write.csv(
    rrms_strict,
    "Results/RRMS_vs_Control/GSE17048_RRMS_FDR05_logFC05_DEGs.csv",
    row.names = FALSE
  )
  
  rrms_exploratory <- subset(
    deg_rrms,
    P.Value < 0.05 &
      abs(logFC) > 0.5
  )
  
  nrow(rrms_exploratory)
  
  write.csv(
    rrms_exploratory,
    "Results/RRMS_vs_Control/GSE17048_RRMS_rawP05_logFC05_DEGs.csv",
    row.names = FALSE
  )
  
  rrms_paper_like <- subset(
    deg_rrms,
    P.Value < 0.05 &
      logFC > 0
  )
  
  nrow(rrms_paper_like)
  write.csv(
    rrms_paper_like,
    "Results/RRMS_vs_Control/GSE17048_RRMS_paper_like_upregulated.csv",
    row.names = FALSE
  )
  
  
  rrms_up <- subset(
    deg_rrms,
    P.Value < 0.05 &
      logFC > 0.5
  )
  
  rrms_down <- subset(
    deg_rrms,
    P.Value < 0.05 &
      logFC < -0.5
  )
  
  nrow(rrms_up)
  
  nrow(rrms_down)
  
  
  write.csv(
    rrms_up,
    "Results/RRMS_vs_Control/GSE17048_RRMS_upregulated_rawP05_logFC05.csv",
    row.names = FALSE
  )
  
  write.csv(
    rrms_down,
    "Results/RRMS_vs_Control/GSE17048_RRMS_downregulated_rawP05_logFC05.csv",
    row.names = FALSE
  )
  
  top50_rrms_up <- head(
    rrms_up[
      order(
        rrms_up$logFC,
        decreasing = TRUE
      ),
    ],
    50
  )
  
  top50_rrms_down <- head(
    rrms_down[
      order(
        rrms_down$logFC,
        decreasing = FALSE
      ),
    ],
    50
  )
  
  write.csv(
    top50_rrms_up,
    "Results/RRMS_vs_Control/GSE17048_RRMS_top50_upregulated.csv",
    row.names = FALSE
  )
  
  write.csv(
    top50_rrms_down,
    "Results/RRMS_vs_Control/GSE17048_RRMS_top50_downregulated.csv",
    row.names = FALSE
  )

  results_summary <- data.frame(
    Comparison = "RRMS vs Control",
    Total_probes = nrow(deg_rrms),
    FDR_05 = nrow(rrms_fdr05),
    FDR_05_logFC_05 = nrow(rrms_strict),
    RawP_05_logFC_05 = nrow(rrms_exploratory),
    Paper_like_upregulated = nrow(rrms_paper_like),
    Upregulated = nrow(rrms_up),
    Downregulated = nrow(rrms_down)
  )
  
  results_summary
  
  write.csv(
    results_summary,
    "Results/RRMS_vs_Control/GSE17048_RRMS_results_summary.csv",
    row.names = FALSE
  )
  rrms_unique <- rrms_paper_like[
    !duplicated(rrms_paper_like$SYMBOL),
  ]
rrms_unique
  
  


