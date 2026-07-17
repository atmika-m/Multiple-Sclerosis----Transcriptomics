# ==========================================
# GSE136411 VALIDATION: CIS vs CONTROL
# ==========================================

# Install packages only if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "GEOquery",
  "limma",
  "Biobase"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      update = FALSE,
      ask = FALSE
    )
  }
}

# Load libraries
library(GEOquery)
library(limma)
library(Biobase)


setwd(
  "C:/Users/Mrinalini/Desktop/MultipleSclerosis"
)

dir.create(
  "GSE136411",
  showWarnings = FALSE
)

gse <- getGEO(
  "GSE136411",
  GSEMatrix = TRUE,
  AnnotGPL = FALSE,
  getGPL = FALSE
)

length(gse)

annotation(gse[[1]])
annotation(gse[[2]])

getGEOSuppFiles(
  "GSE136411",
  baseDir = getwd(),
  makeDirectory = TRUE
)

list.files(
  "GSE136411",
  recursive = TRUE
)

valid_raw <- read.delim(
  "GSE136411/GSE136411_Matrix-merged-normalized-batch.corrected.txt.gz",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(valid_raw)

colnames(valid_raw)[1:10]

head(valid_raw$ID_REF)

rownames(valid_raw) <- valid_raw$ID_REF

valid_raw$ID_REF <- NULL

dim(valid_raw)

head(rownames(valid_raw))

head(colnames(valid_raw))

valid_raw[] <- lapply(
  valid_raw,
  function(x) {
    as.numeric(
      gsub(",", ".", x, fixed = TRUE)
    )
  }
)

str(valid_raw[, 1:5])

sum(is.na(valid_raw))

summary(
  as.vector(
    as.matrix(valid_raw)
  )
)

metadata1 <- pData(gse[[1]])
metadata2 <- pData(gse[[2]])

metadata <- rbind(
  metadata1,
  metadata2
)

nrow(metadata)

ncol(valid_raw)

gsm_number <- as.numeric(
  sub(
    "GSM",
    "",
    rownames(metadata)
  )
)

metadata <- metadata[
  order(gsm_number),
  ,
  drop = FALSE
]

head(
  data.frame(
    Matrix_column = colnames(valid_raw)[1:20],
    Metadata_title = metadata$title[1:20]
  ),
  20
)

metadata$group <- NA_character_

metadata$group[
  grepl(
    "^PBMC_HC",
    metadata$title,
    ignore.case = TRUE
  )
] <- "Control"

metadata$group[
  grepl(
    "^PBMC_CIS",
    metadata$title,
    ignore.case = TRUE
  )
] <- "CIS"

metadata$group[
  grepl(
    "^PBMC_RR",
    metadata$title,
    ignore.case = TRUE
  )
] <- "RRMS"

metadata$group[
  grepl(
    "^PBMC_PP",
    metadata$title,
    ignore.case = TRUE
  )
] <- "PPMS"

metadata$group[
  grepl(
    "^PBMC_SP",
    metadata$title,
    ignore.case = TRUE
  )
] <- "SPMS"

metadata$group[
  grepl(
    "^PBMC_OND",
    metadata$title,
    ignore.case = TRUE
  )
] <- "OND"

table(
  metadata$group,
  useNA = "ifany"
)

validation_expr <- as.matrix(valid_raw)

storage.mode(validation_expr) <- "numeric"

metadata$Matrix_column <- colnames(validation_expr)

metadata$Sample_ID <- sub(
  "[._].*$",
  "",
  metadata$Matrix_column
)

group_check <- tapply(
  metadata$group,
  metadata$Sample_ID,
  function(x) {
    length(
      unique(
        na.omit(x)
      )
    )
  }
)

table(group_check)

sample_ids <- unique(
  metadata$Sample_ID
)

validation_expr_collapsed <- sapply(
  sample_ids,
  function(id) {
    rowMeans(
      validation_expr[
        ,
        metadata$Sample_ID == id,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }
)

colnames(validation_expr_collapsed) <- sample_ids

metadata_collapsed <- metadata[
  !duplicated(metadata$Sample_ID),
  ,
  drop = FALSE
]

rownames(metadata_collapsed) <-
  metadata_collapsed$Sample_ID

metadata_collapsed <- metadata_collapsed[
  colnames(validation_expr_collapsed),
  ,
  drop = FALSE
]

all(
  colnames(validation_expr_collapsed) ==
    rownames(metadata_collapsed)
)

dim(validation_expr_collapsed)

table(metadata_collapsed$group)

metadata_cis <- metadata_collapsed[
  metadata_collapsed$group %in%
    c("Control", "CIS"),
  ,
  drop = FALSE
]

table(metadata_cis$group)

expr_cis <- validation_expr_collapsed[
  ,
  rownames(metadata_cis),
  drop = FALSE
]
all(
  colnames(expr_cis) ==
    rownames(metadata_cis)
)

dim(expr_cis)

group_cis <- factor(
  metadata_cis$group,
  levels = c(
    "Control",
    "CIS"
  )
)

design_cis <- model.matrix(
  ~0 + group_cis
)

colnames(design_cis) <- c(
  "Control",
  "CIS"
)

fit_cis <- lmFit(
  expr_cis,
  design_cis
)

contrast_cis <- makeContrasts(
  CIS_vs_Control = CIS - Control,
  levels = design_cis
)

fit_cis2 <- contrasts.fit(
  fit_cis,
  contrast_cis
)

fit_cis2 <- eBayes(
  fit_cis2
)


validation_deg_probe <- topTable(
  fit_cis2,
  coef = "CIS_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

validation_deg_probe$ID_REF <-
  rownames(validation_deg_probe)

validation_deg_probe <-
  validation_deg_probe[
    ,
    c(
      "ID_REF",
      setdiff(
        colnames(validation_deg_probe),
        "ID_REF"
      )
    )
  ]

dim(validation_deg_probe)

head(validation_deg_probe)

write.csv(
  validation_deg_probe,
  "GSE136411/Validation_CIS_vs_Control_Probe_Level.csv",
  row.names = FALSE
)

upregulated_cis <- validation_deg_probe[
  validation_deg_probe$adj.P.Val < 0.05 &
    validation_deg_probe$logFC > 0,
]
write.csv(
  upregulated_cis,
  "GSE136411/Validation_CIS_vs_Control_Upregulated_Probes.csv",
  row.names = FALSE
)

downregulated_cis <- validation_deg_probe[
  validation_deg_probe$adj.P.Val < 0.05 &
    validation_deg_probe$logFC < 0,
]
write.csv(
  downregulated_cis,
  "GSE136411/Validation_CIS_vs_Control_Downregulated_Probes.csv",
  row.names = FALSE
)

validation_summary <- data.frame(
  Comparison = "CIS vs Control",
  Total_probes = nrow(validation_deg_probe),
  Significant_FDR_05 = sum(validation_deg_probe$adj.P.Val < 0.05),
  Upregulated_FDR_05 = nrow(upregulated_cis),
  Downregulated_FDR_05 = nrow(downregulated_cis)
)

validation_summary

write.csv(
  validation_summary,
  "GSE136411/Validation_CIS_vs_Control_Summary.csv",
  row.names = FALSE
)

upregulated_rawp <- validation_deg_probe[
  validation_deg_probe$P.Value < 0.05 &
    validation_deg_probe$logFC > 0,
]

downregulated_rawp <- validation_deg_probe[
  validation_deg_probe$P.Value < 0.05 &
    validation_deg_probe$logFC < 0,
]

nrow(upregulated_rawp)
nrow(downregulated_rawp)



metadata_rrms <- metadata_collapsed[
  metadata_collapsed$group %in% c("Control", "RRMS"),
  ,
  drop = FALSE
]

table(metadata_rrms$group)

expr_rrms <- validation_expr_collapsed[
  ,
  rownames(metadata_rrms),
  drop = FALSE
]

all(colnames(expr_rrms) == rownames(metadata_rrms))

group_rrms <- factor(
  metadata_rrms$group,
  levels = c("Control", "RRMS")
)

design_rrms <- model.matrix(~0 + group_rrms)

colnames(design_rrms) <- c(
  "Control",
  "RRMS"
)
fit <- lmFit(
  expr_rrms,
  design_rrms
)

contrast <- makeContrasts(
  RRMS_vs_Control = RRMS - Control,
  levels = design_rrms
)

fit2 <- contrasts.fit(
  fit,
  contrast
)

fit2 <- eBayes(fit2)

fit2,


validation_rrms <- topTable(
  fit2,
  coef = "RRMS_vs_Control",
  number = Inf,
  adjust.method = "BH"
)

validation_rrms$ID_REF <- rownames(validation_rrms)

head(validation_rrms)

write.csv(
  validation_rrms,
  "Validation_RRMS_vs_Control.csv",
  row.names = FALSE
)

rrms_up <- subset(
  validation_rrms,
  adj.P.Val < 0.05 &
    logFC > 0
)
rrms_down <- subset(
  validation_rrms,
  adj.P.Val < 0.05 &
    logFC < 0
)

write.csv(
  rrms_up,
  "RRMS_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  rrms_down,
  "RRMS_Downregulated.csv",
  row.names = FALSE
)

gpl10558 <- getGEO("GPL10558")
gpl6104 <- getGEO("GPL6104")
colnames(gpl6104)

ann10558 <- Table(gpl10558)[,c("ID", "Symbol")]
ann6104 <- Table(gpl6104)[,c("ID", "Symbol")]

colnames(ann10558) <- c("ID_REF", "Gene.Symbol")
colnames(ann6104)  <- c("ID_REF", "Gene.Symbol")

annotation_all <- rbind(
  ann10558,
  ann6104
)

annotation_all <- annotation_all[
  !is.na(annotation_all$Gene.Symbol) &
    annotation_all$Gene.Symbol != "",
]

annotation_all <- unique(annotation_all)
getwd()
write.csv(
  annotation_all,
  "GSE136411_GPL_Annotation.csv",
  row.names = FALSE
)
