# BLOMOPARKED - (BLOod Multi-Omics for PARKinson's disease Early Detection)

# Project summary: 
# Multi-Omics data integration of RNA-Seq and smallRNA-Seq data
# from blood samples of Parkinson's disease patients at different stages of the disease,
# using Differential Expression Analysis, correlation analysis, and Machine Learning.
# with the aim of identifying early stage biomarkers or key molecular players in the disease progression.

# Dataset: PPMI (Parkinson's Progression Markers Initiative)

# Table of Contents:
# 00. Load necessary libraries
# 01. Genetic Data Preparation (RNA-seq & miRNA raw counts)
# 02. Metadata Preparation (clinical data, blood chemistry, phenotypes)
# 03. Data Cleaning & Harmonization (sample alignment, filtering)
# 04. Genetic Data Preprocessing (annotation, normalization)
# 05. Differential Gene Expression Analysis (DGE: DEGs & DEmiRNAs, Volcano plots, UpSet plot)
# 06. miRNA-mRNA Integration (batch correction, correlation analysis, target prediction)
# 07. Machine Learning (8 classifiers: RF, SVM, GLMNet, GBM, kNN, NB, LDA, NNet)

## 00. Load necessary libraries ----
library(readr)
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(R.utils)
library(ggplot2)
library(AnnotationHub)
library(ensembldb)
library(AnnotationDbi)
library(HGNChelper)
library(MIRit)
library(UpSetR)
library(doParallel)
library(BiocParallel)
library(gwasrapidd)
library(igraph)
library(patchwork)
library(caret)
library(randomForest)
library(e1071)
library(glmnet)
library(gbm)
library(class)
library(klaR)
library(MASS)
library(kernlab)
library(nnet)
library(pROC)
library(gridExtra)
## 01. Genetic Data Preparation ----

counts_dir_RNAseq <- "PPMI_RNAseq_IR3_Analysis/counts"
# Get the list of featureCounts files
fc_files <- list.files(
  counts_dir_RNAseq,
  pattern = "featureCounts",
  full.names = TRUE
)
# Extract Unique_ID from filename (PATNO + EVENT_ID separated by "_")
get_unique_id <- function(path) {
  parts <- strsplit(basename(path), "\\.")[[1]]
  patno <- parts[2]
  event <- parts[3]
  paste0(patno, "_", event)
}

unique_ids_rna <- vapply(fc_files, get_unique_id, character(1))
names(fc_files) <- unique_ids_rna
# Read the gene list from a file
genes <- read.table(
  fc_files[1],
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  quote = ""
)[[1]]
# Read the counts column
read_last <- function(f) {
  fread(
    f,
    select = 7,
    skip = "#",
    colClasses = "integer"
  )[[1]]
}
# Combining the counts columns from all files
counts_matrix_RNAseq <- do.call(cbind, lapply(fc_files, read_last))
# Prepend the Geneid column
counts_matrix_RNAseq <- cbind(Geneid = genes, counts_matrix_RNAseq)

# Optionally, write to CSV
# write.csv(counts_matrix_RNAseq, file = "RNAseq_counts.csv", row.names = FALSE)


# Edit the single sncRNA-seq countmatrix file keeping only the Unique_ID pattern

counts_dir_sncRNA <- "PPMI_sncRNAcounts/counts"

counts_matrix_miRNA <- fread(file.path(counts_dir_sncRNA, "mirna_quantification_matrix_raw.csv"))

# Extract Unique_IDs
colnames_miRNA <- colnames(counts_matrix_miRNA)[-1]
get_unique_id_miRNA <- function(name) {
  parts <- strsplit(name, "\\.")[[1]]
  patno <- parts[2]
  event <- parts[3]
  paste0(patno, "_", event)
}

unique_ids_miRNA <- vapply(colnames_miRNA, get_unique_id_miRNA, character(1))
# Update column names
colnames(counts_matrix_miRNA)[-1] <- unique_ids_miRNA

counts_matrix_miRNA <- as.data.frame(counts_matrix_miRNA)

# Optionally, write to CSV
# fwrite(counts_matrix_miRNA, file = "miRNA_counts.csv")

## 02. Metadata preparation ----

metadata <- readxl::read_excel("PPMI_Curated_Data_Cut_Public_20250714.xlsx")
# Add 'Unique_ID' column 
metadata <- metadata %>%
  mutate(Unique_ID = paste0(PATNO, "_", EVENT_ID))
# Set 'Unique_ID' as row names
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Unique_ID
# Optionally, write to CSV
# write.csv(metadata, file = "metadata.csv", row.names = FALSE)

# Filter out observations with "Not NSD" in NSD_STAGE, and keep NAs only for COHORT 2 or 4 (Healthy Control or Prodromal)
metadata_nsd <- dplyr::filter(metadata, NSD_STAGE != "Not NSD" | (is.na(NSD_STAGE) & COHORT %in% c(2, 4)))

# Optionslly, save the filtered metadata
# write.csv(as.data.frame(nsd_stage_counts), file = "PPMI_Metadata_NSD_Stage_Counts.csv", row.names = TRUE)

# Collapse NSD_STAGE into broader categories
metadata_nsd$New_Stage <- NA
metadata_nsd$New_Stage[metadata_nsd$COHORT == 2] <- "Healthy"
metadata_nsd$New_Stage[metadata_nsd$COHORT == 4] <- "Pre-clinical"
metadata_nsd$New_Stage[metadata_nsd$NSD_STAGE %in% c("1a", "1b", "2a")] <- "Pre-clinical"
metadata_nsd$New_Stage[metadata_nsd$NSD_STAGE == "2b"] <- "Mild"
metadata_nsd$New_Stage[metadata_nsd$NSD_STAGE == "3"] <- "Moderate"
metadata_nsd$New_Stage[metadata_nsd$NSD_STAGE %in% c("4", "5", "6")] <- "Severe"

# Optionally, save the new metadata
# write.csv(metadata_nsd, file = "metadata_newstage.csv", row.names = FALSE)

# Blood chemistry data preparation
blood_chem <- read.delim('Blood_Chemistry___Hematology-Archived_05Sep2025.csv', sep = ",", header = TRUE)


# Create Unique_ID and filter for GI/L units, keeping relevant columns
blood_chem$Unique_ID <- paste0(blood_chem$PATNO, "_", blood_chem$EVENT_ID)
blood_chem_sub <- blood_chem[blood_chem$LSIUNIT == "GI/L", c("Unique_ID", "PATNO", "EVENT_ID", "LTSTNAME", "LSIRES")]

blood_chem_wide <- pivot_wider(blood_chem_sub, names_from = LTSTNAME, values_from = LSIRES)
# Add the new columns (WBC, Platelets, Neutrophils, Lymphocytes) to the metadata table
meta_expanded <- left_join(metadata_nsd, blood_chem_wide[, c("Unique_ID", "WBC", "Platelets", "Neutrophils", "Lymphocytes")], by = c("Unique_ID"))

# Convert the new columns to numeric, handling lists and NULLs
meta_expanded <- meta_expanded %>%
  mutate(
    WBC = as.numeric(sapply(WBC, function(x) if (is.null(x) || length(x) == 0) NA else x[[1]])),
    Platelets = as.numeric(sapply(Platelets, function(x) if (is.null(x) || length(x) == 0) NA else x[[1]])),
    Neutrophils = as.numeric(sapply(Neutrophils, function(x) if (is.null(x) || length(x) == 0) NA else x[[1]])),
    Lymphocytes = as.numeric(sapply(Lymphocytes, function(x) if (is.null(x) || length(x) == 0) NA else x[[1]]))
  )

# Impute missing values hierarchically: by PATNO, then by New_Stage, then globally
impute_column <- function(x, group = NULL) {
  idx_na <- is.na(x)
  if (!any(idx_na)) return(x)
  x[idx_na] <- mean(x, na.rm = TRUE)
  return(x)
}

meta_imputed <- meta_expanded %>%
  group_by(PATNO) %>%
  mutate(across(c(WBC, Platelets, Neutrophils, Lymphocytes, BMI), impute_column)) %>%
  ungroup() %>%
  group_by(New_Stage) %>%
  mutate(across(c(WBC, Platelets, Neutrophils, Lymphocytes, BMI), impute_column)) %>%
  ungroup() %>%
  mutate(across(c(WBC, Platelets, Neutrophils, Lymphocytes, BMI), impute_column))

# Add a column with Neutrophil to Lymphocyte Ratio (NLR)
meta_imputed$NLR <- meta_imputed$Neutrophils / meta_imputed$Lymphocytes

# Optionally, write to CSV
# write.csv(meta_imputed, file = "metadata_imputed.csv", row.names = FALSE)


## 03. Data cleaning and harmonization ----

# Keep only data from common samples
common_samples <- Reduce(intersect, list(colnames(counts_matrix_RNAseq)[-1], colnames(counts_matrix_miRNA)[-1], meta_imputed$Unique_ID))

counts_matrix_RNAseq_filtered <- as.data.frame(counts_matrix_RNAseq[, c("Geneid", common_samples)])
counts_matrix_miRNA_filtered <- as.data.frame(counts_matrix_miRNA[, c("miRNA", common_samples)])
metadata_filtered <- as.data.frame(meta_imputed) %>% dplyr::filter(Unique_ID %in% common_samples)

#Generate a subset of the metadata_filtered table with only the columns of interest for the analysis
metadata_useful <- metadata_filtered[, c("Unique_ID", "PATNO", "EVENT_ID", "NSD_STAGE", "New_Stage", "COHORT", "SEX", "age", "WBC", "Platelets", "Neutrophils", "Lymphocytes", "BMI", "NLR")]


#Write cleaned data to CSV
write.csv(metadata_filtered, file = "metadata_filtered.csv", row.names = FALSE)
write.csv(metadata_useful, file = "metadata_useful.csv", row.names = FALSE)
fwrite(as.data.frame(counts_matrix_RNAseq_filtered), file = "counts_matrix_RNAseq_filtered.csv")
fwrite(as.data.frame(counts_matrix_miRNA_filtered), file = "counts_matrix_miRNA_filtered.csv")


## 04. Genetic data preprocessing ----

#Round the counts values of miRNA featurecounts results downward

numeric_cols <- sapply(counts_matrix_miRNA_filtered, is.numeric)
counts_matrix_miRNA_filtered[numeric_cols] <- floor(counts_matrix_miRNA_filtered[numeric_cols])

# Prepare count matrix for processing (use RNA-seq data)
counts_data <- counts_matrix_RNAseq_filtered[, -1]  # Remove first column (Geneid)

# Ensure all columns are numeric
counts_data[] <- lapply(counts_data, function(x) {
  if (is.character(x) || is.factor(x)) {
    as.numeric(as.character(x))
  } else {
    as.numeric(x)
  }
})

# Convert to matrix and set Geneid as rownames
counts <- as.matrix(counts_data)
rownames(counts) <- counts_matrix_RNAseq_filtered$Geneid

# Extract Ensembl gene IDs (remove version suffix)
ensg <- sub("\\..*$", "", rownames(counts))

# Retrieve gene's mapping symbols and biotype from EnsDb
ah <- AnnotationHub()
q <- query(ah, c("EnsDb", "Homo sapiens", "GRCh38"))
stopifnot(length(q) > 0)
edb <- q[[length(q)]]

gene_annotations <- as.data.frame(genes(edb, return.type = "DataFrame")[, c("gene_id", "gene_name", "gene_biotype")])
colnames(gene_annotations) <- c("ensembl_gene_id", "gene_name", "gene_biotype")

# Align with your matrix
idx <- match(ensg, gene_annotations$ensembl_gene_id)
sym_raw <- gene_annotations$gene_name[idx]
biotype <- gene_annotations$gene_biotype[idx]

# Correct gene symbols with HGNChelper (aliases HGNC approved symbols)
check <- HGNChelper::checkGeneSymbols(sym_raw, unmapped.as.na = TRUE)
sym_clean <- check$Suggested.Symbol
sym_clean[!nzchar(sym_clean)] <- NA  # empty strings → NA

n_corrected <- sum(check$x != check$Suggested.Symbol & !is.na(check$Suggested.Symbol), na.rm = TRUE)

# Build a working table and calculate mean expression
mean_expr <- rowMeans(counts, na.rm = TRUE)
tab <- data.frame(
  ensg = ensg,
  symbol = sym_clean,
  biotype = biotype,
  mean_expr = mean_expr,
  stringsAsFactors = FALSE
)

# Discard rows without symbol
tab$keep_map <- !is.na(tab$symbol)
tab2 <- tab[tab$keep_map, ]
counts_filtered <- counts[tab$keep_map, , drop = FALSE]

# Resolve duplicate gene symbols by keeping highest expression
o <- order(tab2$symbol, -tab2$mean_expr)
tab2_ord <- tab2[o, ]
counts_filtered_ordered <- counts_filtered[o, , drop = FALSE]

keep_top <- !duplicated(tab2_ord$symbol) 
tab_final <- tab2_ord[keep_top, ]
counts_symbol <- counts_filtered_ordered[keep_top, , drop = FALSE]
rownames(counts_symbol) <- tab_final$symbol


# Create a gene annotation mapping table
mapping_table <- tab_final[, c("ensg", "symbol", "biotype", "mean_expr")] # mapping reference
#Save mapping table to CSV
write.csv(mapping_table, file = "gene_annotation_mapping_table.csv", row.names = FALSE)

# Convert back to data.frame format for consistency with pipeline
counts_matrix_RNAseq_symbols <- data.frame(
  Geneid = rownames(counts_symbol),
  counts_symbol,
  stringsAsFactors = FALSE
)

# Save count matrix to CSV
fwrite(counts_matrix_RNAseq_symbols, file = "counts_matrix_RNAseq_symbols.csv")

## 05. Differential Gene Expression Analysis ----

# Convert count matrices to matrix format for MIRit
rownames(counts_matrix_miRNA_filtered) <- counts_matrix_miRNA_filtered$miRNA

counts_matrix_miRNA_filtered <- counts_matrix_miRNA_filtered[, -1]
counts_matrix_miRNA_filtered <- as.matrix(counts_matrix_miRNA_filtered)

counts_matrix_RNAseq_symbols <- as.matrix(counts_matrix_RNAseq_symbols[, -1])
rownames(counts_matrix_RNAseq_symbols) <- counts_matrix_RNAseq_symbols$Geneid
colnames(counts_matrix_RNAseq_symbols) <- sub("^X", "", colnames(counts_matrix_RNAseq_symbols))

# Define sample metadata for MIRit
meta <- data.frame(
  primary = colnames(counts_matrix_miRNA_filtered),
  mirnaCol = colnames(counts_matrix_miRNA_filtered),
  geneCol = colnames(counts_matrix_RNAseq_symbols),
  disease = gsub("-", "_", metadata_useful$New_Stage[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)]),
  age = metadata_useful$age[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)],
  SEX = as.factor(metadata_useful$SEX[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)]),
  BMI = metadata_useful$BMI[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)],
  WBC = metadata_useful$WBC[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)],
  Platelets = metadata_useful$Platelets[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)],
  NLR = metadata_useful$NLR[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)],
  patient = metadata_useful$Unique_ID[match(colnames(counts_matrix_miRNA_filtered), metadata_useful$Unique_ID)]
)

# Create the MirnaExperiment object
experiment <- MirnaExperiment(
  mirnaExpr = counts_matrix_miRNA_filtered,
  geneExpr = counts_matrix_RNAseq_symbols,
  samplesMetadata = meta,
  pairedSamples = TRUE
)

# Visualize expression variability
geneMDS <- plotDimensions(experiment, assay = "genes", condition = "SEX", title = "MDS plot for mRNA")
mirnaMDS <- plotDimensions(experiment, assay = "microRNA", condition = "SEX", title = "MDS plot for miRNA")

plotmds <- ggpubr::ggarrange(geneMDS, mirnaMDS,
                                  nrow = 1, labels = "AUTO", common.legend = TRUE)

ggsave("MDS_plot_mRNA&miRNA.png", plot = plotmds, width = 15, height = 6, dpi = 300)
ggsave("MDS_plot_mRNA&miRNA.svg", plot = plotmds, width = 15, height = 6)

# Design the linear model for both miRNA and mRNA
model <- ~ disease + age + SEX + BMI + WBC + Platelets + NLR

# Run Exact tests for each comparison of interest, separately for mRNA and miRNA
run_comparison <- function(expt, cond, design, ...,
                           pCutoff = 0.05, logFC = 0.58,
                           method = "edgeR",
                           filterByExpr.args = list(min.count = 5,
                                                    min.total.count = 10,
                                                    large.n = 10,
                                                    min.prop = 0.3)) {
  contrast <- paste0(cond, "-Healthy")
  list(
    genes = performGeneDE(expt,
                          group   = "disease",
                          contrast= contrast,
                          design  = design,
                          pCutoff = pCutoff,
                          logFC   = logFC,
                          method  = method,
                          filterByExpr.args = filterByExpr.args,
                          ...),
    mirna = performMirnaDE(expt,
                           group   = "disease",
                           contrast= contrast,
                           design  = design,
                           pCutoff = pCutoff,
                           logFC   = logFC,
                           method  = method,
                           filterByExpr.args = filterByExpr.args,
                           ...)
  )
}

conds <- c("Pre_clinical", "Mild", "Moderate", "Severe")


de_results <- setNames(lapply(conds, run_comparison,
                              expt   = experiment,
                              design = model),
                       conds)

# extract the gene/miRNA tables in one go
deGenes  <- lapply(de_results, function(x) geneDE(x$genes))
deMirnas <- lapply(de_results, function(x) mirnaDE(x$mirna))


# extract the DEG/DEmiRNA tables
deGenes  <- lapply(de_results, function(x) geneDE(x$genes))
deMirnas <- lapply(de_results, function(x) mirnaDE(x$mirna))

# Save DE results to CSV files
for (cond in names(deGenes)) {
  fwrite(as.data.frame(deGenes[[cond]]),
         file = sprintf("DE_genes_%s_vs_Healthy.csv", cond))
  fwrite(as.data.frame(deMirnas[[cond]]),
         file = sprintf("DE_miRNAs_%s_vs_Healthy.csv", cond))
}

# Save the DGE results (not only the DEGs, but all the genes/miRNAs tested)
dgeGenes_all  <- lapply(de_results,
                        function(res) geneDE(res$genes,
                                             onlySignificant = FALSE))
dgeMirnas_all <- lapply(de_results,
                        function(res) mirnaDE(res$mirna,
                                              onlySignificant = FALSE))

names(dgeGenes_all)  <- names(dgeMirnas_all) <- names(de_results)


for (cond in names(dgeGenes_all)) {
  fwrite(as.data.frame(dgeGenes_all[[cond]]),
         file = sprintf("DGE_results_all_genes_%s_vs_Healthy.csv",
                        cond))
  fwrite(as.data.frame(dgeMirnas_all[[cond]]),
         file = sprintf("DGE_results_all_miRNAs_%s_vs_Healthy.csv",
                        cond))
}

#Make volcano plots for genes, miRNAs and conditions
conditions <- c("Preclinical", "Mild", "Moderate", "Severe")

logfc_thr <- 0.58
fc_thr    <- 2^logfc_thr        # ~1.5
padj_thr  <- 0.05


read_deg <- function(file, condition, omics) {
  df <- read_csv(file, show_col_types = FALSE)

  needed <- c("ID", "logFC", "adj.P.Val")
  miss <- setdiff(needed, names(df))
  if (length(miss) > 0) stop(sprintf("Missing columns in %s: %s", file, paste(miss, collapse = ", ")))

  df %>%
    transmute(
      ID = as.character(ID),
      logFC = as.numeric(logFC),
      adj.P.Val = as.numeric(`adj.P.Val`),
      condition = condition,
      omics = omics
    ) %>%
    dplyr::filter(!is.na(logFC), !is.na(adj.P.Val))
}

add_stats <- function(df) {
  df %>%
    mutate(
      FC = logFC,
      neglog10_padj = -log10(adj.P.Val),
      category = case_when(
        adj.P.Val < padj_thr & logFC >  logfc_thr ~ "Up",
        adj.P.Val < padj_thr & logFC < -logfc_thr ~ "Down",
        TRUE ~ "NotSig"
      )
    )
}

volcano <- function(df, title = NULL) {
  ggplot(df, aes(x = FC, y = neglog10_padj, color = category)) +
    geom_point(alpha = 0.65, size = 1.1) +
    geom_vline(xintercept = c(-logfc_thr, logfc_thr), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", linewidth = 0.4) +
    scale_color_manual(values = c(NotSig = "grey70", Down = "blue", Up = "red")) +
    labs(
      title = title,
      x = expression(log[2]("Fold change")),
      y = expression(-log[10]("adj.P.Val"))
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

mrna_files  <- paste0("DGE_results_genes",  conditions, "_vs_Healthy.csv")
mirna_files <- paste0("DGE_results_miRNAs", conditions, "_vs_Healthy.csv")


mrna_plots <- lapply(seq_along(conditions), function(i) {
  read_deg(mrna_files[i], conditions[i], "mRNA") %>%
    add_stats() %>%
    volcano(title = conditions[i])
})

mirna_plots <- lapply(seq_along(conditions), function(i) {
  read_deg(mirna_files[i], conditions[i], "miRNA") %>%
    add_stats() %>%
    volcano(title = conditions[i])
})

grid_plot <- (wrap_plots(mrna_plots, nrow = 1) / wrap_plots(mirna_plots, nrow = 1))

grid_plot

# Save the volcano grid plot
ggsave("volcano_grid_mrna_mirna.png", grid_plot, width = 14, height = 6.5, dpi = 300)
ggsave("volcano_grid_mrna_mirna.svg", grid_plot, width = 14, height = 6.5)


# Make an Upset plots of overlapping DEGs between stages

degs_list <- setNames(lapply(conds, function(cond) rownames(deGenes[[cond]])), conds)

names(degs_list) <- sub("_", "-", names(degs_list))

upset_data <- fromList(degs_list)
upset_plot <- upset(upset_data,
                    nsets = 4,
                    nintersects = NA,
                    order.by = "freq",
                    main.bar.color = "steelblue",
                    sets.bar.color = "skyblue",
                    text.scale = c(2, 2, 2, 1.5, 2, 2))
upset_plot

# Save UpSet plot
png("UpSet_plot_Genes_DEGs.png", width = 10, height = 6, units = "in", res = 300)
print(upset_plot)
dev.off()

# Save in SVG
svg("UpSet_plot_Genes_DEGs.svg", width = 10, height = 6)
print(upset_plot)
dev.off()

## 06. miRNA-mRNA integration ----

# Retrieve targets for all miRNA experiments
mirna_experiments <- lapply(de_results, function(x) getTargets(x$mirna))
names(mirna_experiments) <- names(de_results)

# Save miRNA target pairs to CSV files
for (cond in names(mirna_experiments)) {
  fwrite(as.data.frame(mirnaTargets(mirna_experiments[[cond]])),
         file = sprintf("miRNA_targets_%s_vs_Healthy.csv", cond))
}

# Register parallel backend
register(MulticoreParam(workers = 16))

# Batch correction
batch_params <- list(
  assays = c("microRNA", "genes"),
  batch = "patient",
  batch2 = "SEX",
  covariates = c("age", "BMI", "WBC", "Platelets", "NLR"),
  includeWsva = TRUE,
  n.sv = 2,
  weight.by.sd = TRUE
)

apply_batch_correction <- function(joint_exp, params) {
  for (assay in params$assays) {
    joint_exp <- batchCorrection(
      joint_exp,
      assay = assay,
      batch = params$batch,
      batch2 = params$batch2,
      covariates = params$covariates,
      includeWsva = params$includeWsva,
      n.sv = params$n.sv,
      weight.by.sd = params$weight.by.sd
    )
  }
  joint_exp
}

run_integration_pipeline <- function(mirna_exp, gene_de, cond, bpparam) {

  joint_exp <- mirna_exp
  slot(joint_exp, "geneDE") <- gene_de
  
  joint_exp <- apply_batch_correction(joint_exp, batch_params)
  
  joint_exp <- mirnaIntegration(joint_exp, test = "correlation", BPPARAM = bpparam)
  
  integration(joint_exp)
}

# Run integration
integration_results <- mapply(
  function(cond) {
    
    bpparam <- MulticoreParam(workers = 16)
    
    run_integration_pipeline(
      mirna_exp = de_results[[cond]]$mirna,
      gene_de = slot(de_results[[cond]]$genes, "geneDE"),
      cond = cond,
      bpparam = bpparam
    )
  },
  cond = conds,
)

names(integration_results) <- conds

# Save integration results to CSV files
for (cond in conds) {
  fwrite(as.data.frame(integration_results[[cond]]),
         file = sprintf("miRNA_mRNA_correlation_%s_vs_Healthy.csv", cond))
}


## 07. Machine Learning Model for Pre-clinical Stage Prediction ----
library(caret)
library(randomForest)
library(e1071)
library(glmnet)
library(gbm)
library(class)
library(klaR)
library(MASS)
library(kernlab)
library(nnet)
library(pROC)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggplot2)

# Setup parallel processing
cl <- makeCluster(16)
registerDoParallel(cl)

# Selected miRNAs for prediction (DEGs from preclinical vs healthy)
selected_mirnas <- rownames(deMirnas[[conds[1]]])[which(deMirnas[[conds[1]]]$adj.P.Val < 0.05 & abs(deMirnas[[conds[1]]]$logFC) > 0.58)]

cat("Number of selected miRNA features:", length(selected_mirnas), "\n\n")

# Extract raw counts for selected miRNAs and transpose to have samples as rows and miRNAs as columns
mirna_counts_selected <- as.data.frame(t(counts_matrix_miRNA_filtered[selected_mirnas, ]))
colnames(mirna_counts_selected) <- selected_mirnas
mirna_counts_selected$Unique_ID <- rownames(mirna_counts_selected)

# Merge with metadata
ml_dataset <- merge(mirna_counts_selected, metadata_useful, by = "Unique_ID")

# Select relevant columns
ml_dataset <- ml_dataset[, c("Unique_ID", selected_mirnas, "New_Stage")]

# Only keep healthy and preclinical samples for binary classification
ml_dataset <- ml_dataset %>% dplyr::filter(New_Stage %in% c("Healthy", "Pre-clinical"))

# Turn New_Stage into a factor
ml_dataset$New_Stage <- as.factor(ml_dataset$New_Stage)

# Rename New_Stage variable into 'Class'
colnames(ml_dataset)[which(names(ml_dataset) == "New_Stage")] <- "Class"

# Change Pre-clinical to Preclinical (required for caret - no hyphens)
levels(ml_dataset$Class) <- gsub("Pre-clinical", "Preclinical", levels(ml_dataset$Class))

# Ensure Healthy is the reference level
ml_dataset$Class <- relevel(ml_dataset$Class, ref = "Healthy")

# Set Unique_ID as row names and remove it from the dataset
rownames(ml_dataset) <- ml_dataset$Unique_ID
ml_dataset$Unique_ID <- NULL

cat("Total samples:", nrow(ml_dataset), "\n")
cat("Class distribution:\n")
print(table(ml_dataset$Class))
cat("\n")

# Check for and remove near-zero variance predictors
nzv <- nearZeroVar(ml_dataset[, -ncol(ml_dataset)], saveMetrics = TRUE)
if(sum(nzv$nzv) > 0) {
  cat("Removing", sum(nzv$nzv), "near-zero variance features\n\n")
  ml_dataset <- ml_dataset[, c(rownames(nzv)[!nzv$nzv], "Class")]
}

# Split the dataset into training (70%) and testing (30%) sets stratified by Class
set.seed(123)
trainIndex <- createDataPartition(ml_dataset$Class, p = 0.7, list = FALSE, times = 1)
ml_train <- ml_dataset[trainIndex, ]
ml_test <- ml_dataset[-trainIndex, ]

cat("Training set size:", nrow(ml_train), "\n")
cat("Training set distribution:\n")
print(table(ml_train$Class))
cat("\nTest set size:", nrow(ml_test), "\n")
cat("Test set distribution:\n")
print(table(ml_test$Class))
cat("\n")

# Enhanced train control with repeated cross-validation
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  allowParallel = TRUE,
  verboseIter = FALSE,
  returnResamp = "all"
)

# Initialize list to store all models
model_list <- list()
model_names <- c("Random Forest", "SVM (Radial)", "GLMNet", "Gradient Boosting", "k-NN", "Naive Bayes", "Linear Discriminant Analysis", "Neural Network")

# Model training

# a. Random Forest
cat("Training Random Forest...\n")
set.seed(123)
model_list$rf <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "rf",
  trControl = train_control,
  metric = "ROC",
  tuneGrid = expand.grid(mtry = c(2, 4, 6, 8, 10)),
  importance = TRUE,
  ntree = 500
)

# b. SVM with Radial Kernel
cat("Training SVM (Radial Kernel)...\n")
set.seed(123)
model_list$svm <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "svmRadial",
  trControl = train_control,
  metric = "ROC",
  tuneLength = 10,
  preProcess = c("center", "scale")
)

# c. GLMNet (Elastic Net)
cat("Training GLMNet (Elastic Net)...\n")
set.seed(123)
model_list$glmnet <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "glmnet",
  trControl = train_control,
  metric = "ROC",
  tuneGrid = expand.grid(
    alpha = seq(0, 1, by = 0.1),
    lambda = exp(seq(-6, 1, length = 20))
  ),
  preProcess = c("center", "scale")
)

# d. Gradient Boosting Machine
cat("Training Gradient Boosting Machine...\n")
set.seed(123)
model_list$gbm <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "gbm",
  trControl = train_control,
  metric = "ROC",
  tuneGrid = expand.grid(
    n.trees = c(100, 200, 300),
    interaction.depth = c(1, 3, 5),
    shrinkage = c(0.01, 0.1),
    n.minobsinnode = 10
  ),
  verbose = FALSE
)

# e. k-Nearest Neighbors
cat("Training k-Nearest Neighbors...\n")
set.seed(123)
model_list$knn <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "knn",
  trControl = train_control,
  metric = "ROC",
  tuneGrid = expand.grid(k = seq(3, 21, by = 2)),
  preProcess = c("center", "scale")
)

# f. Naive Bayes
cat("Training Naive Bayes...\n")
set.seed(123)
model_list$nb <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "naive_bayes",
  trControl = train_control,
  metric = "ROC",
  tuneGrid = expand.grid(
    laplace = c(0, 0.5, 1),
    usekernel = c(TRUE, FALSE),
    adjust = c(0.5, 1, 1.5)
  )
)

# g. Linear Discriminant Analysis
cat("Training Linear Discriminant Analysis...\n")
set.seed(123)
model_list$lda <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "lda",
  trControl = train_control,
  metric = "ROC",
  preProcess = c("center", "scale")
)

# h. Neural Network
cat("Training Neural Network...\n")
set.seed(123)
model_list$nnet <- caret::train(
  Class ~ .,
  data = ml_train,
  method = "nnet",
  trControl = train_control,
  metric = "ROC",
  tuneGrid = expand.grid(
    size = c(5, 10, 15),
    decay = c(0, 0.01, 0.1)
  ),
  preProcess = c("center", "scale"),
  trace = FALSE,
  MaxNWts = 5000,
  maxit = 200
)

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

# Evaluate models on the test set

# Initialize results storage
predictions_list <- list()
probability_list <- list()
confusion_matrices <- list()
roc_objects <- list()

# Make predictions for all models
for(i in seq_along(model_list)) {
  model_name <- names(model_list)[i]
  model <- model_list[[i]]
  
  cat("Evaluating", model_names[i], "...\n")
  
  # Predictions
  predictions_list[[model_name]] <- predict(model, ml_test)
  probability_list[[model_name]] <- predict(model, ml_test, type = "prob")
  
  # Confusion matrix
  confusion_matrices[[model_name]] <- confusionMatrix(
    predictions_list[[model_name]], 
    ml_test$Class,
    positive = "Preclinical"
  )
  
  # ROC curve
  roc_objects[[model_name]] <- roc(
    ml_test$Class,
    probability_list[[model_name]]$Preclinical,
    levels = levels(ml_test$Class),
    direction = "<"
  )
}

#Obtain test set results

# Print confusion matrices
for(i in seq_along(confusion_matrices)) {
  cat("\n--- ", model_names[i], " Confusion Matrix ---\n")
  print(confusion_matrices[[i]])
  cat("\n")
}

# Extract comprehensive metrics
extract_comprehensive_metrics <- function(cm, roc_obj, model_name) {
  by_class <- cm$byClass
  overall <- cm$overall
  
  data.frame(
    Model = model_name,
    Accuracy = overall['Accuracy'],
    Kappa = overall['Kappa'],
    Sensitivity = by_class['Sensitivity'],
    Specificity = by_class['Specificity'],
    PPV = by_class['Pos Pred Value'],
    NPV = by_class['Neg Pred Value'],
    Precision = by_class['Precision'],
    Recall = by_class['Recall'],
    F1 = by_class['F1'],
    Balanced_Accuracy = by_class['Balanced Accuracy'],
    AUC = as.numeric(auc(roc_obj)),
    row.names = NULL
  )
}

# Compile all metrics
performance_comparison <- do.call(rbind, lapply(seq_along(confusion_matrices), function(i) {
  extract_comprehensive_metrics(
    confusion_matrices[[i]], 
    roc_objects[[i]], 
    model_names[i]
  )
}))

# Save test set model performance comparison
fwrite(performance_comparison, file = "ML_Model_Performance_Comparison.csv")

# Save individual confusion matrices
for(i in seq_along(confusion_matrices)) {
  model_name_clean <- gsub(" ", "_", gsub("[()]", "", model_names[i]))
  fwrite(
    as.data.frame(confusion_matrices[[i]]$table), 
    file = paste0("ML_CM_", model_name_clean, "_Test.csv")
  )
}

# Create comprehensive ROC curve plot
cat("Generating ROC curves...\n")

roc_plot_data <- data.frame()
for(i in seq_along(roc_objects)) {
  temp_df <- data.frame(
    Sensitivity = roc_objects[[i]]$sensitivities,
    Specificity = roc_objects[[i]]$specificities,
    Model = model_names[i],
    AUC = as.numeric(auc(roc_objects[[i]]))
  )
  roc_plot_data <- rbind(roc_plot_data, temp_df)
}

# Create ROC plot
roc_plot <- ggplot(roc_plot_data, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_line(size = 0.5) +
  geom_abline(linetype = "dashed", color = "gray50") +
  labs(
    title = "ROC Curves - All Models (Test Set)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Model (AUC)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  scale_color_brewer(palette = "Set1") +
  coord_equal()

# Add AUC to legend labels
roc_plot <- roc_plot + 
  scale_color_manual(
    values = rainbow(length(model_names)),
    labels = paste0(model_names, " (", sprintf("%.3f", performance_comparison$AUC), ")")
  )
roc_plot
ggsave("ML_ROC_Curves_All_Models_Test.png", plot = roc_plot, width = 12, height = 8, dpi = 300)
ggsave("ML_ROC_Curves_All_Models_Test.pdf", plot = roc_plot, width = 12, height = 8)
ggsave("ML_ROC_Curves_All_Models_Test.svg", plot = roc_plot, width = 12, height = 8)

# Feature importance for applicable models
cat("\nExtracting feature importance...\n")

varImp_rf <- varImp(model_list$rf, scale = FALSE)
varImp_gbm <- varImp(model_list$gbm, scale = FALSE)
varImp_glmnet <- varImp(model_list$glmnet, scale = FALSE)

# Save variable importance tables
fwrite(as.data.frame(varImp_rf$importance), file = "ML_Variable_Importance_RandomForest.csv")
fwrite(as.data.frame(varImp_gbm$importance), file = "ML_Variable_Importance_GBM.csv")
fwrite(as.data.frame(varImp_glmnet$importance), file = "ML_Variable_Importance_GLMNet.csv")

# Overfitting diagnostics: Compare train vs test performance
compute_train_metrics <- function(model, train_df, positive = "Preclinical") {
  pred <- predict(model, train_df)
  prob <- predict(model, train_df, type = "prob")[, positive]

  cm <- confusionMatrix(pred, train_df$Class, positive = positive)

  roc_obj <- roc(
    response = train_df$Class,
    predictor = prob,
    levels = levels(train_df$Class),
    direction = "<",
    quiet = TRUE
  )

  data.frame(
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Kappa = as.numeric(cm$overall["Kappa"]),
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    F1 = as.numeric(cm$byClass["F1"]),
    Balanced_Accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    AUC = as.numeric(auc(roc_obj))
  )
}

extract_cv_best_metrics <- function(model, positive_metric = "ROC") {
  bt <- model$bestTune
  res <- model$results

  tune_cols <- intersect(names(bt), names(res))
  if (length(tune_cols) > 0) {
    for (tc in tune_cols) res <- res[res[[tc]] == bt[[tc]], , drop = FALSE]
  }

  out <- res[1, , drop = FALSE]

  data.frame(
    CV_ROC = if ("ROC" %in% names(out)) as.numeric(out$ROC) else NA_real_,
    CV_Sens = if ("Sens" %in% names(out)) as.numeric(out$Sens) else NA_real_,
    CV_Spec = if ("Spec" %in% names(out)) as.numeric(out$Spec) else NA_real_
  )
}

extract_cv_resample_distribution <- function(model, model_label) {
  rp <- model$pred
  if (is.null(rp) || nrow(rp) == 0) return(NULL)

  bt <- model$bestTune
  tune_cols <- intersect(names(bt), names(rp))
  if (length(tune_cols) > 0) {
    for (tc in tune_cols) rp <- rp[rp[[tc]] == bt[[tc]], , drop = FALSE]
  }

  if (!("Preclinical" %in% names(rp))) return(NULL)

  roc_by_resample <- rp %>%
    group_by(Resample) %>%
    summarise(
      ROC = as.numeric(auc(roc(obs, Preclinical, levels = levels(obs), direction = "<", quiet = TRUE))),
      .groups = "drop"
    ) %>%
    mutate(Model = model_label)

  roc_by_resample
}

# Build Train metrics for each model
model_keys <- names(model_list)

train_metrics <- do.call(rbind, lapply(model_keys, function(k) {
  m <- model_list[[k]]
  tm <- compute_train_metrics(m, ml_train, positive = "Preclinical")
  tm$ModelKey <- k
  tm
}))

key_to_name <- data.frame(
  ModelKey = model_keys,
  Model = model_names[seq_along(model_keys)],
  stringsAsFactors = FALSE
)

train_metrics <- left_join(train_metrics, key_to_name, by = "ModelKey")

# Extract CV metrics at best tuning
cv_best <- do.call(rbind, lapply(model_keys, function(k) {
  m <- model_list[[k]]
  cvm <- extract_cv_best_metrics(m)
  cvm$ModelKey <- k
  cvm
}))

cv_best <- left_join(cv_best, key_to_name, by = "ModelKey")

# Merge with test metrics
test_metrics <- performance_comparison %>%
  transmute(
    Model = as.character(Model),
    Test_Accuracy = as.numeric(Accuracy),
    Test_F1 = as.numeric(F1),
    Test_Balanced_Accuracy = as.numeric(Balanced_Accuracy),
    Test_AUC = as.numeric(AUC)
  )

train_metrics2 <- train_metrics %>%
  transmute(
    Model,
    Train_Accuracy = Accuracy,
    Train_F1 = F1,
    Train_Balanced_Accuracy = Balanced_Accuracy,
    Train_AUC = AUC
  )

cv_best2 <- cv_best %>%
  transmute(
    Model,
    CV_AUC = CV_ROC
  )

overfit_table <- train_metrics2 %>%
  left_join(cv_best2, by = "Model") %>%
  left_join(test_metrics, by = "Model") %>%
  mutate(
    Gap_AUC_Train_CV = Train_AUC - CV_AUC,
    Gap_AUC_Train_Test = Train_AUC - Test_AUC,
    Gap_AUC_CV_Test = CV_AUC - Test_AUC,
    Gap_F1_Train_Test = Train_F1 - Test_F1,
    Gap_BAcc_Train_Test = Train_Balanced_Accuracy - Test_Balanced_Accuracy
  ) %>%
  arrange(desc(Gap_AUC_Train_Test))

print(overfit_table)

fwrite(overfit_table, file = "ML_Overfitting_Diagnostics_Table.csv")


# Plot Overfitting gap barplot (Train - Test)
p_gap <- overfit_table %>%
  mutate(Model = factor(Model, levels = rev(unique(Model)))) %>%
  ggplot(aes(x = Model, y = Gap_AUC_Train_Test)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Overfitting indicator (AUC gap: Train - Test)", x = NULL, y = "AUC gap")

p_gap
ggsave("ML_Overfitting_Gap_AUC_TrainMinusTest.png", p_gap, width = 10, height = 6, dpi = 300)
ggsave("ML_Overfitting_Gap_AUC_TrainMinusTest.pdf", p_gap, width = 10, height = 6)

# Save final workspace
save.image(file = "PD_BT_analysis.RData")
