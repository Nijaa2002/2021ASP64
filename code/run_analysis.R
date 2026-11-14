# Complete analysis pipeline for ICA3 Bioinformatics 2026
# Run from the project root: source("code/run_analysis.R")

options(stringsAsFactors = FALSE)
set.seed(2026)

required <- c("GEOquery", "Biobase", "affy", "limma", "ggplot2",
              "dplyr", "pheatmap", "hgu133plus2.db", "AnnotationDbi")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Run source('code/install_packages.R') first.")
}

dirs <- c("data/raw_CEL", "data/processed", "figures", "results")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# 1. Validate and process the nine supplied raw CEL files
cel_files <- list.files("data/raw_CEL", pattern = "\\.CEL\\.gz$", full.names = TRUE,
                        ignore.case = TRUE)
if (length(cel_files) != 9L) {
  stop("Expected 9 supplied .CEL.gz files in data/raw_CEL; found ", length(cel_files), ".")
}
if (any(file.info(cel_files)$size <= 0)) stop("At least one CEL file is empty.")

message("Importing supplied CEL files and performing RMA preprocessing...")
# affy expects uncompressed CEL files on some installations. Decompress to a
# temporary directory; the original .CEL.gz files remain unchanged in the repo.
cel_tmp <- tempfile("cel_uncompressed_")
dir.create(cel_tmp)
cel_uncompressed <- file.path(cel_tmp, sub("\\.gz$", "", basename(cel_files),
                                           ignore.case = TRUE))
for (i in seq_along(cel_files)) {
  con_in <- gzfile(cel_files[i], open = "rb")
  con_out <- file(cel_uncompressed[i], open = "wb")
  repeat {
    chunk <- readBin(con_in, what = "raw", n = 1024 * 1024)
    if (!length(chunk)) break
    writeBin(chunk, con_out)
  }
  close(con_in)
  close(con_out)
}
raw_affy <- affy::ReadAffy(filenames = cel_uncompressed)
eset_cel <- affy::rma(raw_affy)
cel_expr <- Biobase::exprs(eset_cel)
write.csv(cel_expr, "data/processed/uploaded_nine_RMA_expression.csv")
write.csv(data.frame(sample = colnames(cel_expr),
                     file = basename(cel_files)),
          "results/uploaded_CEL_sample_manifest.csv", row.names = FALSE)

png("figures/uploaded_CEL_boxplot.png", width = 1400, height = 900, res = 150)
boxplot(cel_expr, las = 2, outline = FALSE,
        main = "RMA-normalised expression: nine supplied adjacent-normal samples",
        ylab = "Log2 expression")
dev.off()

cel_cor <- cor(cel_expr, use = "pairwise.complete.obs")
pheatmap::pheatmap(cel_cor, filename = "figures/uploaded_CEL_sample_correlation.png",
                   main = "Correlation between supplied CEL samples")
pc_cel <- prcomp(t(cel_expr), scale. = FALSE)
pca_cel <- data.frame(sample = rownames(pc_cel$x),
                      PC1 = pc_cel$x[, 1], PC2 = pc_cel$x[, 2])
write.csv(pca_cel, "results/uploaded_CEL_PCA.csv", row.names = FALSE)
p <- ggplot2::ggplot(pca_cel, ggplot2::aes(PC1, PC2, label = sample)) +
  ggplot2::geom_point(size = 3) + ggplot2::geom_text(nudge_y = 0.5, size = 3) +
  ggplot2::labs(title = "PCA of the nine supplied CEL samples",
                subtitle = "All nine supplied samples are GEO-labelled adjacent normal lung",
                x = "PC1", y = "PC2") + ggplot2::theme_minimal()
ggplot2::ggsave("figures/uploaded_CEL_PCA.png", p, width = 8, height = 6, dpi = 150)

# 2. Retrieve parent public GEO series matrix and phenotype data
message("Retrieving GSE19804 from NCBI GEO. Internet access is required...")
geo <- GEOquery::getGEO("GSE19804", GSEMatrix = TRUE, getGPL = FALSE)
if (length(geo) < 1L) stop("GEOquery returned no Series Matrix for GSE19804.")
# Select the GPL570 expression set if multiple matrices are returned.
platforms <- vapply(geo, function(x) as.character(Biobase::annotation(x)), character(1))
idx <- which(grepl("GPL570", platforms, ignore.case = TRUE))
eset <- if (length(idx)) geo[[idx[1]]] else geo[[1]]
expr <- Biobase::exprs(eset)
pheno <- Biobase::pData(eset)
if (is.null(colnames(expr)) || ncol(expr) != nrow(pheno)) {
  stop("Expression matrix and phenotype table are not aligned.")
}
if (!all(colnames(expr) == rownames(pheno))) {
  pheno <- pheno[colnames(expr), , drop = FALSE]
}
if (any(!is.finite(expr))) stop("Expression matrix contains missing or non-finite values.")

# 3. Identify tissue groups conservatively from GEO sample annotations
annotation_text <- apply(pheno, 1, function(x) paste(x, collapse = " | "))
annotation_text <- tolower(annotation_text)
group <- ifelse(grepl("lung cancer|cancer tissue|tumou?r tissue|lung tumou?r", annotation_text),
                "Tumour",
                ifelse(grepl("lung normal|adjacent normal|normal adjacent|normal tissue", annotation_text),
                       "Adjacent_Normal", NA_character_))
keep <- !is.na(group)
expr <- expr[, keep, drop = FALSE]
pheno <- pheno[keep, , drop = FALSE]
group <- factor(group[keep], levels = c("Adjacent_Normal", "Tumour"))
if (length(unique(group)) != 2L) {
  stop("Could not safely identify both tumour and adjacent-normal groups from GSE19804 metadata. ",
       "Inspect the GEO phenotype table before proceeding; no group labels have been invented.")
}
if (any(table(group) < 2L)) stop("At least one group has fewer than two samples.")

group_counts <- as.data.frame(table(group))
names(group_counts) <- c("group", "n")
write.csv(group_counts, "results/group_counts.csv", row.names = FALSE)
write.csv(data.frame(sample = rownames(pheno), group = group,
                     title = if ("title" %in% names(pheno)) pheno$title else rownames(pheno)),
          "results/analysis_sample_manifest.csv", row.names = FALSE)

# Expression distribution and PCA for full parent series comparison
png("figures/full_cohort_boxplot.png", width = 1600, height = 900, res = 150)
boxplot(expr, las = 2, outline = FALSE, col = as.numeric(group),
        main = "GSE19804 processed expression values",
        ylab = "Expression value", xlab = "Samples")
legend("topright", legend = levels(group), fill = seq_along(levels(group)), cex = 0.8)
dev.off()

pc <- prcomp(t(expr), scale. = FALSE)
pca_df <- data.frame(sample = rownames(pc$x), PC1 = pc$x[, 1],
                     PC2 = pc$x[, 2], group = group)
write.csv(pca_df, "results/full_cohort_PCA.csv", row.names = FALSE)
p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, colour = group)) +
  ggplot2::geom_point(size = 3, alpha = 0.85) +
  ggplot2::labs(title = "PCA of GSE19804 lung tissue samples",
                x = "PC1", y = "PC2", colour = "Tissue group") +
  ggplot2::theme_minimal()
ggplot2::ggsave("figures/full_cohort_PCA.png", p, width = 8, height = 6, dpi = 150)

# 4. Differential expression using limma
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- limma::lmFit(expr, design)
contrast <- limma::makeContrasts(Tumour_vs_AdjacentNormal = Tumour - Adjacent_Normal,
                                 levels = design)
fit2 <- limma::contrasts.fit(fit, contrast)
fit2 <- limma::eBayes(fit2)
de <- limma::topTable(fit2, coef = "Tumour_vs_AdjacentNormal",
                      number = Inf, adjust.method = "BH", sort.by = "P")
de$probe_id <- rownames(de)

# Add gene symbols when the platform annotation package is available
de$gene_symbol <- NA_character_
if (requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  ann <- AnnotationDbi::select(hgu133plus2.db::hgu133plus2.db,
                               keys = unique(de$probe_id),
                               keytype = "PROBEID",
                               columns = c("SYMBOL", "GENENAME"))
  ann <- ann[!duplicated(ann$PROBEID), ]
  de$gene_symbol <- ann$SYMBOL[match(de$probe_id, ann$PROBEID)]
  de$gene_name <- ann$GENENAME[match(de$probe_id, ann$PROBEID)]
}
de$significant_FDR_0_05 <- de$adj.P.Val < 0.05
write.csv(de, "results/differential_expression_all.csv", row.names = FALSE)
write.csv(utils::head(de, 50), "results/top_differentially_expressed.csv", row.names = FALSE)

# 5. Volcano plot
de$neglog10FDR <- -log10(pmax(de$adj.P.Val, .Machine$double.xmin))
de$direction <- ifelse(de$adj.P.Val < 0.05 & de$logFC > 0, "Higher_in_Tumour",
                       ifelse(de$adj.P.Val < 0.05 & de$logFC < 0,
                              "Lower_in_Tumour", "Not_FDR_significant"))
p <- ggplot2::ggplot(de, ggplot2::aes(logFC, neglog10FDR, colour = direction)) +
  ggplot2::geom_point(alpha = 0.55, size = 1.2) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
  ggplot2::labs(title = "Differential expression: tumour vs adjacent normal",
                x = "Log2 fold change", y = "-Log10 adjusted p-value",
                colour = "Result") + ggplot2::theme_minimal()
ggplot2::ggsave("figures/volcano_plot.png", p, width = 8, height = 6, dpi = 160)

# 6. Heatmap of top 30 probes by adjusted p-value
top_ids <- utils::head(de$probe_id, min(30L, nrow(de)))
heat <- expr[top_ids, , drop = FALSE]
heat_z <- t(scale(t(heat)))
heat_z[!is.finite(heat_z)] <- 0
ann_col <- data.frame(Tissue_Group = group)
rownames(ann_col) <- colnames(heat_z)
pheatmap::pheatmap(heat_z, annotation_col = ann_col,
                   filename = "figures/top30_expression_heatmap.png",
                   main = "Top 30 ranked probe sets")

# 7. Save concise summary and software environment
summary_lines <- c(
  "GSE19804 differential-expression analysis",
  paste("Samples analysed:", ncol(expr)),
  paste("Adjacent normal:", sum(group == "Adjacent_Normal")),
  paste("Tumour:", sum(group == "Tumour")),
  paste("Probe sets tested:", nrow(de)),
  paste("Probe sets with BH-adjusted p-value < 0.05:", sum(de$adj.P.Val < 0.05, na.rm = TRUE)),
  paste("GEO retrieval date:", Sys.Date()),
  "Contrast: Tumour minus Adjacent_Normal; logFC > 0 indicates higher expression in tumour."
)
writeLines(summary_lines, "results/analysis_summary.txt")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message(paste(summary_lines, collapse = "\n"))
message("Analysis complete. Render report/Bioinformatics_Report.Rmd for the report.")
