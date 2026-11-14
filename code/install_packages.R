# Install CRAN and Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
cran <- c("ggplot2", "dplyr", "tidyr", "readr", "pheatmap",
          "rmarkdown", "knitr", "RColorBrewer")
bioc <- c("GEOquery", "Biobase", "affy", "limma", "hgu133plus2.db")
for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}
message("Package installation finished.")
