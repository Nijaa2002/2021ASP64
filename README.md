# Transcriptomic Analysis of Lung Cancer and Adjacent Normal Lung Tissue

## 1. Project overview

**Research question:** Which genes show statistically significant differences in expression between lung cancer tissue and adjacent normal lung tissue in the publicly available Taiwanese non-smoking female lung cancer cohort?

This project analyses Affymetrix Human Genome U133 Plus 2.0 microarray data from the Gene Expression Omnibus (GEO) series **GSE19804**, titled *Genome-wide screening of transcriptional modulation in non-smoking female lung cancer in Taiwan*. The uploaded archive contains nine raw CEL files (GSM494667–GSM494675). GEO identifies these samples as adjacent normal lung tissue. The uploaded files are therefore genuine raw microarray data, but they do not by themselves contain a cancer-versus-normal comparison group.

To make the biological question answerable, the reproducible workflow retains the nine uploaded CEL files and retrieves the full public GSE19804 processed expression dataset and sample metadata from NCBI GEO at run time. The differential-expression comparison uses the full series matrix when both tumour and normal groups are available. The nine supplied CEL files are also processed separately using RMA for raw-data quality control and exploratory inspection. This distinction is important: the nine supplied files must not be incorrectly labelled as tumour samples.

## 2. Abstract

Lung cancer is a major health problem, and understanding molecular changes associated with tumour tissue can support biological research into disease mechanisms and candidate biomarkers. This project investigates gene-expression differences between lung cancer and adjacent normal lung tissue using a public Affymetrix microarray cohort from non-smoking female patients in Taiwan. Nine raw CEL files supplied with the assessment are included as a traceable subset of the study; GEO metadata identify them as adjacent normal samples. Since this subset alone cannot support a tumour-versus-normal test, the analysis retrieves the parent GEO series GSE19804 to obtain the wider cohort and its sample annotations. The workflow includes data validation, RMA processing of supplied CEL files, exploratory visualisation, sample-group inspection, and differential-expression analysis with the limma framework. Results are generated from the retrieved expression matrix when the project is run. No numerical findings are asserted in advance of execution.

## 3. Introduction

Gene expression describes the activity of genes in a biological sample. Changes in expression can reflect differences in cell state, tissue composition, environmental exposure, or disease. Microarray technology measures the abundance of many transcript targets in parallel. Although RNA sequencing is now widely used, microarray datasets remain valuable because they provide accessible, well-documented public cohorts that can be reanalysed using established statistical methods.

The present project focuses on lung cancer tissue and adjacent normal lung tissue. The dataset is associated with a study of transcriptional modulation in non-smoking female lung cancer patients in Taiwan. Restricting the source cohort to the population described by the original study is important when interpreting the results: findings should not automatically be generalised to men, smokers, other populations, or every lung-cancer subtype.

The assessment archive contains nine files with GEO sample accessions GSM494667 through GSM494675. NCBI GEO describes these as “Lung Normal” samples and identifies their tissue as paired adjacent normal tissue. They were generated on the Affymetrix Human Genome U133 Plus 2.0 platform (GPL570). Thus, the supplied files are appropriate for demonstrating raw microarray import, preprocessing, and quality control, but they do not provide both biological groups. Treating them as a cancer group would be an invalid analysis. This project addresses that limitation transparently by using the parent GSE19804 series matrix for the group comparison while keeping the supplied raw CEL files in the project.

**Aim:** To identify and interpret genes whose expression differs between lung cancer and adjacent normal lung tissue in GSE19804.

**Working hypothesis:** The tumour and adjacent-normal groups contain measurable transcriptomic differences, reflected by a set of genes with statistically significant differential expression. This is a testable hypothesis, not a pre-assumed result.

## 4. Data source and description

The primary public data source is NCBI Gene Expression Omnibus:

- Series: [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804)
- Platform: [GPL570 – Affymetrix Human Genome U133 Plus 2.0 Array](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570)
- Supplied raw samples: GSM494667–GSM494675
- Organism: *Homo sapiens*
- Supplied sample tissue: adjacent normal lung tissue
- Raw file format: compressed Affymetrix CEL (`.CEL.gz`)

The nine CEL files are stored under `data/raw_CEL/`. The analysis script retrieves the series matrix and phenotype metadata using Bioconductor's GEOquery package. This makes the project reproducible without distributing the complete parent cohort inside the ZIP. The GEO series matrix is a processed expression dataset; the script uses its provided expression values for the tumour-versus-normal limma comparison. The uploaded CEL files are separately read and normalised with the `affy` package using Robust Multi-array Average (RMA).

The data are public research data. The original GEO record and associated publication should be consulted for the complete experimental design, clinical annotations, and any reuse conditions. This project is a secondary computational analysis and does not attempt to identify individual patients. The supplied sample annotations include clinical characteristics; only variables required to identify the tissue groups should be used in this assessment.

## 5. Methods

### 5.1 Software and packages

All analysis is performed in R. The workflow uses `GEOquery` and `Biobase` to retrieve and handle GEO expression data, `affy` for CEL import and RMA preprocessing, `limma` for linear modelling and empirical Bayes moderation, `ggplot2` for visualisation, and `hgu133plus2.db` for probe annotation where available. `pheatmap`, `dplyr`, and `tibble` support exploratory summaries and result presentation.

### 5.2 Data import and validation

The pipeline first checks that the nine expected CEL files exist and are non-empty. It then downloads GSE19804 through GEOquery and selects the available series matrix. The expression matrix and sample metadata are aligned by sample identifiers. The script checks that expression values are numeric, sample identifiers are unique, and both tissue groups can be identified from GEO sample titles and characteristics. If the public series matrix does not expose both groups or the group labels cannot be identified safely, the script stops with an informative error rather than inventing labels.

### 5.3 Raw CEL preprocessing and quality control

The supplied CEL files are imported with `ReadAffy`. RMA preprocessing includes background correction, quantile normalisation, and summarisation of probe-level measurements to probe-set expression values. The resulting log2 expression matrix is used to produce sample distribution plots, a sample-correlation heatmap, and principal component analysis (PCA). These plots help reveal unusual samples, broad sample similarities, and possible technical variation. Nine adjacent-normal samples are too few to represent the entire parent cohort, so this QC subset is interpreted descriptively.

### 5.4 Exploratory analysis of the full cohort

For the parent GEO series matrix, the workflow checks sample counts by tissue group and generates a principal component plot and expression distribution plots. PCA is an unsupervised method that projects high-dimensional expression profiles into a small number of components. It is used here to explore whether samples show group-related structure, not as proof of a biological mechanism.

### 5.5 Differential-expression analysis

The tumour-versus-normal comparison is fitted using `limma`. A design matrix encodes the two tissue groups, and a linear model is fitted separately for each probe set. Empirical Bayes moderation stabilises variance estimates across probe sets. The output includes log2 fold change, average expression, moderated t-statistic, raw p-value, and Benjamini–Hochberg adjusted p-value. Adjusted p-values are used to account for multiple testing. The default significance threshold is adjusted p-value below 0.05; the threshold is explicitly recorded in the script.

The analysis is an unpaired group comparison unless the available GEO metadata and a defensible patient-pair identifier support a paired design. The script does not assume pairing merely because the samples are described as adjacent normal tissue. If patient pairing is not reliably represented in the downloaded metadata, the results must be interpreted as an unpaired comparison.

### 5.6 Annotation and visualisation

Probe-set identifiers are mapped to gene symbols using the platform annotation package where possible. A volcano plot displays estimated fold changes against statistical evidence. A heatmap displays a selected set of top differentially expressed probe sets across samples. A table of top results is written to `results/`. Multiple probe sets can map to the same gene, and not every probe set necessarily has a current gene-symbol annotation; this is reported rather than silently discarded.

## 6. Results

**Results are generated when the R workflow is executed.** This repository does not fabricate differential-expression counts, p-values, gene names, PCA separation, or figures. After running the analysis, inspect:

- `results/group_counts.csv` for the number of samples in each group.
- `results/differential_expression_all.csv` for the complete limma table.
- `results/top_differentially_expressed.csv` for the ranked results.
- `figures/` for QC, PCA, volcano, and heatmap plots.
- `report/ICA3_Bioinformatics_Report.html` for the rendered report.

The R Markdown report includes the actual sample counts, top-ranked probe sets, adjusted p-values, and plots generated from the GEO data retrieved at execution time. A result should only be described as statistically significant if its adjusted p-value meets the stated threshold. Lack of significant results would also be a valid finding and should be reported honestly.

## 7. Discussion and limitations

The project demonstrates a complete reproducible workflow for public microarray data and directly uses the nine raw CEL files supplied with the assessment. The most important data limitation is that all nine uploaded samples are adjacent normal tissue. They cannot independently answer a tumour-versus-normal question. The parent GSE19804 series is therefore retrieved to obtain the comparison group. The analysis report should state the exact number of samples retained after metadata filtering and any exclusions made by the script.

Microarray data have platform-specific limitations, including probe cross-hybridisation, limited transcript coverage, and dependence on preprocessing and annotation. The series matrix contains processed values, so its preprocessing history should be checked in the GEO record and original publication. The CEL-based RMA output is a separate analysis of the supplied nine-sample subset and should not be numerically mixed with the processed series matrix unless the data are processed together under a consistent workflow.

The cohort is described as non-smoking female lung cancer patients from Taiwan. This limits generalisability. Adjacent normal tissue is also not necessarily equivalent to lung tissue from a healthy person with no cancer; it may contain field effects or other changes associated with the tumour-bearing lung. Differential expression is an association and does not establish that a gene causes cancer. Findings would require independent validation, ideally using a separate cohort and complementary methods such as RNA sequencing or laboratory experiments.

## 8. Conclusion

This project provides an R-based and reproducible analysis of the GSE19804 lung tissue microarray study. It preserves the supplied nine normal CEL files, validates their biological identity, and retrieves the parent GEO expression matrix to support a scientifically meaningful tumour-versus-adjacent-normal comparison. The workflow covers import, preprocessing, exploratory analysis, statistical testing, visualisation, annotation, and interpretation. The final biological conclusions must be written from the generated report after the code has run successfully.

## 9. Reproducibility: how to run

1. Install R (recommended: a current R release) and RStudio.
2. Open `Lung_Cancer_Bioinformatics_Project.Rproj`.
3. Run `code/run_analysis.R` from the project root, or render `report/ICA3_Bioinformatics_Report.Rmd`.
4. The first run requires internet access to retrieve GSE19804 and Bioconductor packages.
5. Generated plots and tables are saved under `figures/` and `results/`.

The package installation script is `code/install_packages.R`. Package versions and the R session information are saved to `results/sessionInfo.txt` after a successful run. Exact output may change if GEO updates its records or annotation packages change; record the retrieval date and session information for the submitted version.

## 10. References

1. National Center for Biotechnology Information. Gene Expression Omnibus: GSE19804, *Genome-wide screening of transcriptional modulation in non-smoking female lung cancer in Taiwan*. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804
2. National Center for Biotechnology Information. GEO sample GSM494667. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM494667
3. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015;43(7):e47. https://doi.org/10.1093/nar/gkv007
4. Irizarry RA, Hobbs B, Collin F, et al. Exploration, normalization, and summaries of high density oligonucleotide array probe level data. *Biostatistics*. 2003;4(2):249–264. https://doi.org/10.1093/biostatistics/4.2.249
5. Bioconductor. GEOquery package documentation. https://bioconductor.org/packages/GEOquery/
6. Bioconductor. affy package documentation. https://bioconductor.org/packages/affy/
7. Bioconductor. limma package documentation. https://bioconductor.org/packages/limma/

## Assessment submission reminder

The assignment asks for a public GitHub repository named with your student registration number and at least five meaningful commits. Do not submit this ZIP unchanged as evidence of five commits: create the repository, add the files, and make your own meaningful commits. Render and inspect the report, check the generated results, and update the interpretation using your actual output before submitting the GitHub URL through LMS.
