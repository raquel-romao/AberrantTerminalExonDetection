#'---
#' title: Fit OUTRIDER Model
#' author: Raquel Romão
#' wb:
#'  log:
#'   - snakemake: '/s/project/first_last_exon/Data/log/{tissue}/StrandDiff/runOUTRIDER_PCA.Rds'
#'  input:
#'   - ods: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_optimal_dim_PCA.Rds'
#'  output:
#'   - ods_fitted: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_fitted_PCA.Rds'
#'  type: script
#'  threads: 20
#'  resources:
#'   - mem_mb: 100000
#'---

#+ echo=F
saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(magrittr)
    library(tools)
})

ods <- readRDS(snakemake@input$ods)
implementation <- snakemake@config$implementation
register(MulticoreParam(snakemake@threads))


## fit OUTRIDER
opt_q <- getBestQ(ods)
message(date(), ": SizeFactor estimation ...")
ods <- estimateSizeFactors(ods)
message(date(), ": Controlling for confounders ...")
implementation <- tolower(implementation)
ods <- controlForConfounders(ods, q=opt_q, implementation=implementation)
if(grepl("^(peer|pca)$", implementation)){
    message(date(), ": Fitting the data ...")
    ods <- fit(ods)
}
message("outrider fitting finished")

saveRDS(ods, snakemake@output$ods_fitted)

