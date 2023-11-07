#'---
#' title: P value calculation for OUTRIDER
#' author: Ines Scheller, Raquel Romão
#' wb:
#'  log:
#'   - snakemake: '/s/project/first_last_exon/Data/log/{tissue}/StrandDiff/pvalsOUTRIDER_PCA.Rds'
#'  input:
#'   - ods_fitted: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_fitted_PCA.Rds'
#'  output:
#'   - ods_with_pvals: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_PCA.Rds'
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
    library(yaml)
})

ods <- readRDS(snakemake@input$ods_fitted)
implementation <- snakemake@config$aberrantExpression$implementation
register(MulticoreParam(snakemake@threads))

# read in gene subsets from sample anno if present (returns NULL if not present)
#source(snakemake@params$parse_subsets_for_FDR)
#outrider_sample_ids <- snakemake@params$ids
#subsets <- parse_subsets_for_FDR(snakemake@params$genes_to_test,
                                 #sampleIDs=outrider_sample_ids)
#subsets_final <- convert_to_geneIDs(subsets, snakemake@input$gene_name_mapping)

# P value calculation
message(date(), ": P-value calculation ...")
ods <- computePvalues(ods)
message(date(), ": Zscore calculation ...")
ods <- computeZscores(ods, 
                      peerResiduals=grepl('^peer$', implementation))

# save ods with pvalues
saveRDS(ods, snakemake@output$ods_with_pvals)
