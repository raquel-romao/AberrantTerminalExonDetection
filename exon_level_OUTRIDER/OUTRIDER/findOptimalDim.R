#'---
#' title: Calculate Optimal Encoding Dimension for OUTRIDER
#' author: Raquel Romão
#' wb:
#'  log:
#'   - snakemake: '/s/project/first_last_exon/Data/log/{tissue}/StrandDiff/optimal_dim_PCA.Rds'
#'  params:
#'   - maxTestedDimensionProportion: '`sm config["maxTestedDimensionProportion"]`'
#'  input:
#'   - ods: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_unfitted_PCA.Rds'
#'  output:
#'   - ods_optimal_dim: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_optimal_dim_PCA.Rds'
#'  type: script
#'  threads: 10
#'  resources:
#'   - mem_mb: 450000
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
mp <- snakemake@params$maxTestedDimensionProportion
implementation <- snakemake@config$implementation
register(MulticoreParam(snakemake@threads))

## subset filtered
ods <- ods[mcols(ods)$passedFilter,]

# add gene ranges to rowData
gr <- unlist(endoapply(rowRanges(ods), range))
if(length(gr) > 0){
    rd <- rowData(ods)
    rowRanges(ods) <- gr
    rowData(ods) <- rd
}

ods <- estimateSizeFactors(ods)

## find optimal encoding dimension
a <- 5 
b <- min(ncol(ods), nrow(ods)) / mp   # N/3

maxSteps <- 15
if(mp < 4){
    maxSteps <- 20
}

Nsteps <- min(maxSteps, b)   # Do at most 20 steps or N/3
# Do unique in case 2 were repeated
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
ods <- findEncodingDim(ods, params = pars_q, implementation = implementation)
#opt_q <- getBestQ(ods)

saveRDS(ods, snakemake@output$ods_optimal_dim)
