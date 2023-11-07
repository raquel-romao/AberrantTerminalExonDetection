#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller, Raquel Romão
#' wb:
#'  log:
#'   - snakemake: '/s/project/first_last_exon/Data/log/{tissue}/StrandDiff/filter_PCA.Rds'
#'  params:
#'   - fpkmCutoff: '`sm config["fpkmCutoff"]`'
#'  input:
#'   - counts: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_matrix_sum.csv'
#'  output:
#'   - ods: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_unfitted_PCA.Rds'
#'  type: script
#'  threads: 1
#'  resources:
#'   - mem_mb: 16000
#'---



saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

counts <- fread(snakemake@input$counts)

##############################################
count_mat <- as.matrix(counts[, -"exon_strand"])
rownames(count_mat) <- counts[, exon_strand]
##############################################

ods <- OutriderDataSet(countData=count_mat)

# filter not expressed genes
fpkmCutoff <- snakemake@params$fpkmCutoff
mcols(ods)['basepairs'] <- 1 #line added to re-define filterExpression
ods <- filterExpression(ods, filter=FALSE,
                        fpkmCutoff=fpkmCutoff, addExpressedGenes=TRUE)

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# External data check
#if (is.null(ods@colData$GENE_COUNTS_FILE)){ #column does not exist in sample annotation table
#    has_external <- FALSE
#}else if(all(is.na(ods@colData$GENE_COUNTS_FILE))){ #column exists but it has no values
#    has_external <- FALSE
#}else if(all(ods@colData$GENE_COUNTS_FILE == "")){ #column exists with non-NA values but this group has all empty strings
#   has_external <- FALSE
#}else{ #column exists with non-NA values and this group has at least 1 non-empty string
#    has_external <- TRUE
#}

#if(has_external){
#    ods@colData$isExternal <- as.factor(ods@colData$GENE_COUNTS_FILE != "")
#}else{
#    ods@colData$isExternal <- as.factor(FALSE)
#}


# Save the ods before filtering to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)
