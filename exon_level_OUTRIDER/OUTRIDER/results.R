#'---
#' title: OUTRIDER Results
#' author: mumichae, Raquel Romão
#' wb:
#'  log:
#'   - snakemake: '/s/project/first_last_exon/Data/log/{tissue}/StrandDiff/OUTRIDER_results_PCA.Rds'
#'  params:
#'   - padjCutoff: '`sm config["padjCutoff"]`' 
#'   - zScoreCutoff: '`sm config["zScoreCutoff"]`'
#'  input:
#'   - ods: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/ods_PCA.Rds'
#'  output:
#'   - results: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/OUTRIDER_results_PCA.tsv'
#'   - results_all: '/s/project/first_last_exon/Data/output_data/outrider/{tissue}/StrandDiff/OUTRIDER_results_all_PCA.Rds'
#'  type: script
#'  threads: 1
#'  resources:
#'   - mem_mb: 100000
#'---

saveRDS(snakemake, snakemake@log$snakemake)
#source(snakemake@input$add_HPO_cols)

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

ods <- readRDS(snakemake@input$ods)
res <- results(ods, padjCutoff = snakemake@params$padjCutoff,
			   zScoreCutoff = snakemake@params$zScoreCutoff, all = TRUE)

# Add fold change
res[, foldChange := round(2^l2fc, 2)]

# Save all the results and significant ones
saveRDS(res, snakemake@output$results_all)

# Subset to significant results
padj_cols <- grep("padjust", colnames(res), value=TRUE)
res <- res[do.call(pmin, c(res[,padj_cols, with=FALSE], list(na.rm = TRUE))) 
                <= snakemake@params$padjCutoff &
            abs(zScore) >= snakemake@params$zScoreCutoff]


#I removed the HPO terms and gene mapping sections (gene mapping + gene input)

#gene_annot_dt <- fread(snakemake@input$gene_name_mapping)
#if(!is.null(gene_annot_dt$gene_name)){
#  if(grepl('ENSG00', res[1,geneID]) & grepl('ENSG00', gene_annot_dt[1,gene_id])){
#    res <- merge(res, gene_annot_dt[, .(gene_id, gene_name)],
#                 by.x = 'geneID', by.y = 'gene_id', sort = FALSE, all.x = TRUE)
#    setnames(res, 'gene_name', 'hgncSymbol')
#    res <- cbind(res[, .(hgncSymbol)], res[, - 'hgncSymbol'])
#  }
#}

# Add HPO terms, requires online connection and for there to be annotated HPO terms
#sa <- fread(snakemake@config$sampleAnnotation, 
#              colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
#if(!is.null(sa$HPO_TERMS) & nrow(res) > 0){
#  if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
#    res <- add_HPO_cols(res, hpo_file = snakemake@params$hpoFile)
#  }
#}


# Save results
fwrite(res, snakemake@output$results, sep = "\t", quote = F)
