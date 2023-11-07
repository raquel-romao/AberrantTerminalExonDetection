#'---
#' title: Results of FRASER analysis
#' author: Christian Mertes, raquelromao
#' wb:
#'  log:
#'    - snakemake: '/s/project/first_last_exon/Data/log/fraser_nofilter/{tissue}--{annotation}/07FL_results.Rds'
#'  params:
#'   - workingDir: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/'
#'   - padjCutoff: '`sm cfg.AS.get("padjCutoff")`'
#'   - deltaPsiCutoff: '`sm cfg.AS.get("deltaPsiCutoff")`'
#'  threads: 10
#'  input:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - fdsin: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}--{annotation}/" +
#'                  "padjBetaBinomial_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'   - txdb: '/s/project/first_last_exon/Data/input_data/preprocess/{annotation}/txdb.db'
#'  output:
#'   - resultTableJunc: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/results/{annotation}/{tissue}/results_per_junction2.tsv'
#'   - resultTableGene_full: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/results/{annotation}/{tissue}/results_gene_all2.tsv'
#'   - resultTableGene_aberrant: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/results/{annotation}/{tissue}/results2.tsv'
#'   - resultTableJuncAll: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/results/{annotation}/{tissue}/results_all.tsv'
#'  type: script
#'  threads: 15
#'  resources:
#'   - mem_mb: 100000
#'---


#snakemake <- readRDS('/s/project/first_last_exon/Data/log/fraser/Liver--v40/07FL_results.Rds')

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@input$setup, echo=FALSE)
#source(snakemake@input$add_HPO_cols)
library(AnnotationDbi)

annotation    <- snakemake@wildcards$annotation
dataset    <- snakemake@wildcards$tissue
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# fix sample names of nonSplitReadCounts
nsrObj <- nonSplicedReads(fds)
colnames(nsrObj) <- colnames(fds)
nonSplicedReads(fds) <- nsrObj

# Extract results per junction
res_junc <- results(fds, psiType=psiTypes,
                    padjCutoff=snakemake@params$padjCutoff,
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff)
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')

# Add features
if(nrow(res_junc_dt) > 0){

    # dcast to have one row per outlier
    subsets <- res_junc_dt[, unique(FDR_set)]
    subsets <- subsets[subsets != "transcriptome-wide"]
    colorder <- colnames(res_junc_dt[, !"FDR_set", with=FALSE])
    res_junc_dt <- dcast(res_junc_dt, ... ~ FDR_set, value.var="padjust")
    setnames(res_junc_dt, "transcriptome-wide", "padjust")
    for(subset_name in subsets){
        setnames(res_junc_dt, subset_name, paste0("padjust_", subset_name))
    }
    setcolorder(res_junc_dt, colorder)
    
    # number of samples per gene and variant
    res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
    res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
    
    # add colData to the results
    res_junc_dt <- merge(res_junc_dt, as.data.table(colData(fds)), by = "sampleID")
    res_junc_dt[, c("bamFile", "pairedEnd", "STRAND", "RNA_BAM_FILE", "DNA_VCF_FILE", "COUNT_MODE", "COUNT_OVERLAPS") := NULL]
} else{
    warning("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset.")
}


#################################################################################################################


# Extract full results (non-significant included)
res_full <- results(fds, psiType='jaccard',
                    all=TRUE)

res_full_dt <- as.data.table(res_full)
print('Non-sig results (all) extracted')
write_tsv(res_full_dt, file=snakemake@output$resultTableJuncAll)




# Extract full results by gene
res_gene <- results(fds, psiType='jaccard',
                    aggregate=TRUE, collapse=FALSE,
                    all=TRUE)

res_genes_dt   <- as.data.table(res_gene)
subsets <- res_genes_dt[, unique(FDR_set)]
subsets <- subsets[subsets != "transcriptome-wide"]
colorder <- colnames(res_genes_dt[, !c("padjust", "FDR_set"), with=FALSE])
res_genes_dt <- dcast(res_genes_dt[,!"padjust", with=FALSE], ... ~ FDR_set, value.var="padjustGene")
setnames(res_genes_dt, "transcriptome-wide", "padjustGene")
for(subset_name in subsets){
    setnames(res_genes_dt, subset_name, paste0("padjustGene_", subset_name))
}
setcolorder(res_genes_dt, colorder)
print('Results per gene extracted')
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene_full)

 # Subset gene results to aberrant
padj_cols <- grep("padjustGene", colnames(res_genes_dt), value=TRUE)
res_genes_dt <- res_genes_dt[do.call(pmin, c(res_genes_dt[,padj_cols, with=FALSE], 
                                                 list(na.rm = TRUE))) <= snakemake@params$padjCutoff &
                                     abs(deltaPsi) >= snakemake@params$deltaPsiCutoff & 
                                     totalCounts >= 5,]
 
if(length(res_gene) > 0){
    res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
    res_genes_dt[, c("bamFile", "pairedEnd", "STRAND", "RNA_BAM_FILE", "DNA_VCF_FILE", "COUNT_MODE", "COUNT_OVERLAPS") := NULL]

     # add HPO overlap information
     #sa <- fread(snakemake@config$sampleAnnotation, 
     #            colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
     #if(!is.null(sa$HPO_TERMS)){
     #    if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
     #        res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = snakemake@params$hpoFile)
     #    }
     #}
} else{
    res_genes_dt <- data.table()
    warning("The aberrant splicing pipeline gave 0 gene-level results for the ", dataset, " dataset.")
}

# Annotate results with spliceEventType and blacklist region overlap
# load reference annotation
library(AnnotationDbi)
txdb <- loadDb(snakemake@input$txdb)

# annotate the type of splice event and UTR overlap
res_junc_dt <- annotatePotentialImpact(result=res_junc_dt, txdb=txdb, fds=fds)
res_genes_dt <- annotatePotentialImpact(result=res_genes_dt, txdb=txdb, fds=fds)

# set genome assembly version to load correct blacklist region BED file (hg19 or hg38)
assemblyVersion <- snakemake@config$genomeAssembly
if(grepl("grch37", assemblyVersion, ignore.case=TRUE)){
    assemblyVersion <- "hg19"
}
if(grepl("grch38", assemblyVersion, ignore.case=TRUE)){
    assemblyVersion <- "hg38"
}

# annotate overlap with blacklist regions
if(assemblyVersion %in% c("hg19", "hg38")){
    res_junc_dt <- flagBlacklistRegions(result=res_junc_dt, 
                                        assemblyVersion=assemblyVersion)
    res_genes_dt <- flagBlacklistRegions(result=res_genes_dt, 
                                         assemblyVersion=assemblyVersion)
} else{
    message(date(), ": cannot annotate blacklist regions as no blacklist region\n", 
            "BED file is available for genome assembly version ", assemblyVersion, 
            " as part of FRASER.")
}


# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene_aberrant)
