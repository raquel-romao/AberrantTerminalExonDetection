#'---
#' title: Annotate introns with gene symbols
#' author: Ines Scheller
#' wb:
#'  log:
#'    - snakemake: '/s/project/first_last_exon/Data/log/fraser_nofilter/{tissue}--{annotation}/05FL_geneAnnotation.Rds'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/'
#'   - outputDir: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/'
#'  threads: 20
#'  input:
#'   - fdsin: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}/predictedMeans_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'   - txdb: '/s/project/first_last_exon/Data/input_data/preprocess/{annotation}/txdb.db'
#'   - gene_name_mapping: '/s/project/first_last_exon/Data/input_data/preprocess/{annotation}/gene_name_mapping_{annotation}.tsv'
#'  output:
#'   - fdsout: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}--{annotation}/predictedMeans_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'   - fds_rds: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}--{annotation}/fds-object.RDS'
#'  type: script
#'  threads: 15
#'  resources:
#'   - mem_mb: 100000
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)
library(AnnotationDbi)

annotation <- snakemake@wildcards$annotation
dataset    <- snakemake@wildcards$tissue
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir
outputDir  <- snakemake@params$outputDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds_input <- loadFraserDataSet(dir=workingDir, name=dataset)

# Read annotation and match the chr style
txdb <- loadDb(snakemake@input$txdb)
orgdb <- fread(snakemake@input$gene_name_mapping)

seqlevels_fds <- seqlevelsStyle(fds_input)[1]
seqlevelsStyle(orgdb$seqnames) <- seqlevels_fds
seqlevelsStyle(txdb) <- seqlevels_fds

# Annotate the fds with gene names and save it as a new object
fds_input <- annotateRangesWithTxDb(fds_input, txdb = txdb, orgDb = orgdb, 
                    feature = 'gene_name', featureName = 'hgnc_symbol', 
                    keytype = 'gene_id')

# add basic annotations for overlap with the reference annotation
# run this function before creating the results table to include it there
fds_input <- annotateIntronReferenceOverlap(fds_input, txdb)

# save fds
fds <- saveFraserDataSet(fds_input, dir=outputDir, name = paste(dataset, annotation, sep = '--'), rewrite = TRUE)

# remove .h5 files from previous runs with other FRASER version
fdsDir <- dirname(snakemake@output$fdsout[1])
for(type in psiTypesNotUsed){
    predMeansFile <- file.path(fdsDir, paste0("predictedMeans_", type, ".h5"))
    if(file.exists(predMeansFile)){
        unlink(predMeansFile)
    }
}
