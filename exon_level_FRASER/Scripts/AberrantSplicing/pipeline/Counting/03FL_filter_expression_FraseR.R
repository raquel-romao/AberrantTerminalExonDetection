#'---
#' title: Filter and clean dataset
#' author: Christian Mertes, raquelromao
#' wb:
#'  log:
#'    - snakemake: '/s/project/first_last_exon/Data/log/fraser_nofilter/{tissue}/filter.Rds'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/'
#'  input:
#'   - jaccard: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/raw-local-{tissue}/jaccard.h5'
#'  output:
#'   - fds: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}/fds-object.RDS'
#'   - done: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}/filter_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
#'  type: script
#'  threads: 3
#'  resources:
#'   - mem_mb: 100000
#'---



saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

opts_chunk$set(fig.width=12, fig.height=8)

# input
dataset    <- snakemake@wildcards$tissue
workingDir <- snakemake@params$workingDir
params     <- snakemake@config$aberrantSplicing
#exCountIDs <- snakemake@params$exCountIDs
#exCountFiles <- snakemake@input$exCounts
sample_anno_file <- snakemake@config$sampleAnnotation
minExpressionInOneSample <- params$minExpressionInOneSample
quantile <- params$quantileForFiltering
quantileMinExpression <- params$quantileMinExpression
minDeltaPsi <- params$minDeltaPsi
filterOnJaccard <- (params$FRASER_version == "FRASER2")

#fds <- loadFraserDataSet(dir=workingDir, name=paste0("raw-local-", dataset)) #or path to a file/jaccard,file=
fds <- loadFraserDataSet(file=snakemake@input$jaccard)

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Add external data if provided by dataset
#if(length(exCountIDs) > 0){
#   message("create new merged fraser object")
#    fds <- saveFraserDataSet(fds,dir = workingDir, name=paste0("raw-", dataset))

#   for(resource in unique(exCountFiles)){
 #       exSampleIDs <- exCountIDs[exCountFiles == resource]
 #       exAnno <- fread(sample_anno_file, key="RNA_ID")[J(exSampleIDs)]
  #      setnames(exAnno, "RNA_ID", "sampleID")
        
   #     ctsNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
    #    ctsFiles <- paste0(dirname(resource), "/", ctsNames, "_counts.tsv.gz")
        
        # Merging external counts restricts the junctions to those that 
        # are only present in both the counted (fromBam) junctions AND the 
        # junctions from the external counts.
     #   fds <- mergeExternalData(fds=fds, countFiles=ctsFiles,
      #          sampleIDs=exSampleIDs, annotation=exAnno)
      #  fds@colData$isExternal <- as.factor(!is.na(fds@colData$SPLICE_COUNTS_DIR))
    #}
#} else {
message("symLink fraser dir")
file.symlink(paste0(workingDir, "savedObjects/","raw-local-", dataset),
            paste0(workingDir, "savedObjects/","raw-", dataset))

fds@colData$isExternal <- as.factor(FALSE)
workingDir(fds) <- workingDir
name(fds) <- paste0("raw-", dataset)
#}

# filter for expression and write it out to disc.
fds <- filterExpressionAndVariability(fds, 
        minExpressionInOneSample = minExpressionInOneSample,
        quantile=quantile,
        quantileMinExpression=quantileMinExpression,
        minDeltaPsi = minDeltaPsi,
        filterOnJaccard=filterOnJaccard,
        filter=FALSE)

devNull <- saveFraserDataSet(fds,dir = workingDir)

# Keep junctions that pass filter
name(fds) <- dataset
if (params$filter == TRUE) {
    filtered <- mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}

seqlevels(fds) <- seqlevelsInUse(fds)
colData(fds)$sampleID <- as.character(colData(fds)$sampleID)
fds <- saveFraserDataSet(fds,dir = workingDir)

# remove previous filter.done files and create new one
outdir <- dirname(snakemake@output$done)
prevFilterFiles <- grep("filter(.*)done", list.files(outdir), value=TRUE)
unlink(file.path(outdir, prevFilterFiles))
file.create(snakemake@output$done)
