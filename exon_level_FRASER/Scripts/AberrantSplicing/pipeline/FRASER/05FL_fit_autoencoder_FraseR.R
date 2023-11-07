#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  log:
#'    - snakemake: '/s/project/first_last_exon/Data/log/fraser_nofilter/{tissue}/04FL_fit.Rds'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/'
#'  threads: 20
#'  input:
#'   - hyper: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}/hyper_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
#'  output:
#'   - fdsout: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}/predictedMeans_{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'  type: script
#'  threads: 15
#'  resources:
#'   - mem_mb: 100000
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

dataset    <- snakemake@wildcards$tissue
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

fds <- loadFraserDataSet(dir=workingDir, name=dataset)

# Fit autoencoder
# run it for every type
implementation <- snakemake@config$aberrantSplicing$implementation

for(type in psiTypes){
    currentType(fds) <- type
    q <- bestQ(fds, type)
    verbose(fds) <- 3   # Add verbosity to the FRASER object
    fds <- fit(fds, q=q, type=type, iterations=15, implementation=implementation)
    fds <- saveFraserDataSet(fds)
}

# remove .h5 files from previous runs with other FRASER version
fdsDir <- dirname(snakemake@output$fdsout[1])
for(type in psiTypesNotUsed){
    predMeansFile <- file.path(fdsDir, paste0("predictedMeans_", type, ".h5"))
    if(file.exists(predMeansFile)){
        unlink(predMeansFile)
    }
}
