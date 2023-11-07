#'---
#' title: Hyper parameter optimization
#' author: Christian Mertes, raquelromao
#' wb:
#'  log:
#'    - snakemake: '/s/project/first_last_exon/Data/log/fraser_nofilter/{tissue}/03FL_hyper.Rds'
#'  params:
#'   - setup: '`sm cfg.AS.getWorkdir() + "/config.R"`'
#'   - workingDir: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/'
#'  threads: 12
#'  input:
#'   - filter: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}/filter_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
#'  output:
#'   - hyper: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/{tissue}/hyper_{version}.done", version=cfg.AS.get("FRASER_version"), allow_missing=True)`'
#'  type: script
#'  threads: 10
#'  resources:
#'   - mem_mb: 100000
#'---

saveRDS(snakemake, snakemake@log$snakemake)
source(snakemake@params$setup, echo=FALSE)

if ("random_seed" %in% names(snakemake@config)){
  rseed <- snakemake@config$random_seed
  if(isTRUE(rseed)){
    set.seed(42)
  } else if (is.numeric(rseed)){
    set.seed(as.integer(rseed))
  }
}

#+ input
dataset    <- snakemake@wildcards$tissue
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load PSI data
fds <- loadFraserDataSet(dir=workingDir, name=dataset)
fitMetrics(fds) <- psiTypes
    
# Run hyper parameter optimization
implementation <- snakemake@config$aberrantSplicing$implementation
mp <- snakemake@config$aberrantSplicing$maxTestedDimensionProportion

# Get range for latent space dimension
a <- 2 
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

maxSteps <- 12
if(mp < 6){
  maxSteps <- 15
}

Nsteps <- min(maxSteps, b)
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique

#for(type in psiTypes){
   #message(date(), ": ", type)
   # fds <- optimHyperParams(fds, type=type, 
   #                         implementation=implementation,
   #                         q_param=pars_q,
   #                         plot = FALSE)
   # fds <- saveFraserDataSet(fds)
#}

#fds <- saveFraserDataSet(fds)

fds <- calculatePSIValues(fds, types = c("psi5", "psi3")) #is it ok to do this?

#only one psi type relevant = jaccard
message(date(), ": jaccard")
fds <- optimHyperParams(fds, type='jaccard', 
                            implementation=implementation,
                            q_param=pars_q,
                            plot = FALSE)
fds <- saveFraserDataSet(fds)


# remove previous hyper.done files and create new one
outdir <- dirname(snakemake@output$hyper)
prevFilterFiles <- grep("hyper(.*)done", list.files(outdir), value=TRUE)
unlink(file.path(outdir, prevFilterFiles))
file.create(snakemake@output$hyper)
