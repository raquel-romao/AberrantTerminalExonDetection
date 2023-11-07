#'---
#' title: Create FraserDataSet from existing count matrices
#' author: raquelromao
#' wb:
#'  log:
#'    - snakemake: '/s/project/first_last_exon/Data/log/fraser_nofilter/{tissue}/01FL_definefds.Rds'
#'  params:
#'   - workingDir: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/'
#'  input:
#'   - up_counts: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_matrix_up.csv'
#'   - down_counts: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_matrix_down.csv'
#'  output:
#'    - fds: '/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/raw-local-{tissue}/fds-object.RDS'
#'    - splice_metrics: '`sm expand("/s/project/first_last_exon/Data/output_data/fraser_nofilter/savedObjects/raw-local-{tissue}/{type}.h5", type=cfg.AS.getPsiTypeAssay(), allow_missing=True)`'
#'  type: script
#'  threads: 1
#'  resources:
#'   - mem_mb: 100000
#'---
 


saveRDS(snakemake, snakemake@log$snakemake)

workingDir <- snakemake@params$workingDir
dataset    <- snakemake@wildcards$tissue

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(tidyr)
    library(FRASER)
    library(GenomicRanges)
})

up_counts <- fread(snakemake@input$up_counts)
down_counts <- fread(snakemake@input$down_counts)

U <- as.matrix(up_counts[, -"exon_strand"])

D <- as.matrix(down_counts[, -"exon_strand"])


# Set up junctionCts table

df_strand <- separate(up_counts, exon_strand, into = c("exon", "strand"), sep = "_")

seqdf <- df_strand %>%
  separate(exon, into = c("seqnames", "positions"), sep = ":") %>%
  separate(positions, into = c("start", "end"), sep = "-") %>%
  mutate(across(c(start, end), as.numeric)) %>%
  mutate(width = end - start + 1)

gr <- makeGRangesFromDataFrame(seqdf, keep.extra.columns = TRUE)

junctionCts <- FRASER:::annotateSpliceSite(gr)


## Set up sampleTable

# Get the metadata columns from junctionCts
metadata <- mcols(junctionCts)

# Get the column names from the metadata
sampleIDs <- colnames(metadata)

# Define the names you want to remove
remove_names <- c("endID", "startID")

# Use setdiff() to remove the names
sampleIDs <- setdiff(sampleIDs, remove_names)

# Create a new sampleTable with the sampleIDs
sampleTable <- data.frame(sampleID = sampleIDs) 

# Make sure the sampleTable is a data.table
sampleTable <- as.data.table(sampleTable)

# To make GTEx sample names not prone to changes
sampleTable$sampleID <- gsub("-", "\\.", sampleTable$sampleID)


# Set up spliceSiteCts

junctionCts_dt <- data.table(as.data.frame(junctionCts)) # here the GTEx sample names get changed from "-" to "."

spliceSiteCts <- copy(junctionCts_dt)

spliceSiteCts[, c("startID", "endID") := NULL]

spliceSiteIDs <- unique(c(junctionCts_dt$startID, junctionCts_dt$endID))

desired_length <- nrow(junctionCts_dt)

# Calculate how many times you need to repeat the unique values
repeat_count <- ceiling(desired_length / length(spliceSiteIDs))

# Repeat the unique values as needed and truncate to the desired length
adjusted_spliceSiteIDs <- rep(spliceSiteIDs, repeat_count)[1:desired_length]

# Add the adjusted list as a new column in dt2
spliceSiteCts[, spliceSiteID := adjusted_spliceSiteIDs]

# Add 'type' column with all values set to 'Donor'
spliceSiteCts[, type := "Donor"]

colorder <- c("seqnames", "start", "end", "width", "strand", "spliceSiteID", "type", colnames(spliceSiteCts)[6:208])
spliceSiteCts[, colorder, with=FALSE]


# Create FraserDataSet - assign to U
fds <- FraserDataSet(colData=sampleTable, junctions=junctionCts_dt,
                     spliceSites=spliceSiteCts, workingDir=workingDir, name = paste0("raw-local-", dataset))

# Make GTEx names concordant with all the other matrixes (when converting to df there's an issue of convertion of the names)
new_names <- gsub("-", "\\.", colnames(U))

colnames(U) <- new_names

colnames(D) <- new_names


# Assign n-k counts
counts(fds, type="jaccard", side="other", withDimnames=FALSE) <- D

# Assign k/n ratio to jaccard assay
assay(fds, type="j", "jaccard", withDimnames=FALSE) <- U / (U+D)

# Calculate delta PSI value
fds <- FRASER:::calculateDeltaPsiValue(fds, psiType="jaccard", assayName="delta_jaccard")

newColnames <- gsub("\\.", "-", sampleTable$sampleID)

colnames(fds) <- newColnames

colData(fds)$sampleID <- newColnames

saveFraserDataSet(fds)
