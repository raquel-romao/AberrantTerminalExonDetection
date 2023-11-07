#'---
#' title: Process Exon Files from HIT to prepare count matrix
#' author: raquelromao
#' wb:
#'  log:
#'   - snakemake: '/s/project/first_last_exon/Data/log/{tissue}/StrandDiff/process_exon_files.log'
#'  params:
#'   - INPUT_PATH: '/s/project/first_last_exon/Data/output_data/HITindex/rna_seq/{tissue}'
#'  input:
#'   - exon_files: '/s/project/first_last_exon/Data/output_data/HITindex/rna_seq/{tissue}/*_HITindex.exon'
#'  output:
#'   - table: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_table.csv'
#'   - matrix_up: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_matrix_up.csv'
#'   - matrix_down: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_matrix_down.csv'
#'  type: script
#'  threads: 1
#'  resources:
#'   - mem_mb: 16000
#'---

saveRDS(snakemake, snakemake@log$snakemake)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyr)
})

# Get list of exon files
exon_files <- list.files(path = snakemake@params$INPUT_PATH, pattern = "_HITindex.exon$", full.names = TRUE)

# Define a function to process exon files
process_exon_files <- function(exon_files) {
  # Initialize an empty data table to store results
  result <- data.table()

  # Iterate over each exon file
  for (file in exon_files) {
    # Extract sample name from file name
    sample_name <- sub(pattern = "_HITindex.exon$", replacement = "", x = basename(file))

    # Load file into a data table
    df <- fread(file, sep = '\t')

    # Filter the necessary columns
    df <- df[, .(exon, gene, strand, nUP, nDOWN, HITindex, ID, ID_position)]


    df[, sample := sample_name]  # Add sample name column
    
    # Create a unique identifier by combining exon and strand
    df[, exon_strand := paste0(exon, "_", strand)]

    # Append to the result data table
    result <- rbind(result, df)
  }

  # Write result to the output file
  fwrite(result, snakemake@output$table, sep = ",", quote = FALSE, row.names = FALSE)

  # Create the count matrices for nUP and nDOWN
  count_matrix_up <- dcast(result, exon_strand ~ sample, value.var = "nUP", fill = 0)
  count_matrix_down <- dcast(result, exon_strand ~ sample, value.var = "nDOWN", fill = 0)

  # Write the matrices to the output files
  fwrite(count_matrix_up, snakemake@output$matrix_up, sep = ",", quote = FALSE, row.names = FALSE)
  fwrite(count_matrix_down, snakemake@output$matrix_down, sep = ",", quote = FALSE, row.names = FALSE)

}

# Execute the function
process_exon_files(exon_files)


