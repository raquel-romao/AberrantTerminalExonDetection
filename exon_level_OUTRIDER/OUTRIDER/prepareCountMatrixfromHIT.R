#'---
#' title: Process Exon Files from HIT to prepare count matrix
#' author: Raquel Romão
#' wb:
#'  log:
#'   - snakemake: '/s/project/first_last_exon/Data/log/{tissue}/StrandDiff/process_exon_files.log'
#'  params:
#'   - INPUT_PATH: '/s/project/first_last_exon/Data/output_data/HITindex/rna_seq/{tissue}'
#'  input:
#'   - exon_files: '/s/project/first_last_exon/Data/output_data/HITindex/rna_seq/{tissue}/*_HITindex.exon'
#'  output:
#'   - table: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_table.csv'
#'   - matrix_sum: '/s/project/first_last_exon/Data/input_data/gtex/{tissue}/StrandDiff/{tissue}_matrix_sum.csv'
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

    # Calculate UplusD
    df[, UplusD := nUP + nDOWN]

    df[, sample := sample_name]  # Add sample name column
    
    # Create a unique identifier by combining exon and strand
    df[, exon_strand := paste0(exon, "_", strand)]

    # Append to the result data table
    result <- rbind(result, df)
  }

  # Write result to the output file
  fwrite(result, snakemake@output$table, sep = ",", quote = FALSE, row.names = FALSE)

  # Create the count matrices for UplusD
  count_matrix_sum <- dcast(result, exon_strand ~ sample, value.var = "UplusD", fill = 0)

  # Write the matrix to the output files
  fwrite(count_matrix_sum, snakemake@output$matrix_sum, sep = ",", quote = FALSE, row.names = FALSE)
}

# Execute the function
process_exon_files(exon_files)


