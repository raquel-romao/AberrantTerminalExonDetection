library(FRASER)

tissue <- commandArgs(trailingOnly = TRUE)[1] #This variable is set using the command-line argument --tissue when the script is run.
fds <- loadFraserDataSet(file = sprintf("/s/project/fraser/fraser2/GTEx_v8/fds/minK20_25_minN10/PCA__pc0.1/savedObjects/%s__optQ__newFilt/fds-object.RDS", tissue))
gr <- rowRanges(fds)
gr_dt <- as.data.table(gr)

fwrite(gr_dt, file = sprintf("/s/project/first_last_exon/Data/fds_files/fds_%s.tsv", tissue), sep = "\t")
