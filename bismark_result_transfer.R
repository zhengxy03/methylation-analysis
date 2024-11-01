#This script is used to transfer the .cov file to facilitate the downstream analysis

library(tidyr)
library(dplyr)

file_name <- commandArgs(T)

func_read_file <- function(file_name){
    dir_vec <- strsplit(file_name, split = "/")[[1]]
    len <- length(dir_vec)
    file_prefix = substring(dir_vec[len], 0, nchar(dir_vec[len]) - 4)
    file_save_path = substring(file_name, 0, nchar(file_name) - nchar(dir_vec[len]))
    print(paste("File", file_name, "is being importing and this may take a while..."), sep = "")
    rawdata_df <- read.table(file_name, header = F, colClasses = "character")
    print("Importing file is finished!")
    colnames(rawdata_df) <- c("chr", "start", "end", "methyl%", "methyled", "unmethyled")
    write.table(rawdata_df, paste(file_save_path, file_prefix, "_transfered.txt", sep = ""), row.names = F)
}
lapply(file_names, func_read_file)