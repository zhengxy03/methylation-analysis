#The script is to divide bismark result .cov file into several subfiles

library(tidyr)
library(dplyr)

folder_path <- commandArgs(T)
file_name <- list.files(folder_path, pattern = "*.txt")

func_divide_by_chr <- function(chrname, file_prefix, data_df){
    print(paste("Processing chromosome", chrname, sep = " "))
    wfile_name <- paste(folder_path, file_prefix, "_chr", chrname, ".txt", sep ="")
    tempdata <- data_df %>% filter(chr == chrname) %>%
        write.table(., wfile_name, row.names = F)
}

func_read_file <- function(file_name){
    file_path <- paste(folder_path, file_name, sep = "/")
    file_prefix <- substring(file_name, 0, nchar(file_name) - 4)
    	print(paste("File", file_name, "is being importing and this may take a while..."), sep = "")
	rawdata_df <- read.table(file_path, header = F, stringsAsFactors = F)
	print("Importing file is finished!")
	colnames(rawdata_df) <- c("chr", "start", "end", "methyl%", "methyled", "unmethyled")
	chrvec <- unique(rawdata_df$chr)
	print(paste("Processing", file_prefix, sep = " "))
	lapply(chrvec, func_divide_by_chr, file_prefix = file_prefix, data_df = rawdata_df)
}