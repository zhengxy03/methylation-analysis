#The script is used to make the differential methylation analysis
library(tidyr)
library(dplyr)
library(DSS)

args <- commandArgs(T)
first_file <- args[1]
second_file <- args[2]
file_prefix <- args[3]
file_save_path <- args[4]

print("Importing data file, this may take a while...")
first_raw_data <- read.table(first_file, header = T)
second_raw_data <- read.table(second_file, header = T)
print("Importing data is finished!")

print("Transforming raw data, this may take a while...")
DSS_first_input_file <- first_raw_data %>%
    mutate(chr = paste("chr", chr, sep = "")) %>%
    mutate(pos = start, N = as.numeric(methyled) + as.numeric(unmethyled), X = as.numeric(methyled)) %>%
    select(chr, pos, N, X)
DSS_second_input_file <- second_raw_data %>%
    mutate(chr = paste("chr", chr, sep = "")) %>%
    mutate(pos = start, N = as.numeric(methyled) + as.numeric(unmethyled), X = as.numeric(methyled)) %>%
    select(chr, pos, N, X)
print("Transforming raw data is finished!")

print("Making DSS analysis, this may take a long time...")
bsobj <- makeBSseqData(list(DSS_first_input_data, DSS_second_input_data), c("S1", "S2"))
dmlTest <- DMLtest(bsobj, group1 = c("S1"), group2 = c("S2"), smoothing = T)

print("Finding DMLs, this may take a while...")
dmls <- callDML(dmlTest, p.threshold = 0.001)

print("Finding DMRs, this may take a while...")
dmrs <- callDMR(dmlTest, p.threshold=0.01)

write.table(dmlTest, paste(file_save_path, file_prefix, "_DSS_test_result.txt", sep = ""), row.names = F)
write.table(dmls, paste(file_save_path, file_prefix, "_DSS_dmls_result.txt", sep = ""), row.names = F)
write.table(dmrs, paste(file_save_path, file_prefix, "_DSS_dmrs_result.txt", sep = ""), row.names = F)