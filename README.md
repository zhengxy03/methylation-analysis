# methylation-analysis
* [download data](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#1-download-data)
    * [reference data](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#11-reference-data)
    * [experiment data](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#12-experiment-data)
* [quality control and trimming](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#2-quality-control-and-trimming)
* [methylation analysis](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#3-methylation-analysis)
    * [bismark download](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#31-bismark-download)
    * [genome indexing](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#32-genome-indexing)
    * [read alignment](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#33-read-alignment)
    * [aligned reads deduplication](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#34-aligned-reads-deduplication)
    * [methylation information extracting](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#35-methylation-information-extracting)
* [downstream analysis](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#4-downstream-analysis)
    * [data preparation](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#41-data-preparation)
    * [DML/DMR detection](https://github.com/zhengxy03/methylation-analysis/blob/main/README.md#42-dmldmr-detection)

## 1 download data
### 1.1 reference data
Ensembl-mouse
```
mkdir -p ~/project/musculus
cd ~/project/musculus
mkdir genome annotation sequence output

cd genome
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.1.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.2.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.3.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.4.fa.gz

aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.5.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.6.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.7.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.8.fa.gz

aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.9.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.10.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.11.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.12.fa.gz

aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.13.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.14.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.15.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.16.fa.gz

aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.17.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.18.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.19.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.X.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.Y.fa.gz
```
### 1.2 experiment data
ncbi-run selector(GSE116016)
```
cd ../sequence
mkdir WT_mESC_rep1 TetTKO_mESC_rep1

nohup prefetch SRR736884{1..2} SRR7368845 -O . &

mkdir srr
array=(SRR736884{1..2} SRR7368845)
for i in "${array[@]}";
do
    dir="$HOME/project/musculus/sequence/$i"
    cd "${dir}"
    mv ${dir}/* $HOME/project/musculus/sequence/srr
done

cd $HOME/project/musculus/sequence/srr
parallel -j 4 "
    fastq-dump --split-3 --gzip {1}
" ::: $(ls *.sra)
```
[EBI-ENA search page](https://www.ebi.ac.uk/ena)-project-Generated FASTQ files: FTP
```
aria2c -d ./WT_mESC_rep1/ -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/001/SRR7368841/SRR7368841.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/002/SRR7368842/SRR7368842.fastq.gz

aria2c -d ./TetTKO_mESC_rep1/ -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR736/005/SRR7368845/SRR7368845.fastq.gz
```
## 2 quality control and trimming
```
fastqc --threads 3 ./WT_mESC_rep1/*.fastq.gz ./TetTKO_mESC_rep1/*.fastq.gz

trim_galore -o ./WT_mESC_rep1/trimmed_data/ --fastqc ./WT_mESC_rep1/*.fastq.gz
trim_galore -o ./TetTKO_mESC_rep1/trimmed_data/ --fastqc ./TetTKO_mESC_rep1/*.fastq.gz
```
## 3 methylation analysis
## 3.1 bismark download
[bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
> Bismark 是一款用于分析 Bisulfite Sequencing（亚硫酸氢盐测序）数据的软件工具。它可以高效地进行读段比对和甲基化探测，能区分 CG（CpG）、CHG 和 CHH 等不同类型的甲基化位点，并且允许用户通过可视化来解释数据。<br>
> 序列比对策略（bowtie2）对于每个读段中的碱基，在与参考基因组比对时，会根据可能的甲基化情况进行匹配。例如，读段中的一个胸腺嘧啶（T）可能对应参考基因组中的胞嘧啶（C），这意味着该位点可能是未甲基化的
## 3.2 genome indexing
input: genome reference data (.fa.gz)<br>
output:Bisulfite_Genome(CT/GA conversion)(*.bt2)
```
bismark_genome_preparation --bowtie2 ~/project/musculus/genome
```
## 3.3 read alignment
input:Bisulfite_Genome & experiment data(.fastq.gz)<br>
output:bismark_result/*.bam
```
genome_path="$HOME/project/musculus/genome/chr1"
cd ~/project/musculus/sequence

bismark -o ./WT_mESC_rep1/bismark_result/ --parallel 4 --genome_folder ${genome_path} ./WT_mESC_rep1/*.fastq.gz
bismark -o ./TetTKO_mESC_rep1/bismark_result/ --parallel 4 --genome_folder ${genome_path} ./TetTKO_mESC_rep1/*.fastq.gz

samtools cat -o SRX4241790 trimmed_bismark_bt2.bam ./WT_mESC_rep1/bismark_result/*.bam
```
## 3.4 aligned reads deduplication
去除PCR扩增过度的影响<br>
input:bismark_result<br>
output:deduplicated_result(.bam)
```
mkdir -p ./WT_mESC_rep1/deduplicated_result/
mkdir -p ./TetTKO_mESC_rep1/deduplicated_result/

deduplicate_bismark --bam --output_dir ./WT_mESC_rep1/deduplicated_result/ ./SRX4241790_trimmed_bismark_bt2.bam

deduplicate_bismark --bam --output_dir ./TetTKO_mESC_rep1/deduplicated_result/ ./TetTKO_mESC_rep1/bismark_result/*.bam
```
## 3.5 methylation information extracting
```
genome_path="$HOME/project/musculus/genome/chr1"
cd $HOEM/project/musculus/sequence

bismark_methylation_extractor  --single-end --gzip --parallel 4 --bedGraph \
--cytosine_report --genome_folder ${genome_path} \
-o ./WT_mESC_rep1/deduplicated_result/ ./WT_mESC_rep1/deduplicated_result/*.bam
bismark_methylation_extractor --single-end --gzip --parallel 4 --bedGraph \
--cytosine_report --genome_folder ${genome_path} \
-o ./TetTKO_mESC_rep1/deduplicated_result/ ./TetTKO_mESC_rep1/deduplicated_result/*.bam
```
# 4 downstream analysis
## 4.1 data preparation
```
mkdir -p ~/project/musculus/R_analysis/WT_data
mkdir -p ~/project/musculus/R_analysis/TetTKO_data
cd ~/project/musculus/R_analysis

WT_path="$HOME/project/musculus/sequence/WT_mESC_rep1/deduplicated_result/SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov.gz"
TetTKO_path="$HOME/project/musculus/sequence/TetTKO_mESC_rep1/deduplicated_result/SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz"

cp $WT_path ./WT_data/
cp $TetTKO_path ./TetTKO_data/

gunzip -d ./WT_data/SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov.gz
gunzip -d ./TetTKO_data/SRR7368845_bismark_bt2.deduplicated.bismark.cov.gz

#convert to a .txt
cp ./WT_data/SRX4241790_trimmed_bismark_bt2.deduplicated.bismark.cov ./WT_data/SRX4241790_methylation_result.txt
cp ./TetTKO_data/SRR7368845_bismark_bt2.deduplicated.bismark.cov ./TetTKO_data/SRR7368845_methylation_result.txt

R
```
```
#R
library(tidyr)
library(dplyr)

file_names <- ("./WT_data/SRX4241790_methylation_result.txt", "./TetTKO_data/SRR7368845_methylation_result.txt")

func_read_file <- function(file_name){
    dir_vec <- strsplit(file_name, split="/")[[1]]
    len <- length(dir_vec)
    file_prefix=substring(dir_vec[len], 0, nchar(dir_vec[len]) - 4)
    file_save_path=substring(file_name, 0, nchar(file_name) - nchar(dir_vec[len]))
    print(paste("File", file_name, "is being importing and this may take a while..."), sep = "")
    rawdata_df <- read.table(file_name, header = F, colClasses = "character")
    print("Importing file is finished!")
    colnames(rawdata_df) <- c("chr", "start", "end", "methyl%", "methyled", "unmethyled")
    write.table(rawdata_df, paste(file_save_path, file_prefix, "_transfered.txt", sep = ""), row.names = F)
}
lapply(file_names, func_read_file)

q()

#Rscript

#Rscript $HOME/project/Script/bismark_result_transfer.R ./WT_data/SRX4241790_methylation_result.txt ./TetTKO_data/SRR7368845_methylation_result.txt
```
## 4.2 DML/DMR detection
```
#R
BiocManager::install("DSS")

library(tidyr)
library(dplyr)
library(DSS)

first_file <- "./WT_data/SRX4241790_methylation_result_transfered.txt"
second_file <- "./TetTKO_data/SRR7368845_methylation_result_transfered.txt"
file_prefix <- "mm_all_chr"
file_save_path <- "./"

setwd("/home/zxy0303/project/musculus/R_analysis")
first_raw_data <- read.table(first_file, header = T)
second_raw_data <- read.table(second_file, header = T)

DSS_first_input_data <- first_raw_data %>%
    mutate(chr = paste("chr", chr, sep="")) %>%
    mutate(pos = start, N = as.numeric(methyled) + as.numeric(unmethyled), X = as.numeric(methyled)) %>%
    select(chr, pos, N, X)
DSS_second_input_data <- first_raw_data %>%
    mutate(chr = paste("chr", chr, sep="")) %>%
    mutate(pos = start, N = as.numeric(methyled) + as.numeric(unmethyled), X = as.numeric(methyled)) %>%
    select(chr, pos, N, X)

bsobj <- makeBSseqData(list(DSS_first_input_data, DSS_second_input_data), c("S1", "S2"))
dmlTest <- DMLtest(bsobj, group1 = c("S1"), group2 = c("S2"), smoothing = T)

dmls <- callDML(dmlTest, p.threshold = 0.001)
dmrs <- callDMR(dmlTest, p.threshold = 0.01)

write.table(dmlTest, paste(file_save_path, file_prefix, "_DSS_test_result.txt", sep=""), row.names = F)
write.table(dmls, paste(file_save_path, file_prefix, "_DSS_dmls_result.txt", sep = ""), row.names = F)
write.table(dmrs, paste(file_save_path, file_prefix, "_DSS_dmrs_result.txt", sep = ""), row.names = F)

#Rscript
#Rscript $HOME/project/Script/bismark_result_transfer.R ./WT_data/SRX4241790_methylation_result.txt ./TetTKO_data/SRR7368845_methylation_result.txt
```

# others
bowtie原理：
> 参考序列预处理
> Bowtie 首先对参考基因组或参考序列进行预处理，构建索引。这个过程类似于给一本书建立索引，目的是为了在后续的比对过程中能够快速地定位短读段可能出现的位置。例如，对于一个包含大量基因序列的参考基因组，Bowtie 会将其分割成小的片段（通常是固定长度），并记录每个片段的位置和相关信息。这些片段就像索引中的词条，用于快速查找。
> 基于 Burrows - Wheeler 变换（BWT）
> 构建索引的核心是 Burrows - Wheeler 变换。BWT 是一种数据变换方法，它对参考序列进行重排，使得相似的字符尽可能地聚集在一起。这种变换后的序列具有特殊的性质，能够大大提高后续查找的效率。以一个简单的字符串 “ACGTACG” 为例，经过 BWT 变换后，序列会被重排成一种更有利于查找的形式。在这个新的序列形式中，相同的碱基（如 “C”）会相对更集中，这有助于快速定位包含这些碱基的短读段在参考序列中的位置。