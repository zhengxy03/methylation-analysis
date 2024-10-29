# methylation-analysis
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
## 3.2 genome indexing
```
bismark_genome_preparation --bowtie2 ~/project/musculus/genome
```
## 3.3 read alignment
```
genome_path="$HOME/project/musculus/genome"
cd ~/project/musculus/sequence

bismark -o ./WT_mESC_rep1/bismark_result/ --parallel 4 --genome_folder ${genome_path} ./WT_mESC_rep1/*.fastq.gz
bismark -o ./TetTKO_mEC_rep1/bismark_result/ --parallel 4 --genome_folder ${genome_path} ./TetTKO_mEC_rep1/*.fastq.gz
```
