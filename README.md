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
## 1.2 experiment data
cd ../sequence
mkdir WT_mESC_rep1 TetTKO_mESC_rep1

nohup prefetch SRR736884{1..9} SRR736885{0..1} -O . &