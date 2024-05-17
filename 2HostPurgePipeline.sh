[TOC]
# 1. Install Software
## 1.1 install kraken2 v2.1.2
```
conda create -y -n kraken2
mamba create -n kraken2 -y -c bioconda kraken2 bracken krakentools krona r-optparse
kraken2 --version
#2.1.2
```
## 1.2 install kneadData v0.12.0
```
conda create -y -n kneaddata
conda activate kneaddata
kneaddata --version
#0.12.0
```
## 2. Construct the index database
### 2.1 Indexing the Kraken2
```
conda activate kraken2
mkdir -p /public/home/liuyongxin/gyy/230903camisim/rice06/database/kraken2/rice
cd /public/home/liuyongxin/gyy/230903camisim/rice06/database/kraken2/rice
##prepared genomes
##prepared names.dmp
#Download the nucleotide gd accession
time kraken2-build --download-taxonomy --threads 8 --db Oryza_sativa  --use-ftp
##it might be cost me several times based on the internet of linux, so we could download them from website and upload them to linux
##download taxonomy by myself, website: https://ftp.ncbi.nih.gov/pub/taxonomy/
#accession2taxid/nucl_gb.accession2taxid.gz
#accession2taxid/nucl_wgs.accession2taxid.gz
#taxdump.tar.gz
#mkdir Oryza_sativa/taxonomy
#gunzip nucl_gb.accession2taxid.gz Oryza_sativa/taxonomy
#gunzip nucl_wgs.accession2taxid.gz Oryza_sativa/taxonomy
#tar zxf taxdump.tar.gz  Oryza_sativa/taxonomy
##Fasta file not downloaded from NCBI may need their taxonomy information assigned explicitly.
##This can be done using the string kraken:taxid|XXX in the sequence ID.
##Such as >GWHBFPX00000001|kraken:taxid|39947  Adapter sequence
##The taxnonomy id could be found in https://www.ncbi.nlm.nih.gov/datasets/taxonomy/.
#Build the Kraken2 Database
kraken2-build --db Oryza_sativa --threads 8 --add-to-library GWHBFPX00000000.genome.fasta
#Once your library is finalized, you need to build the database. This can be done with the command(2min)
time kraken2-build --build  --db Oryza_sativa --threads 8
##Note: Sometimes, there are plenty of fasta files need modify their  taxonomy information, so we could use the following code to deal with them.
#awk '/^>/ { sub(">", "", $1); $0 = ">" $1 "|kraken:taxid|39947  Adapter sequence" } 1' GWHBFPX00000000.genome.fasta>GWHBFPX00000000.genome_1.fasta
#mv GWHBFPX00000000.genome_1.fasta GWHBFPX00000000.genome.fasta
#awk '/^>/ { sub(">", "", $1); $0 = ">" $1 "|kraken:taxid|9606  Adapter sequence" } 1' human_genome.fasta>human_genome_1.fasta
#mv human_genome_1.fasta human_genome.fasta
```
### 2.2 Indexing the KneadData
### 1)Indexing a Kneaddata v0.12.0 reference genome
```
conda activate kneaddata
cd /public/home/liuyongxin/gyy/230903camisim/rice06/database/kneaddata/rice
bowtie2-build -f GWHBFPX00000000.genome.fasta Oryza_sativa --thread 8 --seed 1
```
# 3 Run with the HostDecontamintaion Software
### 3.1 Run kraken2
```
#For your taxomyid, replace 4530 as your taxonomyid
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
mkdir -p ${path}/kraken2
mkdir -p ${path}/kraken2_qc
for i in {0..4};do
d=anonymous_read${i}
kraken2 --db /public/home/liuyongxin/gyy/230903camisim/rice06/database/kraken2/rice   \
	-1 ${path}/seq/${d}_1.fq \
	-2 ${path}/seq/${d}_2.fq \
	--threads 8 --use-names --report-zero-counts  \
	--report ${path}/kraken2/${d}.report  \
	--output ${path}/kraken2/${d}.output  
extract_kraken_reads.py -k ${path}/kraken2/${d}.output -r ${path}/kraken2/${d}.report \
	-1 ${path}/seq/${d}_1.fq \
	-2 ${path}/seq/${d}_2.fq \
	-t 4530  --include-children \
	--exclude --fastq-output \
	-o ${path}/kraken2_qc/4530remain${i}_1.fq -o2 ${path}/kraken2_qc/4530remain${i}_2.fq
```
### 3.2 Run KneadData
```
conda activate kneaddata
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
mkdir -p ${path}/kraken2_kneaddata/temp
mv ${path}/kraken2_qc/*remain* ${path}/kraken2_kneaddata/temp
file=/public/home/liuyongxin/gyy/230903camisim/rice06/database/kneaddata/rice
for i in {0..4};do
d=4530remain${i}
time memusg kneaddata -i1 ${path}/kraken2_kneaddata/temp/${d}_1.fq -i2 ${path}/kraken2_kneaddata/temp/${d}_2.fq  \
	-o  ${path}/kraken2_kneaddata \
	-db ${file}/Oryza_sativa \
	--bypass-trim  --bypass-trf \
	-t 8
done
```
