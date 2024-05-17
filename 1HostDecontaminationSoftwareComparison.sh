[TOC]

# 1. 软件安装

## 1.1 Install camisim v1.3
```
#python 3.11, 3.10 不适合camisim
conda create -n camisim python=3.7 -y
conda activate camisim
git clone https://github.com/CAMI-challenge/CAMISIM
#Here I just download zip file form https://github.com/CAMI-challenge/CAMISIM, and unzip it
unzip CAMISIM-master.zip
conda activate camisim
```
## 1.2 install memusg
```
wget https://github.com/jhclark/memusg
wget https://github.com/jhclark/memusg/blob/master/memusg
#运行memusg之前需要将其放到环境中，并检查是否可以调用
export PATH=/public/home/liuyongxin/db/soft/memusg-master:$PATH
```
## 1.3 install kneaddata v0.12.0
```
conda create -y -n kneaddata
conda activate kneaddata
kneaddata --version
#0.12.0
```
## 1.4 install bowtie2 v2.5.1
```
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-linux-x86_64.zip/download
unzip bowtie2-2.5.1-linux-x86_64.zip
cd bowtie2-2.5.1-linux-x86_64
./bowtie2 --version
export PATH="/public/home/liuyongxin/db/soft/bowtie2-2.5.1-linux-x86_64:$PATH"
source ~/.bashrc
```
## 1.5 install bwa v0.7.17
```
conda activate gaoyunyun
cd /public/home/liuyongxin/db/soft
#download bwa-0.7.17 from git-hub https://github.com/lh3/bwa
wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.17.zip
unzip v0.7.17.zip
cd bwa-0.7.17/
makes
```
## 1.6 install kraken2 v2.1.2
```
conda create -y -n kraken2
mamba create -n kraken2 -y -c bioconda kraken2 bracken krakentools krona r-optparse
kraken2 --version
#2.1.2
```
## 1.7 Install KMCP v0.9.3
```
conda activate gaoyunyun
conda install -c bioconda kmcp
kmcp --version
#v0.9.3
```
## 1.8 install KrakenUniq v1.0.4
```
conda create -y -n krakenuniq2 python=3.12.0
conda activate krakenuniq2
wget https://github.com/fbreitwieser/krakenuniq/archive/refs/tags/v1.0.4.tar.gz
tar -xzvf v1.0.4.tar.gz
cd krakenuniq-1.0.4
./install_krakenuniq.sh /public/home/liuyongxin/miniconda3/envs/krakenuniq/bin
krakenuniq --version
#1.0.4
#Remember to check your jellysish, jellyfish1 is good for krakenuniq, others not work
conda install -c bioconda jellyfish=1.1.12
```
# 2. CAMISM模拟数据
```
#修改mini_config.ini
[Main]
dataset_id=OTU
output_directory=public/home/liuyongxin/gyy/camisim/out0725
temp_directory=public/home/liuyongxin/gyy/camisim/tmp0725

[CommunityDesign]
number_of_samples=1

[community0]
metadata=public/home/liuyongxin/gyy/camisim/metadata.tsv
id_to_genome_file=public/home/liuyongxin/gyy/camisim/genome_to_id.tsv
genomes_total=1
num_real_genomes=1
```
# 3 准备数据
```
#copy camisim simulate data (fq.gz fiesll) to new files mkdir 1_ModifiedData
```
## 1)修改数据格式
```
#1_ModifiedData, added label for future mixed data
mkdir -p 1modified
 for i in {0..4}; do
    # Apply sed to modify microbiome files
    awk '{if (NR % 4 == 1) {gsub(/\/[12]/, "microbiome&")} print}' 0RawGenomes/microbiome${i}.fq > 1modified/modified_microbiome${i}.fq
    sed '1d; s/\//microbiome\//' 0RawGenomes/reads_microbiome${i}.tsv > 1modified/modified_microbiome${i}.tsv
    # Apply sed to modify host files
    awk '{if (NR % 4 == 1) {gsub(/\/[12]/, "host&")} print}' 0RawGenomes/host${i}.fq > 1modified/modified_host${i}.fq
    sed '1d; s/\//host\//' 0RawGenomes/reads_host${i}.tsv > 1modified/modified_host${i}.tsv
    cat 1modified/modified_microbiome${i}.tsv 1modified/modified_host${i}.tsv > 1modified/modified_reads${i}.tsv
done
```
## 2)检查数据大小
```
#check the row of different raw data
cd /public/home/liuyongxin/gyy/230903camisim/host06/genomes
wc -l microbiome0.fq
#79998176 microbiome0.fq
wc -l host0.fq
#719934752 host0.fq
echo "Oryza_sativa $(grep -w -c "Oryza_sativa" 3extract${seq}/extracted_reads${i}_2.tsv)" > 3extract${seq}/meta${i}.txt
echo "Pseudomonas_fluorescens $(grep -w -c "Pseudomonas_fluorescens" 3extract${seq}/extracted_reads${i}_2.tsv)" >> 3extract${seq}/meta${i}.txt
```
## 3)选择数据
### A. 水稻单菌
```
#2_SelectedData, Selected 6G data from the total data, and host:microbiomel=9:1, and merged them
seq=06
mkdir -p 2select${seq}
mkdir -p 3extract${seq}
for i in {0..4};do
    head -n "15999632" 1modified/modified_microbiome${i}.fq > 2select${seq}/selected_microbiome${i}.fq
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' 2select${seq}/selected_microbiome${i}.fq | sed 's/^@//' > 2select${seq}/selected_microbiome${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' 2select${seq}/selected_microbiome${i}.fq | sed 's/^@//' > 2select${seq}/selected_microbiome${i}_2.txt
    head -n "143996688" 1modified/modified_host${i}.fq > 2select${seq}/selected_host${i}.fq
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' 2select${seq}/selected_host${i}.fq | sed 's/^@//' > 2select${seq}/selected_host${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' 2select${seq}/selected_host${i}.fq | sed 's/^@//' > 2select${seq}/selected_host${i}_2.txt	   
    cat 2select${seq}/selected_microbiome${i}.fq 2select${seq}/selected_host${i}.fq > 3extract${seq}/anonymous_reads${i}.fq
    cat 2select${seq}/selected_microbiome${i}_1.txt 2select${seq}/selected_host${i}_1.txt> 3extract${seq}/reads${i}_1.txt
    cat 2select${seq}/selected_microbiome${i}_2.txt 2select${seq}/selected_host${i}_2.txt > 3extract${seq}/reads${i}_2.txt
    grep -F -f 3extract${seq}/reads${i}_1.txt 1modified/modified_reads${i}.tsv > 3extract${seq}/extracted_reads${i}_1.tsv
    grep -F -f 3extract${seq}/reads${i}_2.txt 1modified/modified_reads${i}.tsv > 3extract${seq}/extracted_reads${i}_2.tsv
    echo "Oryza_sativa $(grep -w -c "39947" 3extract${seq}/extracted_reads${i}_1.tsv)" > 3extract${seq}/meta${i}.txt
    echo "Pseudomonas_fluorescens $(grep -w -c "294" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
done
#2_SelectedData,6G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,6G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,6G-host:microbiomel=1:9,seq=06-3;15999632:143996688
#2_SelectedData,15G-host:microbiomel=9:1,seq=15-1;345592080:38399120
#2_SelectedData,15G-host:microbiomel=1:1,seq=15-2;191995600;191995600
#2_SelectedData,15G-host:microbiomel=1:9,seq=15-3;38399120:345592080
#2_SelectedData,30G-host:microbiomel=9:1,seq=30-1;691184160:76798240
#2_SelectedData,30G-host:microbiomel=1:1,seq=30-2;383991200:383991200
#2_SelectedData,30G-host:microbiomel=1:9,seq=30-3;76798240:691184160
```
### B. 人类单菌
```
#2_SelectedData, Selected 6G data from the total data, and host:microbiomel=9:1, and merged them
seq=06-1
mkdir -p 2select${seq}
mkdir -p 3extract${seq}
for i in {0..4};do
    head -n "15999632" 1modified/modified_microbiome${i}.fq > 2select${seq}/selected_microbiome${i}.fq
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' 2select${seq}/selected_microbiome${i}.fq | sed 's/^@//' > 2select${seq}/selected_microbiome${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' 2select${seq}/selected_microbiome${i}.fq | sed 's/^@//' > 2select${seq}/selected_microbiome${i}_2.txt
    head -n "143996688" 1modified/modified_host${i}.fq > 2select${seq}/selected_host${i}.fq
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' 2select${seq}/selected_host${i}.fq | sed 's/^@//' > 2select${seq}/selected_host${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' 2select${seq}/selected_host${i}.fq | sed 's/^@//' > 2select${seq}/selected_host${i}_2.txt	   
    cat 2select${seq}/selected_microbiome${i}.fq 2select${seq}/selected_host${i}.fq > 3extract${seq}/anonymous_reads${i}.fq
    cat 2select${seq}/selected_microbiome${i}_1.txt 2select${seq}/selected_host${i}_1.txt> 3extract${seq}/reads${i}_1.txt
    cat 2select${seq}/selected_microbiome${i}_2.txt 2select${seq}/selected_host${i}_2.txt > 3extract${seq}/reads${i}_2.txt
    grep -F -f 3extract${seq}/reads${i}_1.txt 1modified/modified_reads${i}.tsv > 3extract${seq}/extracted_reads${i}_1.tsv
    grep -F -f 3extract${seq}/reads${i}_2.txt 1modified/modified_reads${i}.tsv > 3extract${seq}/extracted_reads${i}_2.tsv
    echo "Homo_sapiens $(grep -w -c "9606" 3extract${seq}/extracted_reads${i}_1.tsv)" > 3extract${seq}/meta${i}.txt
    echo "Escherichia_coli $(grep -w -c "562" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
    done
#2_SelectedData,6G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,6G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,6G-host:microbiomel=1:9,seq=06-3;15999632:143996688
#2_SelectedData,15G-host:microbiomel=9:1,seq=15-1;345592080:38399120
#2_SelectedData,15G-host:microbiomel=1:1,seq=15-2;191995600;191995600
#2_SelectedData,15G-host:microbiomel=1:9,seq=15-3;38399120:345592080
#2_SelectedData,30G-host:microbiomel=9:1,seq=30-1;691184160:76798240
#2_SelectedData,30G-host:microbiomel=1:1,seq=30-2;383991200:383991200
#2_SelectedData,30G-host:microbiomel=1:9,seq=30-3;76798240:691184160
```
### C. 水稻多菌
```
seq=30-2
mkdir -p 2select${seq}
mkdir -p 3extract${seq}
for i in {0..4};do
    head -n "383991200" 1modified/modified_microbiome${i}.fq > 2select${seq}/selected_microbiome${i}.fq
	grep -Eo '^@S[0-9A-Za-z_-]+/1$' 2select${seq}/selected_microbiome${i}.fq | sed 's/^@//' > 2select${seq}/selected_microbiome${i}_1.txt
	grep -Eo '^@S[0-9A-Za-z_-]+/2$' 2select${seq}/selected_microbiome${i}.fq | sed 's/^@//' > 2select${seq}/selected_microbiome${i}_2.txt
    head -n "383991200" 1modified/modified_host${i}.fq > 2select${seq}/selected_host${i}.fq
	grep -Eo '^@S[0-9A-Za-z_-]+/1$' 2select${seq}/selected_host${i}.fq | sed 's/^@//' > 2select${seq}/selected_host${i}_1.txt
	grep -Eo '^@S[0-9A-Za-z_-]+/2$' 2select${seq}/selected_host${i}.fq | sed 's/^@//' > 2select${seq}/selected_host${i}_2.txt	   
	cat 2select${seq}/selected_microbiome${i}.fq 2select${seq}/selected_host${i}.fq > 3extract${seq}/anonymous_reads${i}.fq
	cat 2select${seq}/selected_microbiome${i}_1.txt 2select${seq}/selected_host${i}_1.txt> 3extract${seq}/reads${i}_1.txt
	cat 2select${seq}/selected_microbiome${i}_2.txt 2select${seq}/selected_host${i}_2.txt > 3extract${seq}/reads${i}_2.txt
	grep -F -f 3extract${seq}/reads${i}_1.txt 1modified/modified_reads${i}.tsv > 3extract${seq}/extracted_reads${i}_1.tsv
	grep -F -f 3extract${seq}/reads${i}_2.txt 1modified/modified_reads${i}.tsv > 3extract${seq}/extracted_reads${i}_2.tsv
	echo "Oryza_sativa $(grep -w -c "Oryza_sativa" 3extract${seq}/extracted_reads${i}_1.tsv)" > 3extract${seq}/meta${i}.txt
	echo "Bacillus_cereus $(grep -w -c "Bacillus_cereus" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Lactobacillus_plantarum $(grep -w -c "Lactobacillus_plantarum" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Clostridium_acetobutylicum $(grep -w -c "Clostridium_acetobutylicum" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Methylobacterium_extorquens $(grep -w -c "Methylobacterium_extorquens" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Enterobacter_cloacae $(grep -w -c "Enterobacter_cloacae" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Escherichia_coli $(grep -w -c "Escherichia_coli" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Salmonella_enterica $(grep -w -c "Salmonella_enterica" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Pseudomonas_fluorescens $(grep -w -c "Pseudomonas_fluorescens" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Pseudomonas_aeruginosa $(grep -w -c "Pseudomonas_aeruginosa" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Xanthomonas_oryzae $(grep -w -c "Xanthomonas_oryzae" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Staphylococcus_aureus $(grep -w -c "Staphylococcus_aureus" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Burkholderia_glumae_1 $(grep -w -c "Burkholderia_glumae_1" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Burkholderia_glumae_2 $(grep -w -c "Burkholderia_glumae_2" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Burkholderia_glumae_3 $(grep -w -c "Burkholderia_glumae_3" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Burkholderia_glumae_4 $(grep -w -c "Burkholderia_glumae_4" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Pyricularia_grisea_1 $(grep -w -c "Pyricularia_grisea_1" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Pyricularia_grisea_2 $(grep -w -c "Pyricularia_grisea_2" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Pyricularia_grisea_3 $(grep -w -c "Pyricularia_grisea_3" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Fusarium_graminearum_1 $(grep -w -c "Fusarium_graminearum_1" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Fusarium_graminearum_2 $(grep -w -c "Fusarium_graminearum_2" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Fusarium_graminearum_3 $(grep -w -c "Fusarium_graminearum_3" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
	echo "Fusarium_graminearum_4 $(grep -w -c "Fusarium_graminearum_4" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
done
```
### D. 人类多菌

## 4)提取数据为双端
```
#3_ExtractedData.sh
seq=06-1
for i in {0..4};do
	# Separate reads starting with "/1"
	seqkit grep -f 3extract${seq}/reads${i}_1.txt -p "/1" 3extract${seq}/anonymous_reads${i}.fq > 3extract${seq}/anonymous_read${i}_1.fq
	# Separate reads starting with "/2"
	seqkit grep -f 3extract${seq}/reads${i}_2.txt -p "/2" 3extract${seq}/anonymous_reads${i}.fq > 3extract${seq}/anonymous_read${i}_2.fq
done
mkdir seq6-1/seq
mv ../genomes/3extract06/anonymous_read*_1.fq ./seq/
mv ../genomes/3extract06/anonymous_read*_2.fq ./seq/
cp ../genomes/0RawGenomes/test.sh ./
cp ../genomes/0RawGenomes/00sbatch_microbiome.sh ./ 
```
# 4. 软件测评
```
export PATH=/public/home/liuyongxin/db/soft/memusg-master:$PATH 
```
## 4.1 Kneaddata v0.12.0测试
### 1)Indexing a Kneaddata v0.12.0 reference genome
```
source activate kneaddata
cd /public/home/liuyongxin/gyy/230903camisim/rice06/database/kneaddata/rice
time  memusg bowtie2-build -f GWHBFPX00000000.genome.fasta Oryza_sativa --thread 8 --seed 1
```
### 2)Running Kneaddata v0.12.0
```
source activate kneaddata
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
file=/public/home/liuyongxin/gyy/230903camisim/rice06/database/kneaddata/rice
for i in {0..4};do
	d=anonymous_read${i}
	time memusg kneaddata -i1 ${path}/seq/${d}_1.fq -i2 ${path}/seq/${d}_2.fq  \
	-o  ${path}/kneaddata_output \
	-db ${file}/Oryza_sativa \
	--bypass-trim  --bypass-trf\
	-t 8
done
```
### 3)质控结果汇总(kneaddata)
```
conda activate kneaddata
kneaddata_read_count_table --input kneaddata_output/ --output kneaddata_output/kneaddata.txt
```
### 4)计算获取的reads数量
```
for i in {0..4};do
	seqkit stat "kneaddata_output/anonymous_read${i}_1_kneaddata_paired_1.fastq" -a >> kneaddata_output/seqkit.txt
done
```
### 5)去宿主数据的准确性kneaddata
#### STEP1： mv reads lable of reads1\&reads2 from genomes/3extract06 directory
```
conda activate kneaddata
rm -rf kneaddata_output/*unmatched*
#Check the row of R1_read.txt&R2_read.txt
wc -l doc/extracted_reads0_1.tsv
wc -l doc/extracted_reads0_2.tsv
mkdir -p kneaddata_output/f1_score
```
#### STEP2：extract the label of all fastq files
```
#extract the label of all fastq files
for i in {0..4};do 
	grep -Eo '^@S[0-9A-Za-z_-]+/1$' kneaddata_output/anonymous_read${i}_1_kneaddata_paired_1.fastq | sed 's/^@//' > kneaddata_output/f1_score/paired${i}_1.txt
	grep -Eo '^@S[0-9A-Za-z_-]+/2$' kneaddata_output/anonymous_read${i}_1_kneaddata_paired_2.fastq | sed 's/^@//' > kneaddata_output/f1_score/paired${i}_2.txt
	grep -Eo '^@S[0-9A-Za-z_-]+/1$' kneaddata_output/anonymous_read${i}_1_kneaddata_Oryza_sativa_bowtie2_paired_contam_1.fastq | sed 's/^@//' > kneaddata_output/f1_score/bowtie2_paired_contam${i}_1.txt
	grep -Eo '^@S[0-9A-Za-z_-]+/2$' kneaddata_output/anonymous_read${i}_1_kneaddata_Oryza_sativa_bowtie2_paired_contam_2.fastq | sed 's/^@//' > kneaddata_output/f1_score/bowtie2_paired_contam${i}_2.txt
done
```
#### STEP3：Prepare different values
```
#filter and separate "Oryza_sativa" and "Pseudomonas_fluorescens"
for i in {0..4}; do
    awk '$2 == "Oryza_sativa" { print > "kneaddata_output/f1_score/Oryza_sativa'"$i"'_R1.txt" } $2 == "Pseudomonas_fluorescens" { print > "kneaddata_output/f1_score/Pseudomonas_fluorescens'"$i"'_R1.txt" }' doc/extracted_reads${i}_1.tsv
    awk '$2 == "Oryza_sativa" { print > "kneaddata_output/f1_score/Oryza_sativa'"$i"'_R2.txt" } $2 == "Pseudomonas_fluorescens" { print > "kneaddata_output/f1_score/Pseudomonas_fluorescens'"$i"'_R2.txt" }' doc/extracted_reads${i}_2.tsv
    # Calculate TP in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/bowtie2_paired_contam${i}_1.txt kneaddata_output/f1_score/Oryza_sativa${i}_R1.txt > kneaddata_output/f1_score/TP${i}_R1.txt
    # Calculate TP in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/bowtie2_paired_contam${i}_2.txt kneaddata_output/f1_score/Oryza_sativa${i}_R2.txt > kneaddata_output/f1_score/TP${i}_R2.txt
    # Calculate FP in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/bowtie2_paired_contam${i}_1.txt kneaddata_output/f1_score/Pseudomonas_fluorescens${i}_R1.txt >kneaddata_output/f1_score/FP${i}_R1.txt
    # Calculate FP in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/bowtie2_paired_contam${i}_2.txt kneaddata_output/f1_score/Pseudomonas_fluorescens${i}_R2.txt > kneaddata_output/f1_score/FP${i}_R2.txt
    # Calculate FN in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/paired${i}_1.txt kneaddata_output/f1_score/Oryza_sativa${i}_R1.txt > kneaddata_output/f1_score/FN${i}_R1.txt
    # Calculate FN in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/paired${i}_2.txt kneaddata_output/f1_score/Oryza_sativa${i}_R2.txt > kneaddata_output/f1_score/FN${i}_R2.txt
    # Calculate TN in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/paired${i}_1.txt kneaddata_output/f1_score/Pseudomonas_fluorescens${i}_R1.txt > kneaddata_output/f1_score/TN${i}_R1.txt
    # Calculate TN in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' kneaddata_output/f1_score/paired${i}_2.txt kneaddata_output/f1_score/Pseudomonas_fluorescens${i}_R2.txt > kneaddata_output/f1_score/TN${i}_R2.txt
done
```
#### STEP4：Calculate different values
```
for i in {0..4}; do
    # Calculate TP, FP, FN, and TN totals
    total_tp_r1=$(wc -l < kneaddata_output/f1_score/"TP${i}_R1.txt")
    total_tp_r2=$(wc -l < kneaddata_output/f1_score/"TP${i}_R2.txt")
    total_tp=$(($total_tp_r1 + $total_tp_r2))
    total_fp_r1=$(wc -l < kneaddata_output/f1_score/"FP${i}_R1.txt")
    total_fp_r2=$(wc -l < kneaddata_output/f1_score/"FP${i}_R2.txt")
    total_fp=$(($total_fp_r1 + $total_fp_r2))
    total_fn_r1=$(wc -l < kneaddata_output/f1_score/"FN${i}_R1.txt")
    total_fn_r2=$(wc -l < kneaddata_output/f1_score/"FN${i}_R2.txt")
    total_fn=$(($total_fn_r1 + $total_fn_r2))
    total_tn_r1=$(wc -l < kneaddata_output/f1_score/"TN${i}_R1.txt")
    total_tn_r2=$(wc -l < kneaddata_output/f1_score/"TN${i}_R2.txt")
    total_tn=$(($total_tn_r1 + $total_tn_r2))
    # Output the results to result.txt
    echo "TP${i}  $total_tp" > kneaddata_output/f1_score/result${i}.txt
    echo "FP${i}  $total_fp" >> kneaddata_output/f1_score/result${i}.txt
    echo "FN${i}  $total_fn" >> kneaddata_output/f1_score/result${i}.txt
    echo "TN${i}  $total_tn" >> kneaddata_output/f1_score/result${i}.txt
    # Calculate precision, recall, accuracy, and F1 score
    precision=$(echo "scale=4; $total_tp / ($total_tp + $total_fp)" | bc)
    recall=$(echo "scale=4; $total_tp / ($total_tp + $total_fn)" | bc)
    accuracy=$(echo "scale=4; ($total_tp + $total_tn) / ($total_tp + $total_fp + $total_fn + $total_tn)" | bc)
    f1_score=$(echo "scale=4; 2 * $total_tp / (2 * $total_tp + $total_fp + $total_fn)" | bc)
    # Output the results to result.txt
    echo "Precision${i}  $precision" >> kneaddata_output/f1_score/result${i}.txt
    echo "Recall${i}  $recall" >> kneaddata_output/f1_score/result${i}.txt
    echo "Accuracy${i}  $accuracy" >> kneaddata_output/f1_score/result${i}.txt
    echo "F1 Score${i}  $f1_score" >> kneaddata_output/f1_score/result${i}.txt
done
```
## 4.2 运行BWA

### 1)Indexing a BWA v0.7.17 reference genome
```
mkdir -p /public/home/liuyongxin/gyy/230903camisim/rice06/database/bwa/rice
cd /public/home/liuyongxin/gyy/230903camisim/rice06/database/bwa/rice
time  memusg /public/home/liuyongxin/db/soft/bwa-0.7.17/bwa index GWHBFPX00000000.genome.fasta
```
### 2)Running v0.7.17
```
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
file=/public/home/liuyongxin/gyy/230903camisim/rice06/database/bwa/rice
mkdir -p ${path}/bwa
for i in {0..4};do
    d=anonymous_read${i}
    time memusg /public/home/liuyongxin/db/soft/bwa-0.7.17/bwa mem -t 8 \
    ${file}/GWHBFPX00000000.genome.fasta \
    ${path}/seq/${d}_1.fq ${path}/seq/${d}_2.fq | samtools sort -o ${path}/bwa/bwaoutput${i}.sorted.bam
    # Index the sorted BAM file 
    samtools index  ${path}/bwa/bwaoutput${i}.sorted.bam
    # Extract aligned reads in BAM format
    samtools view -b -F 4 ${path}/bwa/bwaoutput${i}.sorted.bam > ${path}/bwa/aligned_reads${i}.bam
    samtools view -b -f 4 ${path}/bwa/bwaoutput${i}.sorted.bam > ${path}/bwa/not_aligned_reads${i}.bam
    # Convert the aligned reads BAM to FASTQ format
    samtools fastq ${path}/bwa/aligned_reads${i}.bam > ${path}/bwa/aligned_reads${i}.fq
    samtools fastq ${path}/bwa/not_aligned_reads${i}.bam > ${path}/bwa/not_aligned_reads${i}.fq
done
```
### 3)准备需要计算的数据
```
mkdir -p bwa/f1_score
mkdir -p bwa/temp
mv bwa/*bam* bwa/temp
mv bwa/*.fq bwa/f1_score
```
### 4)去宿主数据的准确性BWA
#### STEP1：extract the label of all fastq files
```
for i in {0..4};do
    #extract the label of all fastq files
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' bwa/f1_score/aligned_reads${i}.fq | sed 's/^@//' > bwa/f1_score/contam${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' bwa/f1_score/aligned_reads${i}.fq | sed 's/^@//' > bwa/f1_score/contam${i}_2.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' bwa/f1_score/not_aligned_reads${i}.fq | sed 's/^@//' > bwa/f1_score/paired${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' bwa/f1_score/not_aligned_reads${i}.fq | sed 's/^@//' > bwa/f1_score/paired${i}_2.txt
done
```
#### STEP2：Prepare different values
```
for i in {0..4};do
    #  filter and separate "Oryza_sativa" and "Pseudomonas_fluorescens"
    awk '$2 == "Oryza_sativa" { print > "bwa/f1_score/Oryza_sativa'"$i"'_R1.txt" } $2 == "Pseudomonas_fluorescens" { print > "bwa/f1_score/Pseudomonas_fluorescens'"$i"'_R1.txt" }' doc/extracted_reads${i}_1.tsv
    awk '$2 == "Oryza_sativa" { print > "bwa/f1_score/Oryza_sativa'"$i"'_R2.txt" } $2 == "Pseudomonas_fluorescens" { print > "bwa/f1_score/Pseudomonas_fluorescens'"$i"'_R2.txt" }' doc/extracted_reads${i}_2.tsv
    # Calculate TP in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/contam${i}_1.txt bwa/f1_score/Oryza_sativa${i}_R1.txt > bwa/f1_score/TP${i}_R1.txt
    # Calculate TP in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/contam${i}_2.txt bwa/f1_score/Oryza_sativa${i}_R2.txt > bwa/f1_score/TP${i}_R2.txt
    # Calculate FP in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/contam${i}_1.txt bwa/f1_score/Pseudomonas_fluorescens${i}_R1.txt > bwa/f1_score/FP${i}_R1.txt
    # Calculate FP in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/contam${i}_2.txt bwa/f1_score/Pseudomonas_fluorescens${i}_R2.txt > bwa/f1_score/FP${i}_R2.txt
    # Calculate FN in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/paired${i}_1.txt bwa/f1_score/Oryza_sativa${i}_R1.txt > bwa/f1_score/FN${i}_R1.txt
    # Calculate FN in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/paired${i}_2.txt bwa/f1_score/Oryza_sativa${i}_R2.txt > bwa/f1_score/FN${i}_R2.txt
    # Calculate TN in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/paired${i}_1.txt bwa/f1_score/Pseudomonas_fluorescens${i}_R1.txt > bwa/f1_score/TN${i}_R1.txt
    # Calculate TN in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bwa/f1_score/paired${i}_2.txt bwa/f1_score/Pseudomonas_fluorescens${i}_R2.txt > bwa/f1_score/TN${i}_R2.txt
done
```
#### STEP3：Calculate different values
```
for i in {0..4}; do
    # Calculate TP, FP, FN, and TN totals
    total_tp_r1=$(wc -l < bwa/f1_score/"TP${i}_R1.txt")
    total_tp_r2=$(wc -l < bwa/f1_score/"TP${i}_R2.txt")
    total_tp=$(($total_tp_r1 + $total_tp_r2))
    total_fp_r1=$(wc -l < bwa/f1_score/"FP${i}_R1.txt")
    total_fp_r2=$(wc -l < bwa/f1_score/"FP${i}_R2.txt")
    total_fp=$(($total_fp_r1 + $total_fp_r2))
    total_fn_r1=$(wc -l < bwa/f1_score/"FN${i}_R1.txt")
    total_fn_r2=$(wc -l < bwa/f1_score/"FN${i}_R2.txt")
    total_fn=$(($total_fn_r1 + $total_fn_r2))
    total_tn_r1=$(wc -l < bwa/f1_score/"TN${i}_R1.txt")
    total_tn_r2=$(wc -l < bwa/f1_score/"TN${i}_R2.txt")
    total_tn=$(($total_tn_r1 + $total_tn_r2))
    # Output the results to result.txt
    echo "TP${i}  $total_tp" > bwa/f1_score/result${i}.txt
    echo "FP${i}  $total_fp" >> bwa/f1_score/result${i}.txt
    echo "FN${i}  $total_fn" >> bwa/f1_score/result${i}.txt
    echo "TN${i}  $total_tn" >> bwa/f1_score/result${i}.txt
    # Calculate precision, recall, accuracy, and F1 score
    precision=$(echo "scale=4; $total_tp / ($total_tp + $total_fp)" | bc)
    recall=$(echo "scale=4; $total_tp / ($total_tp + $total_fn)" | bc)
    accuracy=$(echo "scale=4; ($total_tp + $total_tn) / ($total_tp + $total_fp + $total_fn + $total_tn)" | bc)
    f1_score=$(echo "scale=4; 2 * $total_tp / (2 * $total_tp + $total_fp + $total_fn)" | bc)
    # Output the results to result.txt
    echo "Precision${i}  $precision" >> bwa/f1_score/result${i}.txt
    echo "Recall${i}  $recall" >> bwa/f1_score/result${i}.txt
    echo "Accuracy${i}  $accuracy" >> bwa/f1_score/result${i}.txt
    echo "F1 Score${i}  $f1_score" >> bwa/f1_score/result${i}.txt
done
```
## 4.3 运行Bowtie2 v2.5.1
### 1)Indexing a Bowtie2 v2.5.1 reference genome
```
mkdir -p /public/home/liuyongxin/gyy/230903camisim/rice06/database/bowtie2/rice
cd /public/home/liuyongxin/gyy/230903camisim/rice06/database/bowtie2/rice
cp ../../bwa/rice/GWHBFPX00000000.genome.fasta ./
time memusg /public/home/liuyongxin/db/soft/bowtie2-2.5.1-linux-x86_64/bowtie2-build ./GWHBFPX00000000.genome.fasta rice
```
### 2)Running Bowtie2 v2.5.1
```
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
file=/public/home/liuyongxin/db/soft/bowtie2-2.5.1-linux-x86_64/
mkdir -p ${path}/bowtie2
for i in {0..4};do
    d=anonymous_read${i}
    time memusg /public/home/liuyongxin/db/soft/bowtie2-2.5.1-linux-x86_64/bowtie2 -x ${file}/rice \
    -1 ${path}/seq/${d}_1.fq -2 ${path}/seq/${d}_2.fq \
    -S ${path}/bowtie2/bowtie2${i}.sam	
    # Convert the BAM file to a sorted BAM file.
    samtools sort ${path}/bowtie2/bowtie2${i}.sam -o ${path}/bowtie2/bowtie2${i}.sorted.bam
    # Index the sorted BAM file 
    samtools index ${path}/bowtie2/bowtie2${i}.sorted.bam
    # Extract aligned reads in BAM format
    samtools view -b -F 4 ${path}/bowtie2/bowtie2${i}.sorted.bam > ${path}/bowtie2/aligned_reads${i}.bam
    samtools view -b -f 4 ${path}/bowtie2/bowtie2${i}.sorted.bam > ${path}/bowtie2/not_aligned_reads${i}.bam
    # Convert the aligned reads BAM to FASTQ format
    samtools fastq ${path}/bowtie2/aligned_reads${i}.bam > ${path}/bowtie2/aligned_reads${i}.fq
    samtools fastq ${path}/bowtie2/not_aligned_reads${i}.bam > ${path}/bowtie2/not_aligned_reads${i}.fq
done
```
### 3)准备需要计算的数据
```
mkdir -p bowtie2/f1_score
mkdir -p bowtie2/temp
mv bowtie2/*bam* bowtie2/temp
mv bowtie2/*sam* bowtie2/temp
mv bowtie2/*.fq bowtie2/f1_score
```
### 4)去宿主数据的准确性Bowtie2
#### STEP1：extract the label of all fastq files
```
for i in {0..4};do
    #extract the label of all fastq files
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' bowtie2/f1_score/aligned_reads${i}.fq | sed 's/^@//' > bowtie2/f1_score/contam${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' bowtie2/f1_score/aligned_reads${i}.fq | sed 's/^@//' > bowtie2/f1_score/contam${i}_2.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/1$' bowtie2/f1_score/not_aligned_reads${i}.fq | sed 's/^@//' > bowtie2/f1_score/paired${i}_1.txt
    grep -Eo '^@S[0-9A-Za-z_-]+/2$' bowtie2/f1_score/not_aligned_reads${i}.fq | sed 's/^@//' > bowtie2/f1_score/paired${i}_2.txt
done
```
#### STEP2：Prepare different values
```
for i in {0..4};do
    #  filter and separate "Oryza_sativa" and "Pseudomonas_fluorescens"
	awk '$2 == "Oryza_sativa" { print > "bowtie2/f1_score/Oryza_sativa'"$i"'_R1.txt" } $2 == "Pseudomonas_fluorescens" { print > "bowtie2/f1_score/Pseudomonas_fluorescens'"$i"'_R1.txt" }' doc/extracted_reads${i}_1.tsv
    awk '$2 == "Oryza_sativa" { print > "bowtie2/f1_score/Oryza_sativa'"$i"'_R2.txt" } $2 == "Pseudomonas_fluorescens" { print > "bowtie2/f1_score/Pseudomonas_fluorescens'"$i"'_R2.txt" }' doc/extracted_reads${i}_2.tsv
    # Calculate TP in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/contam${i}_1.txt bowtie2/f1_score/Oryza_sativa${i}_R1.txt > bowtie2/f1_score/TP${i}_R1.txt
    # Calculate TP in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/contam${i}_2.txt bowtie2/f1_score/Oryza_sativa${i}_R2.txt > bowtie2/f1_score/TP${i}_R2.txt
    # Calculate FP in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/contam${i}_1.txt bowtie2/f1_score/Pseudomonas_fluorescens${i}_R1.txt > bowtie2/f1_score/FP${i}_R1.txt
    # Calculate FP in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/contam${i}_2.txt bowtie2/f1_score/Pseudomonas_fluorescens${i}_R2.txt > bowtie2/f1_score/FP${i}_R2.txt
    # Calculate FN in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/paired${i}_1.txt bowtie2/f1_score/Oryza_sativa${i}_R1.txt > bowtie2/f1_score/FN${i}_R1.txt
    # Calculate FN in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/paired${i}_2.txt bowtie2/f1_score/Oryza_sativa${i}_R2.txt > bowtie2/f1_score/FN${i}_R2.txt
    # Calculate TN in R1
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/paired${i}_1.txt bowtie2/f1_score/Pseudomonas_fluorescens${i}_R1.txt > bowtie2/f1_score/TN${i}_R1.txt
    # Calculate TN in R2
    awk 'FNR==NR{a[$1]; next} $1 in a' bowtie2/f1_score/paired${i}_2.txt bowtie2/f1_score/Pseudomonas_fluorescens${i}_R2.txt > bowtie2/f1_score/TN${i}_R2.txt
done
```
#### STEP3：Calculate different values
```
for i in {0..4}; do
    # Calculate TP, FP, FN, and TN totals
    total_tp_r1=$(wc -l < bowtie2/f1_score/"TP${i}_R1.txt")
    total_tp_r2=$(wc -l < bowtie2/f1_score/"TP${i}_R2.txt")
    total_tp=$(($total_tp_r1 + $total_tp_r2))
    total_fp_r1=$(wc -l < bowtie2/f1_score/"FP${i}_R1.txt")
    total_fp_r2=$(wc -l < bowtie2/f1_score/"FP${i}_R2.txt")
    total_fp=$(($total_fp_r1 + $total_fp_r2))
    total_fn_r1=$(wc -l < bowtie2/f1_score/"FN${i}_R1.txt")
    total_fn_r2=$(wc -l < bowtie2/f1_score/"FN${i}_R2.txt")
    total_fn=$(($total_fn_r1 + $total_fn_r2))
    total_tn_r1=$(wc -l < bowtie2/f1_score/"TN${i}_R1.txt")
    total_tn_r2=$(wc -l < bowtie2/f1_score/"TN${i}_R2.txt")
    total_tn=$(($total_tn_r1 + $total_tn_r2))
    # Output the results to result.txt
    echo "TP${i}  $total_tp" > bowtie2/f1_score/result${i}.txt
    echo "FP${i}  $total_fp" >> bowtie2/f1_score/result${i}.txt
    echo "FN${i}  $total_fn" >> bowtie2/f1_score/result${i}.txt
    echo "TN${i}  $total_tn" >> bowtie2/f1_score/result${i}.txt
    # Calculate precision, recall, accuracy, and F1 score
    precision=$(echo "scale=4; $total_tp / ($total_tp + $total_fp)" | bc)
    recall=$(echo "scale=4; $total_tp / ($total_tp + $total_fn)" | bc)
    accuracy=$(echo "scale=4; ($total_tp + $total_tn) / ($total_tp + $total_fp + $total_fn + $total_tn)" | bc)
    f1_score=$(echo "scale=4; 2 * $total_tp / (2 * $total_tp + $total_fp + $total_fn)" | bc)
    # Output the results to result.txt
    echo "Precision${i}  $precision" >> bowtie2/f1_score/result${i}.txt
    echo "Recall${i}  $recall" >> bowtie2/f1_score/result${i}.txt
    echo "Accuracy${i}  $accuracy" >> bowtie2/f1_score/result${i}.txt
    echo "F1 Score${i}  $f1_score" >> bowtie2/f1_score/result${i}.txt
done
```
## 4.4 运行kraken2 v2.1.2

### 1)构建kraken2 v2.1.2的数据
```
conda activate kraken2
mkdir -p /public/home/liuyongxin/gyy/230903camisim/rice06/database/kraken2/rice
cd /public/home/liuyongxin/gyy/230903camisim/rice06/database/kraken2/rice
##prepared genomes
##prepared names.dmp
#Download the nucleotide gd accession(depend on the internet)
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
#Build the Kraken2 Database
time memusg kraken2-build --db Oryza_sativa --threads 8 --add-to-library GWHBFPX00000000.genome.fasta
#Once your library is finalized, you need to build the database. This can be done with the command(2min)
time kraken2-build --build  --db Oryza_sativa --threads 8
##Note: Sometimes, there are plenty of fasta files need modify their  taxonomy information, so we could use the following code to deal with them.
#awk '/^>/ { sub(">", "", $1); $0 = ">" $1 "|kraken:taxid|39947  Adapter sequence" } 1' GWHBFPX00000000.genome.fasta>GWHBFPX00000000.genome_1.fasta
#mv GWHBFPX00000000.genome_1.fasta GWHBFPX00000000.genome.fasta
#awk '/^>/ { sub(">", "", $1); $0 = ">" $1 "|kraken:taxid|9606  Adapter sequence" } 1' human_genome.fasta>human_genome_1.fasta
#mv human_genome_1.fasta human_genome.fasta
```
### 2)运行kraken2 v2.1.2
```
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
mkdir -p ${path}/kraken2
mkdir -p ${path}/kraken2_qc
for i in {0..4};do
d=anonymous_read${i}
time memusg kraken2 --db /public/home/liuyongxin/gyy/230903camisim/rice06/database/kraken2/rice   \
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
extract_kraken_reads.py -k ${path}/kraken2/${d}.output -r ${path}/kraken2/${d}.report \
        -1 ${path}/seq/${d}_1.fq \
        -2 ${path}/seq/${d}_2.fq \
        -t 4530  --include-children --fastq-output \
        -o ${path}/kraken2_qc/4530${i}_1.fq -o2 ${path}/kraken2_qc/4530${i}_2.fq
done
#For human data, replace 4530 as 9606
```
### 3)准备需要计算的数据
```
mkdir -p kraken2_qc/f1_score
mkdir -p kraken2_qc/temp
mv kraken2_qc/*.fq kraken2_qc/temp
```
### 4)去宿主数据的准确性kraken2
#### STEP1：extract the label of all fastq files
```
 for i in {0..4};do
    #extract the label of all fastq files
	grep -Eo '^@[0-9A-Za-z_-]+/1$' kraken2_qc/temp/4530${i}_1.fq | sed 's/^@//' > kraken2_qc/4530read${i}_1.txt
	grep -Eo '^@[0-9A-Za-z_-]+/2$' kraken2_qc/temp/4530${i}_2.fq | sed 's/^@//' > kraken2_qc/4530read${i}_2.txt
	grep -Eo '^@[0-9A-Za-z_-]+/1$' kraken2_qc/temp/4530remain${i}_1.fq | sed 's/^@//' > kraken2_qc/4530remain_read${i}_1.txt
	grep -Eo '^@[0-9A-Za-z_-]+/2$' kraken2_qc/temp/4530remain${i}_2.fq | sed 's/^@//' > kraken2_qc/4530remain_read${i}_2.txt
done
```
#### STEP2：Prepare different values
```
for i in {0..4};do
    #  filter and separate "Oryza_sativa" and "Pseudomonas_fluorescens"
	awk '$2 == "Oryza_sativa" { print > "kraken2_qc/f1_score/Oryza_sativa'"$i"'_R1.txt" } $2 == "Pseudomonas_fluorescens" { print > "kraken2_qc/f1_score/Pseudomonas_fluorescens'"$i"'_R1.txt" }' doc/extracted_reads${i}_1.tsv
    awk '$2 == "Oryza_sativa" { print > "kraken2_qc/f1_score/Oryza_sativa'"$i"'_R2.txt" } $2 == "Pseudomonas_fluorescens" { print > "kraken2_qc/f1_score/Pseudomonas_fluorescens'"$i"'_R2.txt" }' doc/extracted_reads${i}_2.tsv
	# Calculate TP in R1
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530read${i}_1.txt kraken2_qc/f1_score/Oryza_sativa${i}_R1.txt > kraken2_qc/f1_score/TP${i}_R1.txt
	# Calculate TP in R2
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530read${i}_2.txt kraken2_qc/f1_score/Oryza_sativa${i}_R2.txt > kraken2_qc/f1_score/TP${i}_R2.txt
	# Calculate FP in R1
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530read${i}_1.txt kraken2_qc/f1_score/Pseudomonas_fluorescens${i}_R1.txt > kraken2_qc/f1_score/FP${i}_R1.txt
	# Calculate FP in R2
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530read${i}_2.txt kraken2_qc/f1_score/Pseudomonas_fluorescens${i}_R2.txt > kraken2_qc/f1_score/FP${i}_R2.txt
	# Calculate FN in R1
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530remain_read${i}_1.txt kraken2_qc/f1_score/Oryza_sativa${i}_R1.txt  > kraken2_qc/f1_score/FN${i}_R1.txt
	# Calculate FN in R2
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530remain_read${i}_2.txt kraken2_qc/f1_score/Oryza_sativa${i}_R2.txt > kraken2_qc/f1_score/FN${i}_R2.txt		
	# Calculate TN in R1
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530remain_read${i}_1.txt kraken2_qc/f1_score/Pseudomonas_fluorescens${i}_R1.txt > kraken2_qc/f1_score/TN${i}_R1.txt
	# Calculate TN in R2
	awk 'FNR==NR{a[$1]; next} $1 in a' kraken2_qc/4530remain_read${i}_2.txt kraken2_qc/f1_score/Pseudomonas_fluorescens${i}_R2.txt > kraken2_qc/f1_score/TN${i}_R2.txt
done
```
#### STEP3：Calculate different values
```
for i in {0..4}; do
    # Calculate TP, FP, FN, and TN totals
    total_tp_r1=$(wc -l < kraken2_qc/f1_score/"TP${i}_R1.txt")
    total_tp_r2=$(wc -l < kraken2_qc/f1_score/"TP${i}_R2.txt")
    total_tp=$(($total_tp_r1 + $total_tp_r2))
    total_fp_r1=$(wc -l < kraken2_qc/f1_score/"FP${i}_R1.txt")
    total_fp_r2=$(wc -l < kraken2_qc/f1_score/"FP${i}_R2.txt")
    total_fp=$(($total_fp_r1 + $total_fp_r2))
    total_fn_r1=$(wc -l < kraken2_qc/f1_score/"FN${i}_R1.txt")
    total_fn_r2=$(wc -l < kraken2_qc/f1_score/"FN${i}_R2.txt")
    total_fn=$(($total_fn_r1 + $total_fn_r2))
    total_tn_r1=$(wc -l < kraken2_qc/f1_score/"TN${i}_R1.txt")
    total_tn_r2=$(wc -l < kraken2_qc/f1_score/"TN${i}_R2.txt")
    total_tn=$(($total_tn_r1 + $total_tn_r2))
    # Output the results to result.txt
    echo "TP${i}  $total_tp" > kraken2_qc/f1_score/result${i}.txt
    echo "FP${i}  $total_fp" >> kraken2_qc/f1_score/result${i}.txt
    echo "FN${i}  $total_fn" >> kraken2_qc/f1_score/result${i}.txt
    echo "TN${i}  $total_tn" >> kraken2_qc/f1_score/result${i}.txt
    # Calculate precision, recall, accuracy, and F1 score
    precision=$(echo "scale=4; $total_tp / ($total_tp + $total_fp)" | bc)
    recall=$(echo "scale=4; $total_tp / ($total_tp + $total_fn)" | bc)
    accuracy=$(echo "scale=4; ($total_tp + $total_tn) / ($total_tp + $total_fp + $total_fn + $total_tn)" | bc)
    f1_score=$(echo "scale=4; 2 * $total_tp / (2 * $total_tp + $total_fp + $total_fn)" | bc)
    # Output the results to result.txt
    echo "Precision${i}  $precision" >> kraken2_qc/f1_score/result${i}.txt
    echo "Recall${i}  $recall" >> kraken2_qc/f1_score/result${i}.txt
    echo "Accuracy${i}  $accuracy" >> kraken2_qc/f1_score/result${i}.txt
    echo "F1 Score${i}  $f1_score" >> kraken2_qc/f1_score/result${i}.txt
done
```
## 4.5 运行kmcp v0.9.3

### 1)构建kmcp v0.9.3的数据
```
conda activate gaoyunyun
cd /public/home/liuyongxin/gyy/230903camisim/rice06/database/kmcp/rice
#compute k-mers
time memusg kmcp compute -k 21 --split-number 10 --split-overlap 150 \
    --in-dir genomes/ --out-dir genomes-k21-n10
# index k-mers
kmcp index --false-positive-rate 0.1 --num-hash 1 --in-dir genomes-k21-n10/ --out-dir genomes.kmcp
# delete temporary files
rm -rf genomes-k21-n10/
#如果提交任务建库无法完成index这步，只能在计算节点上计算（速度挺快的,5min）
```
### 2)运行kmcp
```
#需提前准备name.map文件
source activate gaoyunyun
for i in {0..4};do    
    time memusg kmcp search --db-dir /public/home/liuyongxin/db/soft/kmcp/genomes.kmcp/ -1 seq/anonymous_read${i}_1.fq -2 seq/anonymous_read${i}_2.fq --out-file kmcp/read${i}.tsv.gz --load-whole-db
done
```
### 3)准备需要计算的数据
```
mkdir -p kmcp/f1_score
mkdir -p kmcp/temp
gunzip kmcp/read*.gz
mv kmcp/read* kmcp/f1_score/
```
### 4)去宿主数据的准确性KMCP
#### STEP1：extract the label of all fastq files
```
for i in {0..4};do
    awk '!/^#/ {print $1}' kmcp/f1_score/result${i}.tsv | sort | uniq > kmcp/f1_score/contam${i}.txt
    awk 'NR==FNR{labels[$1]; next} !($1 in labels){print $1}' kmcp/f1_score/contam${i}.txt doc/extracted_reads${i}_1.tsv | cut -f 1 > kmcp/f1_score/paired${i}.txt
done
```
#### STEP2：Prepare different values
```
for i in {0..4};do
    #  filter and separate "Oryza_sativa" and "Pseudomonas_fluorescens"
    awk '$2 == "Oryza_sativa" { print > "kmcp/f1_score/Oryza_sativa'"$i"'_R1.txt" } $2 == "Pseudomonas_fluorescens" { print > "kmcp/f1_score/Pseudomonas_fluorescens'"$i"'_R1.txt" }' doc/extracted_reads${i}_1.tsv
    # Calculate TP
    awk 'FNR==NR{a[$1]; next} $1 in a' kmcp/f1_score/contam${i}.txt kmcp/f1_score/Oryza_sativa${i}_R1.txt > kmcp/f1_score/TP${i}_R1.txt
    # Calculate FP
    awk 'FNR==NR{a[$1]; next} $1 in a' kmcp/f1_score/contam${i}.txt kmcp/f1_score/Pseudomonas_fluorescens${i}_R1.txt > kmcp/f1_score/FP${i}_R1.txt
    # Calculate FN
    awk 'FNR==NR{a[$1]; next} $1 in a' kmcp/f1_score/paired${i}.txt kmcp/f1_score/Oryza_sativa${i}_R1.txt > kmcp/f1_score/FN${i}_R1.txt
    # Calculate TN
    awk 'FNR==NR{a[$1]; next} $1 in a' kmcp/f1_score/paired${i}.txt kmcp/f1_score/Pseudomonas_fluorescens${i}_R1.txt > kmcp/f1_score/TN${i}_R1.txt
done

```
#### STEP3：Calculate different values
```
for i in {0..4}; do
    # Calculate TP, FP, FN, and TN totals
    total_tp=$(wc -l < kmcp/f1_score/"TP${i}_R1.txt")
    total_fp=$(wc -l < kmcp/f1_score/"FP${i}_R1.txt")
    total_fn=$(wc -l < kmcp/f1_score/"FN${i}_R1.txt")
    total_tn=$(wc -l < kmcp/f1_score/"TN${i}_R1.txt")
    # Output the results to result.txt
    echo "TP${i}  $total_tp" > kmcp/f1_score/result${i}.txt
    echo "FP${i}  $total_fp" >> kmcp/f1_score/result${i}.txt
    echo "FN${i}  $total_fn" >> kmcp/f1_score/result${i}.txt
    echo "TN${i}  $total_tn" >> kmcp/f1_score/result${i}.txt
    # Calculate precision, recall, accuracy, and F1 score
    precision=$(echo "scale=4; $total_tp / ($total_tp + $total_fp)" | bc)
    recall=$(echo "scale=4; $total_tp / ($total_tp + $total_fn)" | bc)
    accuracy=$(echo "scale=4; ($total_tp + $total_tn) / ($total_tp + $total_fp + $total_fn + $total_tn)" | bc)
    f1_score=$(echo "scale=4; 2 * $total_tp / (2 * $total_tp + $total_fp + $total_fn)" | bc)
    # Output the results to result.txt
    echo "Precision${i}  $precision" >> kmcp/f1_score/result${i}.txt
    echo "Recall${i}  $recall" >> kmcp/f1_score/result${i}.txt
    echo "Accuracy${i}  $accuracy" >> kmcp/f1_score/result${i}.txt
    echo "F1 Score${i}  $f1_score" >> kmcp/f1_score/result${i}.txt
done
```
## 4.6 运行KrakenUniq 1.0.4
### 1)构建KrakenUniq的数据库
```
conda activate krakenuniq
#Note that KrakenUniq natively supports Kraken 1 databases (however not Kraken 2)
krakenuniq-download --db rice taxonomy
#you can also download it by youself and upload it to rice/taxonomy
#tar -C rice/taxonomy -zxvf rice/taxonomy/taxdump.tar.gz
#prepare basic files for constrcuting the database
cp fasta rice/library
#prepare rice.map file, and put it into rice/library
time memusg krakenuniq-build --db rice --kmer-len 31  --taxids-for-genomes --taxids-for-sequences --build --jellyfish-bin /public/home/liuyongxin/miniconda3/envs/krakenuniq/bin/jellyfish
#here choose the jellyfish 1 which you installed before
##Fasta file not downloaded from NCBI may need their taxonomy information assigned explicitly.
##This can be done using the string kraken:taxid|XXX in the sequence ID.
##Such as >GWHBFPX00000001|kraken:taxid|39947  Adapter sequence[methods see kraken2]
```
### 2)运行KrakenUniq
```
conda activate krakenuniq
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
file=/public/home/liuyongxin/gyy/230903camisim/rice06/database/krakenuniq/rice
mkdir -p ${path}/krakenuniq
mkdir -p ${path}/krakenuniq_qc
for i in {0..4};do
d=anonymous_read${i}
time memusg krakenuniq --db ${file}  \
	--threads 8 --report-file ${path}/krakenuniq/anonymous${i}.report      \
	--output ${path}/krakenuniq/anonymous${i}.output --paired \
	${path}/seq/${d}_1.fq ${path}/seq/${d}_2.fq  --hll-precision 12 
	krakenuniq-extract-reads -p 0 ${path}/krakenuniq/anonymous${i}.output \
	${path}/seq/${d}_%.fq > krakenuniq_qc/4530remain${i}.fq
done


path=/data2/gaoyunyun/project/camisim/230903camisim/human06/seq15-2
mkdir -p ${path}/krakenuniq
mkdir -p ${path}/krakenuniq_qc
file=/data2/gaoyunyun/project/camisim/230903camisim/rice06/database/human/
for i in {0..4};do
d=anonymous_read${i}
time memusg krakenuniq --db ${file}  \
--threads 8 --report-file ${path}/krakenuniq/anonymous${i}.report      \
--output ${path}/krakenuniq/anonymous${i}.output --paired \
${path}/seq/${d}_1.fq ${path}/seq/${d}_2.fq  --hll-precision 12
	krakenuniq-extract-reads -p 0 ${path}/krakenuniq/anonymous${i}.output \
	${path}/seq/${d}_%.fq > krakenuniq_qc/4530remain${i}.fq
done



```
### 3)准备需要计算的数据
```
mkdir -p krakenuniq_qc/f1_score
mkdir -p krakenuniq_qc/temp
mv krakenuniq_qc/*.fq krakenuniq_qc/temp
```
### 4)去宿主数据的准确性KrakenUniq
#### STEP1：extract the label of all fastq files
```
for i in {0..4};do
	grep -Eo '^@[0-9A-Za-z_-]+/1$' krakenuniq_qc/temp/4530remain${i}.fq | sed 's/^@//' > krakenuniq_qc/4530remain_read${i}_1.txt
	grep -Eo '^@[0-9A-Za-z_-]+/2$' krakenuniq_qc/temp/4530remain${i}.fq | sed 's/^@//' > krakenuniq_qc/4530remain_read${i}_2.txt
	awk 'NR==FNR{labels[$1]; next} !($1 in labels){print $1}' krakenuniq_qc/4530remain_read${i}_1.txt doc/extracted_reads${i}_1.tsv | cut -f 1 > krakenuniq_qc/4530read${i}_1.txt
	awk 'NR==FNR{labels[$1]; next} !($1 in labels){print $1}' krakenuniq_qc/4530remain_read${i}_2.txt doc/extracted_reads${i}_2.tsv | cut -f 1 > krakenuniq_qc/4530read${i}_2.txt
done

```
#### STEP2：Prepare different values
```
for i in {0..4};do
    #  filter and separate "hostdata" and "microbiomedata"
        awk '$3 == "39947" { print > "krakenuniq_qc/f1_score/hostdata'"$i"'_R1.txt" } $3 != "39947" { print > "krakenuniq_qc/f1_score/microbiomedata'"$i"'_R1.txt" }' doc/extracted_reads${i}_1.tsv
    awk '$3 == "39947" { print > "krakenuniq_qc/f1_score/hostdata'"$i"'_R2.txt" } $3 != "39947" { print > "krakenuniq_qc/f1_score/microbiomedata'"$i"'_R2.txt" }' doc/extracted_reads${i}_2.tsv 
        # Calculate TP in R1
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530read${i}_1.txt krakenuniq_qc/f1_score/hostdata${i}_R1.txt > krakenuniq_qc/f1_score/TP${i}_R1.txt
        # Calculate TP in R2
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530read${i}_2.txt krakenuniq_qc/f1_score/hostdata${i}_R2.txt > krakenuniq_qc/f1_score/TP${i}_R2.txt
        # Calculate FP in R1
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530read${i}_1.txt krakenuniq_qc/f1_score/microbiomedata${i}_R1.txt > krakenuniq_qc/f1_score/FP${i}_R1.txt
        # Calculate FP in R2
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530read${i}_2.txt krakenuniq_qc/f1_score/microbiomedata${i}_R2.txt > krakenuniq_qc/f1_score/FP${i}_R2.txt
        # Calculate FN in R1
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530remain_read${i}_1.txt krakenuniq_qc/f1_score/hostdata${i}_R1.txt  > krakenuniq_qc/f1_score/FN${i}_R1.txt
        # Calculate FN in R2
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530remain_read${i}_2.txt krakenuniq_qc/f1_score/hostdata${i}_R2.txt > krakenuniq_qc/f1_score/FN${i}_R2.txt          
        # Calculate TN in R1
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530remain_read${i}_1.txt krakenuniq_qc/f1_score/microbiomedata${i}_R1.txt > krakenuniq_qc/f1_score/TN${i}_R1.txt
        # Calculate TN in R2
        awk 'FNR==NR{a[$1]; next} $1 in a' krakenuniq_qc/4530remain_read${i}_2.txt krakenuniq_qc/f1_score/microbiomedata${i}_R2.txt > krakenuniq_qc/f1_score/TN${i}_R2.txt
done
```
#### STEP3：Calculate different values
```
for i in {0..4}; do
    # Calculate TP, FP, FN, and TN totals
    total_tp_r1=$(wc -l < krakenuniq_qc/f1_score/"TP${i}_R1.txt")
    total_tp_r2=$(wc -l < krakenuniq_qc/f1_score/"TP${i}_R2.txt")
    total_tp=$(($total_tp_r1 + $total_tp_r2))
    total_fp_r1=$(wc -l < krakenuniq_qc/f1_score/"FP${i}_R1.txt")
    total_fp_r2=$(wc -l < krakenuniq_qc/f1_score/"FP${i}_R2.txt")
    total_fp=$(($total_fp_r1 + $total_fp_r2))
    total_fn_r1=$(wc -l < krakenuniq_qc/f1_score/"FN${i}_R1.txt")
    total_fn_r2=$(wc -l < krakenuniq_qc/f1_score/"FN${i}_R2.txt")
    total_fn=$(($total_fn_r1 + $total_fn_r2))
    total_tn_r1=$(wc -l < krakenuniq_qc/f1_score/"TN${i}_R1.txt")
    total_tn_r2=$(wc -l < krakenuniq_qc/f1_score/"TN${i}_R2.txt")
    total_tn=$(($total_tn_r1 + $total_tn_r2))
    # Output the results to result.txt
    echo "TP${i}  $total_tp" > krakenuniq_qc/f1_score/result${i}.txt
    echo "FP${i}  $total_fp" >> krakenuniq_qc/f1_score/result${i}.txt
    echo "FN${i}  $total_fn" >> krakenuniq_qc/f1_score/result${i}.txt
    echo "TN${i}  $total_tn" >> krakenuniq_qc/f1_score/result${i}.txt
    # Calculate precision, recall, accuracy, and F1 score
    precision=$(echo "scale=4; $total_tp / ($total_tp + $total_fp)" | bc)
    recall=$(echo "scale=4; $total_tp / ($total_tp + $total_fn)" | bc)
    accuracy=$(echo "scale=4; ($total_tp + $total_tn) / ($total_tp + $total_fp + $total_fn + $total_tn)" | bc)
    f1_score=$(echo "scale=4; 2 * $total_tp / (2 * $total_tp + $total_fp + $total_fn)" | bc)
    # Output the results to result.txt
    echo "Precision${i}  $precision" >> krakenuniq_qc/f1_score/result${i}.txt
    echo "Recall${i}  $recall" >> krakenuniq_qc/f1_score/result${i}.txt
    echo "Accuracy${i}  $accuracy" >> krakenuniq_qc/f1_score/result${i}.txt
    echo "F1 Score${i}  $f1_score" >> krakenuniq_qc/f1_score/result${i}.txt
done
```
