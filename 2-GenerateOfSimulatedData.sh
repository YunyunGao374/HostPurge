# 1. The generateion of simulated data for figure1 (9 datasets with Raw data, Removed data, and Microbiome data)
## (1)Three raw data can be found in NCBI with SAMN45841733, SAMN45841733, SAMN45841734, SAMN45841735 (S1, S2 and S3)
## (2)Three Removed data are generated from Raw data after Kneaddata
```
source activate kneaddata
kneaddata --version
#0.12.0
path=/public/home/liuyongxin/gyy/230903camisim/rice06/seq6-1
file=/public/home/liuyongxin/gyy/230903camisim/rice06/database/kneaddata/human
for i in {1..3};do
	d=S${i}
	time memusg kneaddata -i1 ${path}/seq/${d}_1.fq -i2 ${path}/seq/${d}_2.fq  \
	-o  ${path}/kneaddata_output \
	-db ${file}/Oryza_sativa \
	--bypass-trim  --bypass-trf\
	-t 8
done
```
## (3) Three Removed data are generated from Raw data based on the microbiome label
```
for i in {1..3};do
	d=S${i}
	grep -A 3 "microbiome" ${d}_1.fq | grep -v "^--$" > microbiome/${d}_1.fq
	grep -A 3 "microbiome" ${d}_2.fq | grep -v "^--$" > microbiome/${d}_2.fq
done
```

# 2.Generation of Simulated data for 180 dataset for figure2
# 2.1 CAMISM模拟数据
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

# 2.2 准备数据
```
#copy camisim simulate data (fq.gz fiesll) to new files mkdir 1_ModifiedData
```
## (1)修改数据格式
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
## (2)检查数据大小
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
# (3) Extracting 10G, 30G and 60G Datasets
## The analysis involves datasets from human and rice sources, with sizes of 10, 30, and 60 Gb and host contamination levels of 90%, 50%, and 10%. Each configuration includes synthetic bacterial (SynBac) and synthetic community (SynCom) samples with 5 replicates, resulting in a total of 180 datasets.

## Data Availability:
### All raw data are available in the NCBI database.
### For human data (host), the accession numbers are SAMN43249789 to SAMN43249793.
### For single bacterial samples from human (microbiome), the accession numbers are SAMN43249794 to SAMN43249798.
### For multiple bacterial samples from human (microbiome), the accession numbers are SAMN43249799 to SAMN43249803.

### For rice data (host), the accession numbers are SAMN43222004 to SAMN43222008.
### For single bacterial samples from rice (microbiome), the accession numbers are SAMN43222009 to SAMN43222013.
### For multiple bacterial samples from rice (microbiome), the accession numbers are SAMN43222014 to SAMN43222018.

#2_SelectedData, Selected 10G data from the total data, and host:microbiomel=9:1, and merged them
### A. Rice-SinBac(Single bacteria of rice)
```
seq=10
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
#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
#2_SelectedData,30G-host:microbiomel=9:1,seq=15-1;345592080:38399120
#2_SelectedData,30G-host:microbiomel=1:1,seq=15-2;191995600;191995600
#2_SelectedData,30G-host:microbiomel=1:9,seq=15-3;38399120:345592080
#2_SelectedData,60G-host:microbiomel=9:1,seq=30-1;691184160:76798240
#2_SelectedData,60G-host:microbiomel=1:1,seq=30-2;383991200:383991200
#2_SelectedData,60G-host:microbiomel=1:9,seq=30-3;76798240:691184160
```
### B. Human-SinBac(Single bacteria of human)
```
#2_SelectedData, Selected 10G data from the total data, and host:microbiomel=9:1, and merged them
seq=10-1
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
#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
#2_SelectedData,30G-host:microbiomel=9:1,seq=15-1;345592080:38399120
#2_SelectedData,30G-host:microbiomel=1:1,seq=15-2;191995600;191995600
#2_SelectedData,30G-host:microbiomel=1:9,seq=15-3;38399120:345592080
#2_SelectedData,60G-host:microbiomel=9:1,seq=30-1;691184160:76798240
#2_SelectedData,60G-host:microbiomel=1:1,seq=30-2;383991200:383991200
#2_SelectedData,60G-host:microbiomel=1:9,seq=30-3;76798240:691184160
```
### C. Rice-SynCom(Multiple bacterial samples from rice) 
```
seq=60-2
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

#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
#2_SelectedData,30G-host:microbiomel=9:1,seq=15-1;345592080:38399120
#2_SelectedData,30G-host:microbiomel=1:1,seq=15-2;191995600;191995600
#2_SelectedData,30G-host:microbiomel=1:9,seq=15-3;38399120:345592080
#2_SelectedData,60G-host:microbiomel=9:1,seq=30-1;691184160:76798240
#2_SelectedData,60G-host:microbiomel=1:1,seq=30-2;383991200:383991200
#2_SelectedData,60G-host:microbiomel=1:9,seq=30-3;76798240:691184160
```

### D. Human-SynCom(Multiple bacterial samples from human) 
```
#2_SelectedData, Selected 10G data from the total data, and host:microbiomel=9:1, and merged them
seq=60-2
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
    echo "Homo_sapiens $(grep -w -c "9606" 3extract${seq}/extracted_reads${i}_1.tsv)" > 3extract${seq}/meta${i}.txt
	echo "Microbiome $(grep -vw "9606" 3extract${seq}/extracted_reads${i}_1.tsv)" >> 3extract${seq}/meta${i}.txt
    done
#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
#2_SelectedData,30G-host:microbiomel=9:1,seq=15-1;345592080:38399120
#2_SelectedData,30G-host:microbiomel=1:1,seq=15-2;191995600;191995600
#2_SelectedData,30G-host:microbiomel=1:9,seq=15-3;38399120:345592080
#2_SelectedData,60G-host:microbiomel=9:1,seq=30-1;691184160:76798240
#2_SelectedData,60G-host:microbiomel=1:1,seq=30-2;383991200:383991200
#2_SelectedData,60G-host:microbiomel=1:9,seq=30-3;76798240:691184160
```


# 3.Generation of Simulated data for 60 dataset for figure 4

## Data Availability:
### All raw data are available in the NCBI database.

### For Refer (host), the accession numbers are SAMN43222004 to SAMN43222008.
### For Osj data (host), the accession numbers are SAMN45841746 to SAMN45841750.
### For Osi data (host), the accession numbers are SAMN45841741 to SAMN45841745.
### For Or data (host), the accession numbers are SAMN45841736 to SAMN45841740.
### For single bacterial samples from rice (microbiome), the accession numbers are SAMN43222009 to SAMN43222013.

### A. Refer-SinBac(Single bacteria of rice)
```
seq=10
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
#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
```

### B. OsjSinBac(Single bacteria of rice)
```
seq=10
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
#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
```

### C. OsiSinBac(Single bacteria of rice)
```
seq=10
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
#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
```


### D. OrSinBac(Single bacteria of rice)
```
seq=10
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
#2_SelectedData,10G-host:microbiomel=9:1,seq=06-1;143996688:15999632
#2_SelectedData,10G-host:microbiomel=1:1,seq=06-2;79998176:79998176
#2_SelectedData,10G-host:microbiomel=1:9,seq=06-3;15999632:143996688
```