[TOC]
# The pipeline of reference-based metagenomic analysis
## Raw data with the host contamination
### HUManN3
#### Preparing files for HUMAnN3
```
cd /public/home/liuyongxin/gyy/camisim/multiplehuman/seq6-1/seq
mkdir -p free
cd free
mkdir -p temp/humann3
mkdir -p temp/concat
for i in {0..2};do
    d=anonymous_read${i}
    time memusg cat ../${d}_?.fastq \
    > temp/concat/${d}.fastq
done
```
#### Calculating the species and function composition using HUMAnN3
```
conda activate humann3
for i in {0..2};do
d=anonymous_read${i}
time memusg humann --input temp/concat/${d}.fastq --output temp/humann3 --threads 8 --metaphlan-options '--bowtie2db /public/home/liuyongxin/db/metaphlan4 --index mpa_vOct22_CHOCOPhlAnSGB_202212 --offline'
done
```
### Kraken2
#### Annotation with Kraken2
```
conda activate kraken2
mkdir -p temp/kraken2
db=~/db
for i in {0..2};do
    d=anonymous_read${i}
time memusg kraken2 --db ${db}/kraken2/pluspfp230605/ \
      --paired ../${d}_?.fastq \
      --threads 8 --use-names --report-zero-counts \
      --report temp/kraken2/${d}.report \
      --output temp/kraken2/${d}.output
done
```
## Raw data without the host contamination 
### Using Kneaddata to remove host contamination
```
conda activate kneaddata
file=~/gyy/camisim/rice06/database/kneaddata/human/
for i in {0..2};do
	d=anonymous_read${i}
    time memusg kneaddata -i1 seq/${d}_1.fastq -i2 seq/${d}_2.fastq \
      -o  ./kneaddata_output  --output-prefix ${d} \
      --bypass-trim --bypass-trf --reorder \
      --bowtie2-options '--very-sensitive --dovetail' \
      -db ${file}/human \
      --remove-intermediate-output -v -t 8
done
cd ./kneaddata_output
rm -rf *contam*
rm -rf *unmatched*
rename 's/_paired//' *_paired_1.fastq
rename 's/_paired//' *_paired_2.fastq
mv *.fastq ../remove
```
### HUManN3
#### Preparing files for HUMAnN3
```
cd /public/home/liuyongxin/gyy/camisim/multiplehuman/seq6-1/seq
mkdir -p free
cd free
mkdir -p temp/humann3
mkdir -p temp/concat
for i in {0..2};do
    d=anonymous_read${i}
    time memusg cat ../${d}_?.fastq \
    > temp/concat/${d}.fastq
done
```
#### Calculating the species and function composition using HUMAnN3
```
conda activate humann3
for i in {0..2};do
d=anonymous_read${i}
time memusg humann --input temp/concat/${d}.fastq --output temp/humann3 --threads 8 --metaphlan-options '--bowtie2db /public/home/liuyongxin/db/metaphlan4 --index mpa_vOct22_CHOCOPhlAnSGB_202212 --offline'
done
```
### Kraken2
#### Annotation with Kraken2
```
conda activate kraken2
mkdir -p temp/kraken2
db=~/db
for i in {0..2};do
    d=anonymous_read${i}
time memusg kraken2 --db ${db}/kraken2/pluspfp230605/ \
      --paired ../${d}_?.fastq \
      --threads 8 --use-names --report-zero-counts \
      --report temp/kraken2/${d}.report \
      --output temp/kraken2/${d}.output
done
```
## Deal with the data
```
# copy all files in a directory
# merge the result of samples

```

# The pipeline of reference-free metagenomic analysis
## Raw data with the host contamination
### Assembly by megahit
```
mkdir -p based/temp/megahit
cd based/
conda activate megahit
for i in {0..2};do
    d=anonymous_read${i}
    time memusg megahit -t 8 \
        -1 ../${d}_1.fastq \
        -2 ../${d}_2.fastq \
        -o temp/megahit/${d} 
done

```
### Bining by metawarp
```
conda activate metawrap
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/binning/${d}
time memusg metawrap binning \
-o temp/binning/${d} \
-t 8 -a temp/megahit/${d}/final.contigs.fa \
      --metabat2 --maxbin2 \
      --concoct  ../${d}*.fastq
done
```
### Bin refinement
```
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/bin_refinement/
time memusg metawrap bin_refinement \
      -o temp/${d}/bin_refinement/ \
      -A temp/binning/${d}/metabat2_bins/ \
      -B temp/binning/${d}/maxbin2_bins/ \
      -c 50 -x 10 -t 8
done
```
### dRep remove Redundant species/strain 
 ```
conda activate drep
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/drep95
time memusg dRep dereplicate temp/${d}/drep95 \
  -g temp/${d}/bin_refinement/metawrap_50_10_bins/*.fa  \
  -sa 0.95 -nc 0.30 -comp 50 -con 10 -p 5
done
``` 
### Annotation of bacteria
```
conda activate gtdbtk2.3
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/gtdb_classify
time memusg gtdbtk classify_wf \
    --genome_dir temp/${d}/drep95/dereplicated_genomes \
    --out_dir temp/${d}/gtdb_classify \
    --extension fa --skip_ani_screen \
    --prefix tax \
    --cpus 6
done
```

### CheckM2
```
conda activate checkm2
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/checkm2
time memusg checkm2 predict --input temp/${d}/drep95/dereplicated_genomes/*   --output-directory temp/${d}/checkm2 --threads 8
done
```

## Raw data without the host contamination
### Getting the clean data from bowtie2
### Assembly by megahit
```
mkdir -p based/temp/megahit
cd based/
conda activate megahit
for i in {0..2};do
    d=anonymous_read${i}
    time memusg megahit -t 8 \
        -1 ../${d}_1.fastq \
        -2 ../${d}_2.fastq \
        -o temp/megahit/${d} 
done

```

### Bining by metawarp
```
conda activate metawrap
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/binning/${d}
time memusg metawrap binning \
-o temp/binning/${d} \
-t 8 -a temp/megahit/${d}/final.contigs.fa \
      --metabat2 --maxbin2 \
      --concoct  ../${d}*.fastq
done
```
### Bin refinement
```
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/bin_refinement/
time memusg metawrap bin_refinement \
      -o temp/${d}/bin_refinement/ \
      -A temp/binning/${d}/metabat2_bins/ \
      -B temp/binning/${d}/maxbin2_bins/ \
      -c 50 -x 10 -t 8
done
```
### dRep remove Redundant species/strain 
 ```
conda activate drep
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/drep95
time memusg dRep dereplicate temp/${d}/drep95 \
  -g temp/${d}/bin_refinement/metawrap_50_10_bins/*.fa  \
  -sa 0.95 -nc 0.30 -comp 50 -con 10 -p 5
done
``` 
### Annotation of bacteria
```
conda activate gtdbtk2.3
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/gtdb_classify
time memusg gtdbtk classify_wf \
    --genome_dir temp/${d}/drep95/dereplicated_genomes \
    --out_dir temp/${d}/gtdb_classify \
    --extension fa --skip_ani_screen \
    --prefix tax \
    --cpus 6
done
```

### CheckM2
```
conda activate checkm2
for i in {0..2};do
d=anonymous_read${i}
mkdir -p temp/${d}/checkm2
time memusg checkm2 predict --input temp/${d}/drep95/dereplicated_genomes/*   --output-directory temp/${d}/checkm2 --threads 8
done
```
