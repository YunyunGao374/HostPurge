[TOC]

# Benchmarking metagenomics tools for purging host contamination

This is the figure code for the research of Benchmarking metagenomics tools for purging host contamination. 

> Yunyun Gao, Hao Luo, Yong-Xin Liu,et al, Benchmarking metagenomics tools for purging host contamination, **Journal**, PublicationDate. <doi>. <span></span>\<script async src="<https://badge.dimensions.ai/badge.js>" charset="utf-8">\</script>

## Abstract

The rapid evolution of metagenomic sequencing technology offers remarkable opportunities to explore the intricate roles of the microbiome in host health and disease, as well as for exploring the unknown structure and functions of microbial communities. However, the swift accumulation of metagenomic data poses substantial challenges for data analysis. Host DNA contamination can significantly affect downstream analysis, and analyzing non-target sequences consumes extra computational resources. In this study, we summarized introduce the commonly used of host decontamination software, and compared the impact of host decontamination on computational resource and community result. Then we evaluated the performance of several existing software tools, including KneadData, Bowtie2, BWA, KMCP, Kraken2 and KrakenUniq, for their efficacy in host genome purge. Finally, we provide a new tool, HostPurge, to enhance the efficiency, reliability, and comparability of large-scale data mining.

The paper contains the following figures, and these scripts contain codes of all statistics and plotting figures.

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure-GA-F.jpg "image")

**Graphical abstract | Benchmarking metagenomics tools for purging host contamination. **

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure1-F.jpg "image")

**Fig. 1 | Comparative analysis of downstream performance with and without host contamination removal.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure2-F.jpg "image")
**Fig. 2 | Schematic of datasets design and different software performance evaluation.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure3%20-F.jpg "image")
**Fig. 3 | Assessing the efficacy of host contamination removal across various software.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure4%20-F.jpg "image")

**Fig. 4 | The comparison analysis of KneadData, Kraken2 and HostPurge.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS1-Metadata.jpg "image")
**Fig. S1 | Literature search for existing studies addressing host contamination**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS2-MemoryUsage.jpg "image")
**Fig. S2 | Memory usage during the host purging process of BWA, Bowtie2, KneadData, KMCP, Kraken2, KrakenUniq. SiBa, Single Bacteria; SyCo, Synthetic Community.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS3-Time%20consumption.jpg "image")
**Fig. S3 | Time consumption during the host purging process of BWA, Bowtie2, KneadData, KMCP, Kraken2, KrakenUniq. SiBa, Single Bacteria; SyCo, Synthetic Community.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS4-radar-summary.jpg "image")
**Fig. S4 | Comparative analysis of computational efficiency and host contamination removal performance.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS5-recall%2C%20precision%2C%20f1.jpg "image")
**Fig. S5 | The comparison of KneadData, Kraken2 and HostPurge.**

## Figure code

All the figure code can be found in the corresponding folder

## Pipelines
There are three pipelines, 0HostDecontaminationImpactiononDownstreamAnalysis, 1HostDecontaminationSoftwareComparison and 2HostPurgePipeline.

## Contact

Prof. Yong-Xin Liu

Yong-Xin Liu's lab, Agricultural Genomics Institute at Shenzhen, Chinese Academy of Agricultural Sciences

No 97, Buxin Road, Dapeng District, Shenzhen 518120, Guangdong, China

E-mail: <liuyongxin@caas.cn>

Wechat: meta-genomics

Cite: Yunyun Gao, Hao Luo, Yong-Xin Liu,et al, Benchmarking metagenomics tools for purging host contamination, **Journal**, PublicationDate. <doi>. <span></span>\<script async src="<https://badge.dimensions.ai/badge.js>" charset="utf-8">\</script>