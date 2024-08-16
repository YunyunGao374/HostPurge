[TOC]

# Benchmarking metagenomics tools for removing host contamination

This is the figure code for the research of Benchmarking metagenomics tools for removing host contamination. 

> Yunyun Gao, Hao Luo, Hujie Lyu, Haifei Yang, Salsabeel Yousuf, Shi Huang, Yong-Xin Liu, Benchmarking metagenomics tools for removing host contamination, **Journal**, PublicationDate. <doi>. <span></span>\<script async src="<https://badge.dimensions.ai/badge.js>" charset="utf-8">\</script>

## Abstract

The rapid evolution of metagenomic sequencing technology offers remarkable opportunities to explore the intricate roles of microbiome in host health and disease, as well as to uncover the unknown structure and functions of microbial communities. However, the swift accumulation of metagenomic data poses substantial challenges for data analysis. Contamination from host DNA can substantially compromise result accuracy, and increase additional computational resources by including non-target sequences. In this study, we assessed the impact of computational host-DNA decontamination on downstream analyses, highlighting its importance in producing accurate results efficiently. We also evaluated the performance of conventional tools like KneadData, Bowtie2, BWA, KMCP, Kraken2, and KrakenUniq, each offering unique advantages for different applications. Furthermore, we highlighted the importance of a host reference genome, noting that its absence negatively affected the decontamination performance across all tools. Our findings underscore the need for careful selection of decontamination tools and reference genomes to enhance the accuracy of metagenomic analyses.

The paper contains the following figures, and these scripts contain codes of all statistics and plotting figures.

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/GA.jpg "image")

**Graphical abstract | Benchmarking metagenomics tools for removing host contamination.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure1.jpg "image")

**Fig. 1 | Host contamination consumed extra computing resources and affected the accuracy of the results in metagenomic analysis.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure2.jpg "image")
**Fig. 2 | Benchmarking calculates resources of six host removal software on simulated human and rice metagenomic data.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure3.jpg "image")
**Fig. 3 | Assessing the accuracy of host contamination removal across various software.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/Figure4.jpg "image")

**Fig. 4 | Impact of lacking a host reference genome on the performance of host decontamination tools.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS1.jpg "image")
**Fig. S1 | Literature search for existing studies addressing host contamination**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS2.jpg "image")
**Fig. S2 | Comparison of host removal performance on the accuracy of the results in metagenomic analysis.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS3.jpg "image")
**Fig. S3 |  Memory usage during the host removing process of BWA, Bowtie2, KneadData, KMCP, Kraken2, KrakenUniq in simulation rice and human metagenome. SinBac, Single Bacteria; SynCom, Synthetic Community.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS4.jpg "image")
**Fig. S4 | Time consumption during the host removing process of BWA, Bowtie2, KneadData, KMCP, Kraken2, KrakenUniq. SinBac, Single bacterium; SynCom, Synthetic community.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS5.jpg "image")
**Fig. S5 | Comparative analysis of composition, computational efficiency and host contamination removal performance across six software.**

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS6.jpg "image")
**Fig. S6 | Assessing the performance of six tools at the genus level.**

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

Cite: Yunyun Gao, Hao Luo, Hujie Lyu, Haifei Yang, Salsabeel Yousuf, Shi Huang, Yong-Xin Liu, Benchmarking metagenomics tools for removing host contamination, **Journal**, PublicationDate. <doi>.
