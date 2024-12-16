[TOC]

# Benchmarking metagenomics tools for removing host contamination

This is the figure code for the research of Benchmarking metagenomics tools for removing host contamination. 

> Yunyun Gao, Hao Luo, Hujie Lyu, Haifei Yang, Salsabeel Yousuf, Shi Huang, Yong-Xin Liu, Benchmarking metagenomics tools for removing host contamination, **Journal**, PublicationDate. <doi>. <span></span>\<script async src="<https://badge.dimensions.ai/badge.js>" charset="utf-8">\</script>

## Abstract

The rapid evolution of metagenomic sequencing technology offers remarkable opportunities to explore the intricate roles of microbiome in host health and disease, as well as to uncover the unknown structure and functions of microbial communities. However, the swift accumulation of metagenomic data poses substantial challenges for data analysis. Contamination from host DNA can substantially compromise result accuracy, and increase additional computational resources by including non-target sequences. In this study, we assessed the impact of computational host-DNA decontamination on downstream analyses, highlighting its importance in producing accurate results efficiently. We also evaluated the performance of conventional tools like KneadData, Bowtie2, BWA, KMCP, Kraken2, and KrakenUniq, each offering unique advantages for different applications. Furthermore, we highlighted the importance of a host reference genome, noting that its absence negatively affected the decontamination performance across all tools. Our findings underscore the need for careful selection of decontamination tools and reference genomes to enhance the accuracy of metagenomic analyses.

The paper contains the following figures, and these scripts contain codes of all statistics and plotting figures.

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/GA.jpg "image")

**Graphical abstract | Benchmarking metagenomics tools for removing host contamination. **

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
**Fig. S5 | Comparative analysis of composition, computational efficiency and host contamination removal performance across six software. **

![image](https://github.com/YunyunGao374/HostPurge/blob/main/MSFigure/FigureS6.jpg "image")
**Fig. S6 | Assessing the performance of six tools at the genus level. **

## Figure code

All the figure code can be found in the corresponding folder

**Figure1**
Figure 1 contains the following components:
Figure1b-ResourceComparison: A scatter plot illustrating the time and memory consumption of different software.
Figure1d-MAGEvaluation: A heatmap showing the evaluation of the number of metagenome-assembled genomes (MAGs).
Figure1e-GeneFunctionalRelevanceAssessment: A correlation analysis presented to assess the relevance of gene ontology.

**Figure2**
Figure 2 contains the following components:
Figure2b-IndexDB: A bar chart comparing the time and memory usage of different software for indexing the host reference genome.
Figure2c-CPU: A box plot illustrating the memory usage of different software.
Figure2d-Time: A box plot showing the running time for different software.
Figure2e-CPU&Time: A heatmap displaying memory usage (top-right diagonal) and execution time (bottom-left diagonal) among different software, analyzed using the Kruskal-Wallis test.

**Figure3**
Figure 3 contains the following components:
Figure3a-Accuracy: A box plot illustrating the accuracy of six software tools.
Figure3b-Precision & Recall: A scatter plot depicting the precision-recall performance of the software.
Figure3c-F1: A scatter plot showing the impact of high host contamination rates and microbiome complexity on the F1-score across six software tools.
Figure3d-Composition: A species composition chart of the metagenomic dataset with a synthetic community, showing the results of host contamination removal by six software tools at a 90% host contamination level.
Figure3e-Radar Plot: A radar plot comparing the computational efficiency and host contamination removal performance across simulated 60 Gbps datasets with 90% host contamination.

**Figure4**
Figure 4 contains the following components:
Figure4b-heatmap: A heatmap displaying the Average Nucleotide Identity (ANI) analysis conducted using FastANI.
Figure4c-all: A box plot summarizing accuracy, precision, recall, and F1-score of six tools on simulated metagenomic data derived from the Oryza genus.
Figure4d-Accuracy: A bar chart illustrating the accuracy index of different software across varying proportions of host genome content.

**FigureS**
FigureS1-Metadata
FigureS1b: A stacked bar plot showing the proportion of publications mentioning the removal of host contamination.
FigureS1c: A frequency plot illustrating the usage frequency of existing software for removing host genomes.

FigureS3-Memory:
A bar plot displaying the memory usage of all software during the host removal process.
FigureS4-Time:
A bar plot showing the time consumption of all software during the host removal process.



## Pipelines
There are Three pipelines, 0HostDecontaminationImpactiononDownstreamAnalysis, 1HostDecontaminationSoftwareComparison and 2GenerateOfSimulatedData.

0. Host Decontamination Impact on Downstream Analysis:
This section addresses the impact of host contamination on downstream metagenomic analysis. It outlines how we process metagenomic data through reference-based methods, non-reference-based methods, and binning strategies.

1. Host Decontamination Impact on Downstream Analysis:
This section compares the effectiveness of existing host-decontamination software, including comparisons of their installation processes, operation methods, and decontamination performance.

2. Generation of Simulated Data:
This section details the process for generating the simulated data used in our study. It includes the methods for obtaining simulated data, information on where the data has been uploaded, and approaches to creating metagenomes of varying sizes.

## Contact

Prof. Yong-Xin Liu

Yong-Xin Liu's lab, Agricultural Genomics Institute at Shenzhen, Chinese Academy of Agricultural Sciences

No 97, Buxin Road, Dapeng District, Shenzhen 518120, Guangdong, China

E-mail: <liuyongxin@caas.cn>

Wechat: meta-genomics

Cite: Yunyun Gao, Hao Luo, Hujie Lyu, Haifei Yang, Salsabeel Yousuf, Shi Huang, Yong-Xin Liu, Benchmarking metagenomics tools for removing host contamination, **Journal**, PublicationDate. <doi>.
