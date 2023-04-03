# A Machine Learning-Based Calibration Method for Enhanced Accuracy and Consistency in m6A Epitranscriptome Mapping

![alt text](./figure/Graphical%20abstract.png "Graphical abstract")

## Table of Contents 
- [Background](#Background)
- [Pepline](#Pepline)
  - [1. Develop a gold standard benchmark dataset](#1-Develop-a-gold-standard-benchmark-dataset)
  - [2. Prediction of RNA modified sites](#1-prediction-of-rna-modified-sites)
- [Dependencies and versions](#Dependencies-and-versions)
- [Citation](#Citation) 
- [Contact](#Contact) 


## Background
N6-methyladenosine (m6A) is the most prevalent and functionally significant mRNA modification in eukaryotes. However, discrepancies in m6A maps between studies have prompted concerns regarding the reliability of their biological validity, primarily attributed to non-specific antibody enrichment during immunoprecipitation (IP), which leads to false positives. To address this challenge, we developed a novel machine learning-based computational method and a gold standard benchmark dataset using mRNA and in vitro transcribed (IVT) samples to calibrate transcriptome-wide m6A maps. By integrating genomic features, we identify and eliminate non-specific antibody enrichment-induced false positives in MeRIP-seq, generating a high-accuracy m6A epitranscriptome map. The model interpretation results revealed that false positives predominantly occur on short exons and mRNAs with similar sequence contexts. Furthermore, our calibration function can be applied to other antibody-dependent base resolution techniques (e.g., miCLIP and m6ACLIP) to improve their consistency with antibody-independent techniques. We recommend incorporating this calibration approach into peak calling processes to standardize putative m6A sites from various antibody-based mapping techniques. Our method provides a systematic solution to the lack of consistency and reproducibility in m6A maps, paving the way for more precise epitranscriptomic studies.

## Pepline 
*Using Abcam dataset as example*

### 1. Develop a gold standard benchmark dataset

### 1.1. 


## Dependencies and versions

Software | Version | Link
--- | --- | ---
sratoolkit | 2.11.3 | https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
cutadapt | 2.8 | https://files.pythonhosted.org/packages/94/e2/de61c38fbe04933045287fc27bfb56eebc388b16ee8e815ef6bf9f68b4ad/cutadapt-2.8.tar.gz
fastqc | 0.11.9 | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
trim_galore | 0.6.6 | https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz
hisat2 | 2.1.0 | ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
samtools | 1.10 | https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2


