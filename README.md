# A Machine Learning-Based Calibration Method for Enhanced Accuracy and Consistency in m6A Epitranscriptome Mapping
Source code for the main results in our paper

![alt text](./figure/Graphical abstract.png "Graphical_abstract")

## Table of Contents 
- [Background](#General-description)
- [Pepline](#General-description)
  - [1. Develop a gold standard benchmark dataset ](#1-Develop-a-gold-standard-benchmark-dataset-)
  - [2. Prediction of RNA modified sites](#1-prediction-of-rna-modified-sites)
- [Dependencies and versions](#Dependencies-and-versions)
- [Citation](#Citation) 
- [Contact](#Contact) 


## Background
N6-methyladenosine (m6A) is the most prevalent and functionally significant mRNA modification in eukaryotes. However, discrepancies in m6A maps between studies have prompted concerns regarding the reliability of their biological validity, primarily attributed to non-specific antibody enrichment during immunoprecipitation (IP), which leads to false positives. To address this challenge, we developed a novel machine learning-based computational method and a gold standard benchmark dataset using mRNA and in vitro transcribed (IVT) samples to calibrate transcriptome-wide m6A maps. By integrating genomic features, we identify and eliminate non-specific antibody enrichment-induced false positives in MeRIP-seq, generating a high-accuracy m6A epitranscriptome map. The model interpretation results revealed that false positives predominantly occur on short exons and mRNAs with similar sequence contexts. Furthermore, our calibration function can be applied to other antibody-dependent base resolution techniques (e.g., miCLIP and m6ACLIP) to improve their consistency with antibody-independent techniques. We recommend incorporating this calibration approach into peak calling processes to standardize putative m6A sites from various antibody-based mapping techniques. Our method provides a systematic solution to the lack of consistency and reproducibility in m6A maps, paving the way for more precise epitranscriptomic studies.

## 1. Develop a gold standard benchmark dataset

### 1.1. 
