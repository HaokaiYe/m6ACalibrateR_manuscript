# M6A-Cali: Machine Learning-Based Calibration of m6A Epitranscriptome Mapping that Corrects Antibody Non-specific Binding

<p align="center">
  <img src="./plots/Geographic%20abstract.png" width="700" alt="Graphical abstract">
</p>

## Table of Contents 
- [Background](#Background)
- [Workflow](#Workflow)
  - [1. Training Data used in m6ACali](#1-Training-Data-used-in-m6ACali)
  - [2. Comparing Performance of ML Models and Feature Sets](#2-Comparing-Performance-of-ML-models-and-feature-sets)
  - [3. Impact of Exon Length and mRNA Length on Identifying False Positive m6A](#3-Impact-of-Exon-Length-and-mRNA-Length-on-Identifying-False-Positive-m6A)
  - [4. m6ACali Accurately Identifies False Positive m6A Sites](#4-m6ACali-Accurately-Identifies-False-Positive-m6A-Sites)
  - [5. m6ACali Generalizes to Independent Datasets and New Techniques](#5-m6ACali-Generalizes-to-Independent-Datasets-and-New-Techniques)
  - [6. m6ACali Achieves Higher Performance under Rigorous Threshold](#6-m6ACali-Achieves-Higher-Performance-under-Rigorous-Threshold)
  - [7. Randomly Capturing High-coverage Consensus Sequences to Reconstruct False Positive m6A Landscapes](#7-Randomly-Capturing-High-coverage-Consensus-Sequences-to-Reconstruct-False-Positive-m6A-Landscapes)
- [Dependencies and versions](#Dependencies-and-versions)
- [Citation](#Citation) 
- [Contact](#Contact) 


## Background
- **N6-methyladenosine** (m6A) is the most prevalent and functionally significant mRNA modification in eukaryotes. 
- However, discrepancies in m6A maps between studies have prompted concerns regarding the reliability of their biological validity, primarily attributed to **non-specific antibody enrichment** during immunoprecipitation (IP), which leads to **false positives**. 
- To address this challenge, we developed a novel **machine learning-based computational method** and a **gold standard benchmark dataset** using mRNA and in vitro transcribed (IVT) samples to **calibrate** transcriptome-wide m6A maps. 
- By integrating **genomic features**, we identify and eliminate non-specific antibody enrichment-induced false positives in MeRIP-seq, generating a high-accuracy m6A epitranscriptome map. 
- The model interpretation results revealed that false positives predominantly occur on **short exons and mRNAs** with **similar sequence contexts**. 
- Furthermore, our calibration function can be applied to other **antibody-dependent base resolution techniques** (e.g., miCLIP and m6ACLIP) to improve their **consistency** with antibody-independent techniques. 
- We recommend incorporating this calibration approach into **peak calling** processes to standardize putative m6A sites from various **antibody-based mapping techniques**. 
- Our method provides a systematic solution to the lack of **consistency** and **reproducibility** in m6A maps, paving the way for more precise epitranscriptomic studies.

## Workflow 
*Using SYSY dataset as example*


### 1. Establishing the Training Data for m6ACali

1.1. Download the raw sequencing data from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151028)

- All SRR info are listed in `SRR_list.txt`.

```{bash}
fastq-dump --split-3 SRR14765584 -O ~/fastq
```

1.2. Eliminate adaptors and nucleotides of low quality using Trim Galore
```{bash}
trim_galore --stringency 3 --paired -o ~/trimmed ~/fastq/SRR14765584_1.fastq ~/fastq/SRR14765584_2.fastq
```

1.3. Match the processed reads to the reference genome UCSC hg38 using HISAT2

- UCSC hg38 genome index was downloaded from http://daehwankimlab.github.io/hisat2/download/

```{bash}
cd ~/index/hg38
wget https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
tar -zxvf hg38_genome.tar.gz

hisat2 -x ~/index/hg38/hg38/genome -1 ~/trimmed/SRR14765584_1_val_1.fq -2 ~/trimmed/SRR14765584_2_val_2.fq -S ~/sam/SRR14765584.sam
```

1.4. Covert the sam file to bam format to save space using samtools
```{bash}
samtools view -S ~/sam/SRR14765584.sam -b > ~/bam/SRR14765584.bam
```

1.5. Perform motif-based peak calling using exomePeak2

- Motif-based peak calling: replacing sliding windows with single base sites of the **DRACH consensus motif**, while keeping the rest of the peak calling procedures the same.

> The code implementation for peak calling can be found in `./code/peak_calling.R`. The resulting peaks are stored in `./rds/peaks_IVT.rds` and `./rds/peaks_mRNA.rds`.

---
<p align="center">
  <img src="./plots/benchmark_dataset.png" alt="Benchmark dataset topology">
</p>

As IVT RNA can be assured to be devoid of any modifications, it can serve as a negative control for the mRNA sample. Therefore, modification sites identified exclusively in the mRNA sample were considered **true positives**, while all sites identified in the IVT sample were deemed **false positives**.


### 2. Evaluating the Efficacy of Various ML Models and Feature Sets
2.1. Choose model

- We considered four popular machine learning models (Deep Learning, GLM, XGBoost, and Random Forest) and selected the one that best performed on the benchmark datasets.

<p align="center">
  <img src="./plots/compare_models.png" alt="Compare models">
</p>

> The code implementation for comparing models can be found in `./code/compare_models.R`. The resulting performances are stored in `./rds/compare_models.rds`.

2.2. Choose feature set

- We compared the performance of three feature sets to select the most suitable ones for training: sequence features, genomic features, and a combination of both. 

2.2.1 Sequence-derived features

- We encoded the sequence of the 50 base pairs upstream and downstream of the m6A site using one-hot encoding (e.g. A – [1, 0, 0, 0], U – [0, 1, 0, 0], C – [0, 0, 1, 0], G – [0, 0, 0, 1]).

2.2.2 Genome-derived features

- We systematically procured an array of genomic characteristics from exon-only as well as intron-incorporated versions of individual genomic regions, encompassing exons, introns, genes, transcripts, 5'UTR, 3'UTR, and coding sequences.
- For each region, the collated genomic traits encompass an overlapping index, the length of the region, the proximity to the 5'/3' ends of the region, and the relative positioning of annotations within the regions, denoted as 0 for the leftmost and 1 for the rightmost position. 
 
<p align="center">
  <img src="./plots/compare_feature_sets.png" alt="Compare feature sets">
</p>

> The code implementation for comparing feature sets can be found in `./code/compare_feature_sets.R`. The resulting performances are stored in `./rds/compare_feature_sets.rds`.


### 3. Influence of Exon Length and mRNA Length in False Positive m6A Identification
3.1 Feature selection

- We employed reverse feature selection to streamline our data and spotlight the most salient genomic features for m6A site calibration.

<p align="center">
  <img src="./plots/feature_selection.png" alt="Feature selection">
</p>

> The code implementation for feature selection can be found in `./code/feature_selection.R`. The resulting performances are stored in `./rds/feature_selection.rds`.

3.2. Feature maps of the top 2 predictors

- The top two features (exon length and mRNA length) consistently explained the most significant portion of the model performances.

<p align="center">
  <img src="./plots/top2features.png" alt="Top 2 features">
</p>

> The code implementation for visualizing top 2 features can be found in `./code/top2features.R`. The resulting performances are stored in `./rds/top2features.rds`.

### 4. m6ACali: An Accurate Identifier for False Positive m6A Sites
4.1. Build up the final models

- We selectively honed in on the top-performing genomic features that yielded the maximum AUC in our Random Forest models.

<p align="center">
  <img src="./plots/final_models.png" alt="Final models">
</p>

> The code implementation for cross validation can be found in `./code/final_models.R`. The resulting performances are stored in `./rds/final_models.rds`.

4.2. Thorough analysis across DRACH motifs

- We executed an exhaustive analysis across the entire spectrum of DRACH consensus motifs.

<p align="center">
  <img src="./plots/thorough_motifs.png" alt="Thorough motifs">
</p>

> The code implementation for thorough DRACH analysis can be found in `./code/color_code.R`.

4.3. Cross validation

- We conducted cross-validation on benchmark datasets to assess the generalizability of the classifiers.

<p align="center">
  <img src="./plots/cross_validation.png" alt="Cross validation">
</p>

> The code implementation for cross validation can be found in `./code/cross_validation.R`. The resulting performances are stored in `./rds/cross_validation.rds`.

### 5. m6ACali's Broad Applicability: Independent Datasets and Novel Techniques

*To evaluate the utility of m6ACali with novel antibody-based datasets, we deployed the model across 24 MeRIP-Seq and 25 single-base resolution samples.*

5.1. Consistency with antibody-independent data

- In response to the likelihood of non-specific antibody binding introducing false positives in antibody-dependent data, our goal was to assess the effectiveness of our calibration model in boosting the consistency of m6A site detection when using antibody-independent data.

<p align="center">
  <img src="./plots/consistency.png" alt="Consistency">
</p>

> The code implementation for consistency can be found in `./code/consistency.R`. The resulting performances are stored in `./rds/consistency.rds`.

5.2. Validation of predicted false positives 

- Most of the false positives predicted in these newly analyzed datasets were validated by our in vitro transcribed (IVT) benchmark datasets.

<p align="center">
  <img src="./plots/FP_validation.png" alt="FP validation">
</p>

> The code implementation for validation of predicted FP can be found in `./code/cali_befo_aft.R` ("Validation of predicted false positives" part).

### 6. m6ACali's Enhanced Performance under Stringent Thresholds

6.1. Distribution of high-confidence m6A sites

- After discarding false positives, we observed a significant enrichment of calibrated antibody-based m6A sites around stop codons — an enrichment that intensifies with stricter calibration parameters.

> The code implementation for distribution can be found in `./code/cali_befo_aft.R` ("Topology gradient" part).

<p align="center">
  <img src="./plots/topology_gradient.png" alt="Topology gradient">
</p>

6.2. Sensitivity of true m6A site identification

- We leveraged antibody-free data, devoid of non-specific antibody binding events, as a benchmark to evaluate our model's prowess in detecting genuine m6A sites.

<p align="center">
  <img src="./plots/sensitivity.png" alt="Sensitivity">
</p>

> The code implementation for sensitivity can be found in `./code/cali_befo_aft.R` ("Sensitivity  gradient" and "Dot plot" parts).

### 7. Reconstructing False Positive m6A Landscapes via Random Capture of High-Coverage Consensus Sequences

**7.1. Similar sequence content amongst true positive and false positive m6A sites**

7.1.1 LR models

- We constructed a logistic regression for nucleotides surrounding the m6A site, represented by one-hot encoding, to calculate the coefficient value for individual nucleotides at specific positions.

<p align="center">
  <img src="./plots/nucleotides.png" alt="Nucleotides">
</p>

7.1.2 Correlation test

- We then employed correlation tests `R = 0.7472, rho = 0.7947` to authenticated the correlation between logistic regression coefficients configured on high-confidence sites and those on false-positive sites.
 
<p align="center">
  <img src="./plots/cor_test.png" alt="Correlation test">
</p>

7.1.3 Chi-squared test

- We conducted a Chi-squared test `p = 2.155e-07` to evaluate the fit between observed frequencies and expected probabilities (1/4 for the same rank, 3/4 for different ranks), assuming a random association.

|         |  same rank  |  diff rank  |
|:-------:|-----------:|-----------:|
| Observed|       17.00|        7.00|
| Expected|        0.25|        0.75|

> The code implementation for analyzing the sequence content can be found in `./code/coefficients.R`.

**7.2. Reconstruction of false positive m6A landscapes**

7.2.1 Quantifying mapped reads

- We extracted DRACH motifs and proceeded to tally the reads overlapping each of these motifs.

7.2.2 Establishing high-coverage motifs

- For a more representative scenario, we focused on DRACH motifs with high-coverage (non-methylated motifs with read count in input samples exceeding the average count of true positive m6A sites).

> The code implementation for establishing high-coverage motifs can be found in `./code/cali_befo_aft.R` ("FP read count & FP prob (gradient)" part).

7.2.3 Metagene plot

- We hypothesize that the emergence of false positives could be attributed to m6A-specific antibodies haphazardly capturing consensus sequences, the likelihood of which is further influenced by read coverage.
 
<p align="center">
  <img src="./plots/FP_topology.png" alt="FP topology">
</p>

> The code implementation for metagene plot can be found in `./code/FP_topologys.R`.





## Dependencies and versions

- Command-line tools
  - sratoolkit 2.11.3: https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
  - cutadapt 2.8: https://files.pythonhosted.org/packages/94/e2/de61c38fbe04933045287fc27bfb56eebc388b16ee8e815ef6bf9f68b4ad/cutadapt-2.8.tar.gz
  - fastqc 0.11.9: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
  - trim_galore 0.6.6: https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz
  - hisat2 2.1.0: ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
  - samtools 1.10: https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2

- R packages
  
    | Package Name                        | Version     | Package Name                        | Version     |
    |-------------------------------------|-------------|-------------------------------------|-------------|
    | exomePeak2                          | 1.9.1       | Guitar                              | 2.8.0       |
    | predictiveFeatures                  | 0.99.94     | pROC                                | 1.18.0      |
    | SummarizedExperiment                | 1.22.0      | PRROC                               | 1.3.1       |
    | phastCons100way.UCSC.hg38           | 3.7.1       | ggplot2                             | 3.4.0       |
    | BSgenome.Hsapiens.UCSC.hg38         | 1.4.3       | ggrepel                             | 0.9.1       |
    | TxDb.Hsapiens.UCSC.hg38.knownGene   | 3.13.0      | ggtree                              | 3.0.4       |
    | stringr                             | 1.4.0       | aplot                               | 0.1.4       |
    | h2o                                 | 3.34.0.3    | reshape2                            | 1.4.4       |
    | randomForest                        | 4.7.1.1     | rtracklayer                         | 1.52.1      |


## Citation: 

## Contact
Please open an issue in the GitHub repo if you have any questions/doubts/suggestions about how to use this software. Thanks!
