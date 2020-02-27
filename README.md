# A workflow for identification and analysis of driver genes: a case study in breast cancer
---
## User Manual
---
#### 1. Introduction
This repository contains source code, original/preprocessed datasets and publication-quality figures for the paper "A workflow for identification and analysis of driver genes: a case study in breast cancer". This workflow aims to systematically integrate state-of-the-art computational tools to identify, characterize the driver genes and then predict their specific subgroups, which establish the basis for the understanding and further validations of driver genes in the future. We divide the work into two parts: (I) IDENTIFICATION and (II) ANALYSIS, in which the specific tasks of each part as follows:

(I) IDENTIFICATION
 - -Omics data is as input to identify driver genes by driver gene identification tools

(II) ANALYSIS
 - Depending on the cancer of interest and the purpose of research, several aspects of driver genes can be selected for analysis

The workflow is also be applied to the identification and analysis of driver genes for breast cancer patients. All statistical analyses were performed using R statistical software.

#### All data can be downloaded [here](https://drive.google.com/drive/folders/1v-W_ILNXbRHSvLlS0XqFGEQqkQ-q9kZm?usp=sharing). The data and source code should in the same folder!

#### 2. Pipeline

![Figure1](https://imgur.com/ujBkHCl.png)
Fig.1: Workflow for identification and analysis of driver genes. The workflow including two main steps is identification & analysis, in which, the former, -omics data is as input to identify driver genes by driver gene identification tools, and the latter, depending on the cancer of interest and the purpose of research, several aspects of driver genes can be selected for analysis. All tools used for each step can be selected depending on the purpose of the analyses and the input data. Abbreviation (for all sub-figures): CNVs, lymph and npi denote Copy number variants, the number of positive lymph nodes and the Nottingham prognostic index, respectively.

Here, we use the METABRIC breast cancer (BRCA) dataset as an example to demonstrate the use of our workflow to (1) identify and (2) analyze the driver genes.

![Figure2](https://imgur.com/tpJceRp.png)
Fig.2: Application of the workflow to the identification and analysis of driver genes for breast cancer patients.
#### 3. Implementation
**Step 1:** Identification of driver genes using the two tools [OncodriveCLUSTL](http://bbglab.irbbarcelona.org/oncodriveclustl/analysis) and [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/analysis#). </br>
**Step 2:** Enrichment analysis using the tool [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) </br>
**Step 3:** Individual gene-clincial feature association analysis
- Pre-processing gene expression data & clinical data: *"1.Clinical-preprocess.R"*
- Finding driver genes significantly associated with survival rate: *"2.Clinical-SA.R"*
- Identifying which driver genes are significantly correlated with the three remaining clinical features (i.e., number of lymph nodes, Nottingham prognostic index and pathologic stage): *"3.Clinical-corr.R"*

**Step 4:** Co-expressed module-clinical feature association analysis
- Constructing a co-expression gene network and analyzing the association between co-expressed module-clinical feature: *"1.WGCA-STEPbySTEP.R"*
- Investigating more the turquoise module to identify genes that have a high significance for the Nottingham prognostic index and then identify hub genes in the turquoise module: *"2.intramodular_analysis.R"*
- Visualizing the weighted gene co-expression network in terms of a heatmap: *"3.NetworkVisualization.R"*

**Step 5:** Stratification
- Pre-processing clinical data & CNVs data: *"1.CNA-preprocess.R"*
- Patient stratification by the identified driver genes: *"2.CNA-hclut.R"*
- Comparison of survival rate between groups of patients: *"3.CNA-SA.R"*
- Observing the differences between the two groups with the other clinical features: *"4.CNA-test.R"*

#### 4. Contact
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) or [Duc-Hau Le](https://github.com/hauldhut) (hauldhut@gmail.com) for any questions about the paper, datasets, code and results.


