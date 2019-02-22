# Metabolomics-Data-Portal
Metabolomics Data Portal R shiny application for the visualization and analysis of untargeted metabolomics datasets.

# Introduction:
Metabolomics is a fast maturing field which has an intimate relationship to the phenotypes observed in the clinic and is easily actionable for prospective treatment regimens. 

<img src=papers/Kennedy-et-al_2018.png width="250" align="center">

Within the field of metabolomics is the distinction between clinical research metabolomics, which follows a case-control cohort design; and clinical testing metabolomics which compares a single patient to a reference population.

<img src=papers/Kennedy-et-al_2018_2.png width="500" align = "center">

Differences in data collection percolate to differences in analysis needs. Currently, for N-of-1 clinical testing metabolomics, state of the art analysis methods rely on pathway enrichment methods.

<img src=papers/Burrage-et-al_2019.png" align="center">

To quantify perturbations observed in pathway knowledgebases, popular set-based methods such as over-representation analysis (ORA) and metabolite-set enrichment analysis (MSEA) are employed. However, these methods have been criticized for their use of gene sampling in lieu of patient sampling to generate p-values, and for their use of competitive null hypotheses in lieu of self-contained null hypotheses, which have shown to be less powerful (due to their less restrictive nature) in comparison (Goeman & Bulhmann, 2007). 

While various tools currently exist for metabolomics data analysis and pathway analysis (e.g., Metabolomics Workbench, PhenoMeNal and Metaboanalyst, Metscape, Mummichog, MetaMapp, and MetDisease), there are many shortcomings to these existing web tools. Some of these platforms employ popular machine learning models to analyze metabolomics data: unsupervised dimensionality reduction methods to view outliers or batch effects, and clustering methods to look for differences between cases and controls. Existing tools are not tailored for single patient analysis (i.e., N-of-1), such as in clinical testing metabolomics, and are more helpful for case-control cohort design data collection methods.

# Enter, topological enrichment methods!
Topological enrichment methods (good review papers found in Braun & Shah, 2005 and Ihnatova, Popovici & Budinska, 2018) have shown to be more sensitive than set-based enrichment analysis methods.

<img src=papers/Ihnatova-Popovici-Budinska_2018.png align="center">

# Problem
Modern day topological enrichment methods are all narrowly implemented for the analysis/interpretation of differentially expressed *gene sets*, and do not extend their functionality to the analysis and interpretation of perturbed metabolite sets.

# Our Solution
We have examined several R package implementation of existing topological enrichment methods and modified them to be useful for the analysis of metabolite sets, and for an N-of-1 metabolomics data.



Features:
1. Datasets included from published papers including clinical subjects with metabolic diseases.
2. Pathway visualization software and statistical interpretation metrics.
3. New topology-based pathway enrichment analysis methods.
4. Private data upload portal to use above tools on private datasets.

## Installation
- Dependencies:
- with [remotes](https://cran.r-project.org/web/packages/remotes/index.html)
```{r}
remotes::install_github("NCBI-Hackathons/Metabolomics-Data-Portal")
```

## Usage

## Data formats
- Input data
  - tabular data with rows as metabolites and columns as samples; data should be transformed Z-scores
- [Example data](https://github.com/NCBI-Hackathons/Metabolomics-Data-Portal/tree/master/data)



## Example Shiny Site
- Configure the docker-compose.yml file to point to your apps
- To spin up the shiny server, use the docker compose file from the command line:
```bash
docker-compose up -d
```
- To shut down the shiny server, use:
```bash
docker-compose down
```
- The dockerfile can also be built on your own:
```bash
docker build .
```

## References
- J.J. Goeman, P. Buhlmann. Analyzing gene expression data in terms of gene sets: methodological issues. Bioinformatics, 2007, 23(8):980-987.
- Miller MJ, Kennedy AD, Eckhart AD, Burrage LC, Wulff JE, Miller LA, et al. Untargeted metabolomic analysis for the clinical screening of inborn errors of metabolism. J Inherit Metab Dis. 2015;38:1029-39.
- M.F. Wangler, L. Hubert, T.R. Donti, M.J. Ventura, M.J. Miller, N. Braverman, K. Gawron, M. Bose,
A.B. Moser, R.O. Jones, W.B. Rizzo, V.R. Sutton, Q. Sun, A.D. Kennedy & S.H. Elsea. A metabolomic map of Zellweger spectrum disorders reveals novel disease biomarkers. Genetics in Medicine, 2018, 00. 
-MetaboLync Pathway Visualizations software, version 1.1.2, Copyright 2014 Metabolon, Inc., Research Triangle Park, NC, USA
- L.C. Burrage, L. Ashmore, B.M. Stroup, Q. Sun, M.J. Miller, S.C.S Nagamani, W. Craigen, F. Scaglia, V.R. Sutton, B. Graham, A.D. Kennedy, A. Milosavljevic, B.H. Lee,  S.H. Elsea. 2019. Untargeted Metabolomic Profiling Reveals Multiple Pathway Perturbations and New Clinical Biomarkers in Urea Cycle Disorders. 2019, Genetics in Medicine. 
- I. Ihnatova, V. Popovici & E. Budinska. A critical comparison of topology-based pathway analysis methods. PLOS One, 2018, 13(1): e0191154.
- R. Braun, S. Shah. Network Methods for Pathway Analysis of Genomic Data. arXiv, 2014, 1411.1993vl.
