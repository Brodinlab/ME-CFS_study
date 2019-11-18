# ME/CFS INMEST Study
ME/CFS study from Swedish clinical cohort (INMEST) for single-level analyses and omics integration.

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Repo description](#repo-description)
* [Single-level analyses](#single-level-analyses)
* [Omics integration](#omics-integration)
* [Figures](#figures)
* [Setup](#setup)

## General info
This project used multiple omics data:
- Plasma protein expression (Olink panels - NPX)
- Cell abundance (CyTOF - Grid)
- mRNA sequencing (VST counts)

Single-level analyses and omics-integration was performed in order to characterize this heterogeneous ME/CFS cohort and the effects of the INMEST treatment throughout the clinical trial.
	
## Dependencies
Project is created with:
* RStudio version: 3.6.0
* Python version: 2.7
* Unix/Linux

## Repo description
- ```DESeq2/``` contains script to prep Kallisto as well as convert Kallisto estimates and run DESeq2
- ```MOFA/``` contains script to train MOFA model and some functions to uncover sources of variation 
- ```Figures/``` contains scripts used for analyses and figures displayed in paper 
- ```MEM_clean.R``` used for mixed-effect modelling  

## Single-level analyses
### Grid
- Grid is a supervised learning algorithm based on t-SNE implementation. It uses the manual classification of cells in CyTOF samples and then the automatic classification of cells in new samples through the use of machine learning techniques based on these manual classifications.
- To run Grid, install it using:
```
$ pip install cellgrid
```
### DESeq2 (Differential Gene Expression Analysis)
- The principles behind DESeq2 is described in [Love et al. (2014)](https://dx.doi.org/10.1186%2Fs13059-014-0550-8)
- Used Kallisto outputs (estimates) that were converted into read counts using ```tximport``` before running DESeq2 with ```deseq_run.R``` 
### GSEA
### Mixed-Effect modelling
- Partial-bayesian mixed-effect modelling to account for covariates 
- ```MEM_clean.R``` extracts confounding variables for downstream use in GSEA for example, after modelling

## Omics integration
### MOFA
- Multi-Omics Factor Analysis (MOFA) was used in this study in order to deconvolute the main sources of variation in the differents sets of data mentioned above. For more information, read their published Methods paper [Argelaguet et al. (2018)](https://www.embopress.org/doi/10.15252/msb.20178124). 
- MOFA is publicly accessible here: https://github.com/bioFAM/MOFA 

## Figures
### Clinical response
### Spearman Correlation 
- ```spearman_corrmatrix.R``` uses cell abundance dataframe built from Grid that is sub-setted by active treatments for correlation
- Option to re-order one matrix in accordance to the other for comparison purposes
### GO plots
- Post GSEA, ```GO_plot.R``` extracts set of genes associated with GO terms to plot
- Original counts used to calculate median expression and log2 ratio for active treatments (covariate)
- Solely, significant features of MEM used
### Volcano MEM 
- ```Volcano_MEM_clean.R``` reads in generated MEM table built in ```MEM_clean.R``` to build volcano plots of features with significant ones highlighted and a top sub-set of them labeled 
### Metabolomic data trend

## Setup
To run this project, install it locally using devtools:

```
$ install.packages('devtools')
$ library(devtools)
$ install_github('rodriluc/ME-CFS_study')
$ library(ME-CFS_study)
```
