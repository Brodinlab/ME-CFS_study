# ME/CFS INMEST Study
ME/CFS study from two clinical cohorts (INMEST) for single-level analyses and omics integration.

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Single-level Analyses](#single-level-analyses)
* [Omics integration](#omics-integration)
* [Setup](#setup)

## General info
This project used multiple data sets:
- Olink panels (plasma protein expression - NPX)
- CyTOF (Grid)
- mRNA sequencing (VST counts)

Single-level analyses and omics-integration was performed in order to characterize this heterogeneous ME/CFS cohort and the effects of the INMEST treatment throughout the clinical trial.
	
## Technologies
Project is created with:
* RStudio version: 3.6.0
* Python version: 2.7

## Single-level Analyses
- Grid
  - Grid is a software for the manual classification of cells in CyTOF samples and then the automatic classification of cells in new samples through the use of machine learning techniques based on these manual classifications.
  -To run Grid, install it using:
```
$ pip install cellgrid
```
- DESeq2 (Differential Gene Expression Analysis)
  - https://dx.doi.org/10.1186%2Fs13059-014-0550-8 
- GSEA (Gene Set Enrichment Analysis)
- Mixed-effect modeling

## Omics integration
Multi-Omics Factor Analysis (MOFA) was used in this study in order to deconvolute the main sources of variation in the differents sets of data mentioned above. For more information, read their published Methods paper (https://www.embopress.org/doi/10.15252/msb.20178124). MOFA is publicly accessible here: https://github.com/bioFAM/MOFA 
	
## Setup
To run this project, install it locally using devtools:

```
$ install.packages('devtools')
$ library(devtools)
$ install_github('rodriluc/ME-CFS_study')
$ library(ME-CFS_study)
```
