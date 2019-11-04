library(MOFA)
library(MOFAdata)
library(MultiAssayExperiment)
library(reticulate)
library(mofapy)
library(rhdf5)

# Using a specific python binary
py_install("mofapy", envname = "r-reticulate", method="auto") 
use_python('/Users/lucie.rodriguez/.virtualenvs/r-reticulate/bin/python', required = TRUE) 
mofapy <- import('mofapy')
py_install("MOFAdata", envname = "r-reticulate", method="auto") 
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")

#Depended on R version, BiocManager may be needed to load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MOFA", version = "3.9")
BiocManager::install("MOFAdata", version = "3.9")
BiocManager::install("mofapy", version = "3.9")
BiocManager::install("ggplot", version = "3.9")

setwd('~/Documents/MOFA/') #set working directories
#Load Olink data
d_olink = read.csv('Olink_ME12.csv', sep=';', row.names = 1, header = TRUE)
d_olink = as.data.frame(t(d_olink))
head(d_olink)
#load Grid data
d_grid = read.csv('grid_ME_MOFA.csv', sep=';', row.names = 1, header = TRUE)
d_grid = as.data.frame(t(d_grid))
head(d_grid)
#Load DESeq2 data (VST)
d_deseq = read.csv('deseq_ME12_mat.csv', sep=';', row.names = 1)
d_deseq = as.data.frame(d_deseq)
head(d_deseq)
#Metadata
d_meta = read.csv('ME_metadata.csv', sep=';', row.names = 1)
d_meta = as.data.frame(d_meta)
head(d_meta)
d_meta

#MOFA object (MAE)
ME_data = list(d_olink, d_grid, d_deseq)
names(ME_data) <- c('Plasma Protein', 'Cell Population', 'mRNA')
mae_ME <- MultiAssayExperiment(experiments = ME_data, colData = d_meta)
MOFAobject <- createMOFAobject(mae_ME) 
MOFAobject
plotDataOverview(MOFAobject)

#Fit the MOFA model
DataOptions <- getDefaultDataOptions()
DataOptions

ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 10
ModelOptions

TrainOptions <- getDefaultTrainOptions()
TrainOptions$DropFactorThreshold <- 0.01
TrainOptions

MOFAobject <- prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)

#option to regress covariates HERE

#Run MOFA
MOFAobject <- runMOFA(MOFAobject)
#To save model use -> outfile = '~/Documents/MOFA/MOFAobject.hdf5'

#To use saved model
#MOFAobject <- h5read('~/Documents/MOFA/', 'MOFAobject.hdf5')

# Calculate the variance explained (R2) per factor in each view 
r2 <- calculateVarianceExplained(MOFAobject)
r2$R2Total
plotVarianceExplained(MOFAobject)

#Extract covariates
condition <- getCovariates(MOFAobject, 'Group')
condition

#Plot a heatmap of the loadings from multiple factors in a given view
plotWeightsHeatmap(
  MOFAobject, 
  view = "mRNA", 
  factors = 1:10,
  show_colnames = TRUE
)

#Plot the top loadings for a given factor and view
plotTopWeights(
  MOFAobject, 
  view = "mRNA", 
  factor = 2
)

#Plot all loadings for a given factor and view
plotWeights(
  MOFAobject, 
  view = "mRNA", 
  factor = 1
) 

#correlation plot between factors -> should be uncorrelated
plotFactorCor(MOFAobject, method = 'pearson')

#scatterplot between two factors, similar to PCA
plotFactorScatters(MOFAobject, factors = 1:5, color_by = 'ID', shape_by='Group') 

set.seed(1234)
clusters <- clusterSamples(
  MOFAobject, 
  k = 3,        # Number of clusters for the k-means function
  factors = 5   # factors to use for the clustering
)

plotFactorScatter(
  MOFAobject, 
  factors = 4:5, 
  color_by = clusters
)

plotFactorBeeswarm(
  MOFAobject,
  factors = 2,
  color_by = "Group"
)

set.seed(1234)
clusters <- clusterSamples(
  MOFAobject, 
  k = 2,        # Number of clusters for the k-means function
  factors = 2   # factors to use for the clustering
)

plotFactorScatter(
  MOFAobject, 
  factors = 1:2, 
  color_by = clusters
)

MOFAweights <- getWeights(
  MOFAobject, 
  views = "all", 
  factors = "all", 
  as.data.frame = TRUE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
)
MOFAweights
#can be used for further downstream analysis with Factors
write.csv(MOFAweights, file='MOFAweights.csv')

factor1 <- sort(getFactors(MOFAobject,"LF2")[,1])
order_samples <- names(factor1)
df <- data.frame(
  row.names = order_samples,
  factor = factor1
)

df1 <- getFactors(MOFAobject, 'all')
#can be used for further downstream analysis with Weights
write.csv(df1, file='AllFactors.csv')

plotDataHeatmap(
  object = MOFAobject, 
  view = "Cell Population", 
  factor = "LF1", 
  features = 20, 
  transpose = TRUE, 
  show_colnames=FALSE, show_rownames=TRUE, # pheatmap options
  cluster_cols = FALSE, annotation_col=df # pheatmap options
)

plotWeights(
  object = MOFAobject,
  view = "Plasma Protein", 
  factor = 3, 
  nfeatures = 0,
  abs = FALSE,
  scale = FALSE
)

#Correlation
library(ggplot2)
library(ggpubr)
ID <-getCovariates(MOFAobject,'ID') #Active.Treatment / Relative.SS
cdr <-getCovariates(MOFAobject,'log2.ratio') #Active.Treatment / Relative.SS

factor2 <- getFactors(MOFAobject,
                      factors=2)
LF2 <- data.frame(factor = as.numeric(factor2), cdr = cdr)
ggplot(LF2, aes_string(x = "factor", y = "cdr")) + 
  geom_point() + xlab("Factor 2") +
  ylab("log2.ratio") +
  stat_smooth(method="lm") +
  geom_line(aes(group=ID), lty = 2, colour = "purple") +
  theme_bw()
 

res.cor <- cor.test(foo$factor, foo$cdr, method = "pearson")
res.cor
