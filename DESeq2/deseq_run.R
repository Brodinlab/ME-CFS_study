library(DESeq2)
library(tximportData)
library(readr)
library(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(GenomicFeatures)
library(sva)
library(makeTxDbFromGFF)
library(dplyr)
library(DESeqDataSet)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.9")
BiocManager::install('Homo_sapiens.GRCh38.90', version='3.9')
BiocManager::install('rtracklayer', version='3.9')
BiocManager::install('GenomicFeatures', version='3.9')
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("tximport")
BiocManager::install("tximportData")
BiocManager::install("sva")
BiocManager::install("makeTxDbFromGFF", version='3.9')
BiocManager::install("DESeqDataSet", version='3.9')


dir <- list.files('~/Documents/DESeq2/Cohort12')
dir

#use Kallisto output for sample txt
samples <- read.table(file.path('~/Documents/DESeq2/Cohort12', 'samples.txt'), header = TRUE)
samples

# Archive folder should include all .tsv files from Kallisto run
files <- file.path('~/Documents/DESeq2/Cohort12','Archive', samples$run) 
files
names(files) <- paste0('sample', 1:6)
all(file.exists(files))

gtffile <- file.path('~/Documents/DESeq2/Cohort12','Homo_sapiens.GRCh38.90.gtf')
txdb <- makeTxDbFromGFF(gtffile, format='gtf')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

txi <- tximport(files, type='kallisto', tx2gene=tx2gene, ignoreTxVersion = TRUE)
names(txi)

head(txi$counts)

txi.tx <- tximport(files, type = "kallisto", txOut = TRUE) #txi.kallisto

txi.sum <- summarizeToGene(txi.tx, tx2gene, ignoreTxVersion = TRUE)
all.equal(txi$counts, txi.sum$counts)

#DESeq2
samples$condition <- factor(rep(c("preKOS","postKOS")))
rownames(samples) <- samples$run
colnames(txi$counts) =samples$run
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition) 

#--------------------------------------------------------
#Filter out genes with low counts
rs <- rowSums(counts(ddsTxi))
minrs = 100
dds <- ddsTxi[ rs >= minrs, ]

#Normalize and VST option
dds <- estimateSizeFactors(dds)
dds <- varianceStabilizingTransformation(dds) 

#--------------------------------------------------------
#Remove batch effect - limma
mat1 <- assay(dds) 
library(limma)
mat2 <- removeBatchEffect(mat1, dds$cohort) #rename mat2 -> dds 

#Run DESeq2
dds$condition <- factor(dds$condition, levels = c("preKOS","postKOS"))
dds <- DESeq(dds)
dds
res=results(dds)
results(dds)

#Export results to CSV files
setwd('~/Documents/MOFA/') 

write.csv(as.data.frame(mat2), file = 'deseq_ME12_mat.csv') #use this in MOFA

write.csv(as.data.frame(res), file = 'deseq_ME12_results.csv')
write.csv(as.data.frame(txi$counts), file = 'deseq_ME12orig_txi_RESVISEDTEST.csv')

#Plot mean of normalized counts
resNorm <- lfcShrink(dds, coef=2, type="normal")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")

#Heatmap of the count matrix
library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition",'sample')]) #run
df
cdata <- colData(dds)
cdata
assay(ntd)
pheatmap(assay(ntd), cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df, show_colnames = FALSE)

#PCA
plotPCA(ntd, intgroup=("condition"))

#MA-plot: fold changes for all genes
topGene <- rownames(res)[which.min(res$padj)]
topGene
plotMA(res, ylim=c(-5,5))
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#Plot counts for individual genes
plotCounts(dds, gene=topGene, intgroup=c("condition"))

#FDR <0.5 
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
summary(res)
res
