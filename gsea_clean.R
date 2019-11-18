library(dplyr)
library(goseq)
library(data.table)
library(GO.db)
library(geneLenDataBase)
library(org.Hs.eg.db)
library(topGO)

setwd('~/Documents/Mixed_Effects_Model/AT_rerun/GSEA_ATnew/')

######################################
#               Input                #
######################################

# Dataframe built from MEM displays ensembl identifiers with associated p-values for AT 
data_X = read.csv('gsea_ATpvalues_new.csv', sep=';', header = FALSE)

# Sort list in decreasing order (second column includes p-values)
gene_list = data_X[,2]
names(gene_list)=as.character(data_X[,1])
gene_list = sort(gene_list, decreasing=TRUE)

######################################
#               GSEA                 #
######################################

# gseGO() run on hirarchical order list solely looking at small gene sizes for Biological
# Processes in order to not look at large parent gene sets (observed previously to setting
# parameters set)
gsecc <- clusterProfiler::gseGO(gene_list, OrgDb='org.Hs.eg.db', keyType = 'ENSEMBL', ont="BP", pAdjustMethod = 'BH', maxGSSize = 80, nPerm=1000)
gsecc3 <- clusterProfiler::simplify(gsecc) #removes redundant terms

# Results extracted for GO plots
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/GSEA_ATnew/')
c <- as.data.frame(gsecc3@result)
write.csv(c,file='gsecc3_result_new80.csv')

# Enrichment plot ordered by p-adjust values
Z <- enrichplot::dotplot(gsecc3, showCategory=40, orderBy='p.adjust', x='p.adjust', split=NULL)

# Coloring of plot with wesanderson palette
library('wesanderson')
mid <- 0.003 
Z + scale_color_gradient2(midpoint=mid, low="#3B9AB2", mid="#EBCC2A", high="#F21A00", space='Lab')
pal[3]
pal <- wes_palette("Zissou1", 5, type = "discrete")
Z + scale_fill_gradientn(colours = pal) + ggtitle('GSEA MEM AT')

