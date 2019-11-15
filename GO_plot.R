setwd('~/Documents/GO_AT_plots/')

######################################
#    Extract genes from GO terms     #
######################################

# gsecc_GO file input
# originates from GSEA run done using clusterProfiler::gseGO()
# results used as readout to retrieve genes associated with specific GO terms

gsecc_GO <- read.csv('gsecc3_result_new80_FINAL.csv', sep = ';')
head(gsecc_GO)
gsecc_GO[5,]
x <- gsecc_GO[5,6]
x <- data.frame(do.call(rbind, strsplit(as.character(x), "/")))
x[1,]
x.vector <- c(t(x))

# Retrieve count information from Kallisto run (estimates) converted to counts using DESeq functions
# Used to specifically retrieve count info for all samples of particular genes associated with specific GO term

genes_counts <- read.csv('deseq_raw.csv', sep = ';')
head(genes_counts)
genes_counts[1,]
AT_row <- as.data.frame(genes_counts[1,])
data_genes <- genes_counts[genes_counts[,1] %in% x.vector,]
head(data_genes)

write.csv(data_genes, file='long-term synaptic potentiation-AT.csv') #without active treatment info

df_all <- rbind(data_genes, AT_row) #binds number of active treatments to each sample ID
write.csv(df_all, file='long-term synaptic potentiation+AT.csv')

######################################
#    HUGO conversion from Ensembl    #
######################################

# Read in file for HUGO identifier conversion
# As Ensembl identifiers used previously for DESeq2
d <- read.csv('long-term synaptic potentiation-AT.csv', sep = ';', row.names = 1)

# Query biomart
library(biomaRt)
ensembl_ids <- rownames(d)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
x <- getBM(
  mart = mart,
  values = ensembl_ids,
  filter = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "external_gene_name")
)

# Create a dictionary from ensembl id to gene name.
ens_to_gene <- as.character(x$external_gene_name)
names(ens_to_gene) <- as.character(x$ensembl_gene_id)

# Some ensembl ids do not map to any gene names.
gene_names <- ens_to_gene[rownames(d)]
gene_names[is.na(gene_names)] <- rownames(d)[is.na(gene_names)]

# Multiple ensembl ids map to a single gene name, so average those ensembl ids.
library(reshape2)
library(data.table)
mean_by <- function(dat, xs) {
  dat <- data.table(dat)
  dat$agg_var <- xs
  dat <- melt(dat, id.vars = "agg_var")
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = mean, na.rm = TRUE
  )
  rownames(dat) <- dat$agg_var
  dat[ , agg_var := NULL]
  dat
}

# Write the output
dd <- mean_by(d, gene_names)

# Dataframe containing solely Ensembl and HUGO as columns
write.table(gene_names, 'long-term synaptic potentiation_HUGO.csv') 

# Transpose and add AT covariate in place of ID
dd <- t(dd)
AT <- (t(AT_row))
dd$ID[match(AT$ID, dd$ID)] <- AT$active_treatment

# This is the moment to specify or filter genes (significant MEM gene match)

# Calculate median expression for each Gene by active treatment,
library(dplyr)
newdd <- dd %>% group_by(active_treatment) %>% summarise_each(funs(median))

# Option here for log2 ratio (this was what is used downstream)

# Dataframe should be rows (# of active treatments) and columns (HUGO identifiers of MEM sign.)

######################################
#   Plot log2(median expr.) vs. AT   #
######################################

library(grid)
library(directlabels)
library(ggplot2)
setwd('~/Documents/GO_AT_plots')
y <- read.csv('long-term synaptic potentiation_MEMsign_log2relB.csv', sep = ';')

d <- melt(y, id.vars = 'active_treatment')
p <- ggplot(d, aes(x=active_treatment, y=value, group=variable, colour=variable)) +  
  ggtitle('Long-term synaptic potentiation (GO:0060291)') +
  geom_path(aes(group = variable)) + labs(x= 'Active Treatments', y = 'log2(median expression)') +
  scale_colour_discrete(guide = 'none') +
  geom_dl(aes(label = variable), method =list('last.bumpup', cex = 0.6, hjust = -0.1)) +
  scale_x_continuous(limits=c(0,17)) +
  scale_y_continuous(limits=c(-0.26,0.18))
