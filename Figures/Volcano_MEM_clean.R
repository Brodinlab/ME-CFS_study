#Volcano Plots from MEM sign. features

setwd('~/Documents/Mixed_Effects_Model/AT_rerun/VolcanoPlots/VP_MEM/')

############## GENES #################

##################################
#             Input              #
##################################

# Read in generated table of MEM includes: ensembl identifiers, hugo identifiers, as well as
# correlation coefficient and p-values for fixed effects
# Focus downstream will be AT covariate
res_genes <- read.csv("VP_genes.csv", sep = ';', header=TRUE)
head(res_genes)

# Labels p-values less than 0.001 as significant for plotting
res_genes$Significant <- ifelse(res_genes$pvalueAT < 0.001, "FDR < 0.001", "Not Sig")

##################################
#              Plot              #
##################################
library(ggrepel)

# Grey dots not significant and red significant as labeled above, can also 
# label a subset of genes on plot 

ggplot(res_genes, aes(x = Corr_coeff_AT, y = -log10(pvalueAT))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(res_genes, pvalueAT < 0.001), 
    aes(label = hugo),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + scale_x_continuous(limits=c(-0.2,0.2))

############ PROTEINS ################

##################################
#             Input              #
##################################

# Read in generated table of MEM includes: protein names, as well as
# correlation coefficient and p-values for fixed effects
# Focus downstream will be AT covariate
res_pp <- read.csv("VP_proteins.csv", sep = ';', header=TRUE)
head(res_pp)

# Labels p-values less than 0.05 as significant for plotting
res_pp$Significant <- ifelse(res_pp$pvalueAT < 0.05, "FDR < 0.05", "Not Sig")

##################################
#              Plot              #
##################################

# Grey dots not significant and red significant as labeled above, can also 
# label a subset of genes on plot 

ggplot(res_pp, aes(x = Corr_coeff_AT, y = -log10(pvalueAT))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(res_pp, pvalueAT < 0.05), 
    aes(label = proteins),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + scale_x_continuous(limits=c(-0.2,0.2))

############## CELLS #################

##################################
#             Input              #
##################################

# Read in generated table of MEM includes: cell populations, as well as
# correlation coefficient and p-values for fixed effects
# Focus downstream will be AT covariate
res_cells <- read.csv("VP_cells.csv", sep = ';', header=TRUE)
head(res_cells)

# Labels p-values less than 0.05 as significant for plotting
res_cells$Significant <- ifelse(res_cells$pvalueAT < 0.05, "FDR < 0.05", "Not Sig")

##################################
#              Plot              #
##################################

# Grey dots not significant and red significant as labeled above, can also 
# label a subset of genes on plot 

ggplot(res_cells, aes(x = Corr_coeff_AT, y = -log10(pvalueAT))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(res_cells, pvalueAT < 0.05), 
    aes(label = cells),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + scale_x_continuous(limits=c(-0.2,0.2))
