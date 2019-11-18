library(gplots)

setwd('~/Documents/MOFA/')

##################################
#        Input and subset        #
##################################

# Read in file, dataframe has first column with # of active treatments and second column as ID
# the rest are cell populations (Grid data)
grid <- read.csv('grid_ME_AT-others.csv', sep='\t', header = TRUE)
head(grid)

# Subset the dataframe to create two different matrices, one for the samples at the time of the 
# last active treatment (AT 16) and one for baseline samples (AT 0)
grid_16 <- subset(grid, grepl("16", grid$active_treatment))
grid_16 <- grid_16[,-(1:2)]             

grid_0 <- subset(grid, grepl("0", grid$active_treatment))
grid_0 <- grid_0[,-(1:2)]  

##################################
#    Spearman Correlation        #
##################################

# Spearman correlation matrix was undertaken and rounded to 2 decimal places
cormat0 <- round(cor(grid_0, method = c("spearman")),2)
cormat16 <- round(cor(grid_16, method = c("spearman")),2)

# Re-order the correlation matrix by using correlation between variables as distance
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

#cormat0 <- reorder_cormat(cormat0)
cormat16 <- reorder_cormat(cormat16)

##################################
#      Prep for plotting         #
##################################

library(reshape2)
melted_cormat16 <- melt(cormat16)
head(melted_cormat16)

melted_cormat0 <- melt(cormat0)
head(melted_cormat0)

##################################
#              Plot              #
##################################

# Plots Corr. matrix of AT 16 re-ordered
library(ggplot2)
mc16 <- ggplot(data = melted_cormat16, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_fill_gradient2(low = "purple", high = "dark orange", mid = "white", 
                                     midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1))+
  coord_fixed() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 

# For comparison and visualization, AT 0 was reordered according to AT 16  corr. matrix

# Set vector of levels you want
mylevels <- melted_cormat16$Var1
# Re-order factors
melted_cormat0$Var1 <- factor(melted_cormat0$Var1,levels=unique(mylevels))
melted_cormat0$Var2 <- factor(melted_cormat0$Var2,levels=unique(mylevels))

# Plots Corr. matrix of AT 0 re-ordered to AT 16
mc0 <- ggplot(data = melted_cormat0, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_fill_gradient2(low = "purple", high = "dark orange", mid = "white", 
                                     midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1))+
  coord_fixed() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 


