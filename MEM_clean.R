library(lme4)
library(ggplot2)
library(stargazer)
library(blme)
library(sjPlot)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

setwd('~/Documents/Mixed_Effects_Model/AT_rerun/')

##################################
#       Test trial for MEM       #
##################################

# Read in file, build dataframe (columns should be covariates and cell populations, 
# row should be ID and subsequent info.)
me_memG = read.csv('grid_ME12_new.csv', sep = ';')
head(me_memG)
me_memG_df = data.frame(me_memG)

## Construct Mixed-Effect model
# Linear MEM
me_lmer.model = lmer(IgD.pos..Memory.B ~ sex + age + symptom_score + Group + active_treatment +
                        (1|id/Group), data=me_memG)
# Partial Bayesian MEM
me_blmer.model = blmer(IgD.pos..Memory.B ~ sex + age + symptom_score + Group + active_treatment +
                        (1|id/Group), data=me_memG)

# Rescale and center continuous parameters
numcols <- grep("^c\\.",names(me_memG_df))
dfs <- me_memG_df
dfs[,numcols] <- scale(dfs[,numcols])
me_lmer.model <- update(me_lmer.model,data=dfs)

# Evaluates whther a fitted mixed model is singular, if singular TRUE then opt for solution 
# One of these solutions being a partial bayesion model especially with complex models
isSingular(me_lmer.model, tol=1e-05)

# Print model output
summary(me_lmer.model)

# Create table of MEM model either with tab_model() or stargazer()
library(sjPlot)
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table') 
tab_model(me_lmer.model,me_mem4.model,me_mem20.model, file = 'tableME_grid_trial.html')

stargazer(me_lmer.model,
          type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "", out = 'tableME_grid_trial.html')

##################################
#   Input and Modelling - GRID   #
##################################

# Read in file, build dataframe (columns should be covariates and cell populations, 
# row should be ID and subsequent info.)
me_memG = read.csv('grid_ME12_new.csv', sep = ';')
head(me_memG)
me_memG_df = data.frame(me_memG)

# List of cell population names in order to loop through, omitting first 6 columns which 
# includes ID and covariates
G_list = colnames(me_memG_df[-(1:6)]) 
head(G_list)

# Mixed-effect modelling (Partial-bayesian)
varlist=G_list 
blups.models_G <- lapply(varlist, function(x) {
  mod2 = try(blmer(substitute(i ~ sex + age + symptom_score + Group + active_treatment + 
                   (1|id/Group), list(i = as.name(x))), 
                   data = me_memG_df, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)}
})

blups.models_nextG = as.list(blups.models_G)
# Remove NULL models that failed to converge
blups.models_nextG[sapply(blups.models_nextG, is.null)] <- NULL 

##################################
# Extract info from model - GRID #
##################################

library(predictmeans)

# Mixed-effect model expression
varlist1 <- lapply(blups.models_nextG, function(f) summary(f)$call[2])
varlist1

# Correlation Coefficient
estimate_AT <- lapply(blups.models_nextG, function(f) summary(f)$coefficients[6,1])
estimate_KOS <- lapply(blups.models_nextG, function(f) summary(f)$coefficients[5,1])
estimate_SS <- lapply(blups.models_nextG, function(f) summary(f)$coefficients[4,1])

# R^2 and adjusted
R2 <- lapply(blups.models_nextG, function(f) summary(f)$r.squared)
R2_adj <- lapply(blups.models_nextG, function(f) summary(f)$adj.r.squared)

# Median
med <- lapply(blups.models_nextG, function(f) summary(f)$residuals)
med_na <- lapply(med, function(f) na.exclude(f))
med_calc <- lapply(med_na, function(f) median(f))

## p-values for covariates (fixed-effect)
# SS
test_pSS <- lapply(blups.models_nextG, function(f) parameters::p_value_wald(f)[4,])
df.pvalueSS <- as.data.frame(test_pSS)
df.pvalueSS <-t(df.pvalueSS)
df.pvalueSS = df.pvalueSS[seq(0, nrow(df.pvalueSS), 2), ]
# KOS
test_pKOS <- lapply(blups.models_nextG, function(f) parameters::p_value_wald(f)[5,])
df.pvalueKOS <- as.data.frame(test_pKOS)
df.pvalueKOS <-t(df.pvalueKOS)
df.pvalueKOS = df.pvalueKOS[seq(0, nrow(df.pvalueKOS), 2), ]
# AT
test_pAT <- lapply(blups.models_nextG, function(f) parameters::p_value_wald(f)[6,])
df.pvalueAT <- as.data.frame(test_pAT)
df.pvalueAT <-t(df.pvalueAT)
df.pvalueAT = df.pvalueAT[seq(0, nrow(df.pvalueAT), 2), ]

# Prepare dataframe with extracted info. for downstream use
test_data = list(as.character(varlist1), med_calc, estimate_AT, estimate_KOS, estimate_SS, as.numeric(t(df.pvalueSS)), as.numeric(t(df.pvalueKOS)), as.numeric(t(df.pvalueAT)))
names(test_data) <- c('cells', 'Median', 'Estimate_AT', 'Estimate_KOS', 'Estimate_SS', 'pvalueSS', 'pvalueKOS', 'pvalueAT')

test_final <- as.data.frame(do.call(rbind, test_data))
t(test_final)
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table')
write.csv2(t(test_final), file='grid_MEMtable_newSS.csv')

##################################
#    Plasma Protein Expression   #
##################################

setwd('~/Documents/Mixed_Effects_Model/AT_rerun') 

me_mem1 = read.csv('Olink_ME12_new.csv', sep = ';')
head(me_mem1)
me_mem1_df = data.frame(me_mem1)

PP_list = colnames(me_mem1_df[-(1:6)])
head(PP_list)

varlist=PP_list 
blups.models_PP <- lapply(varlist, function(x) {
  mod2 = try(blmer(substitute(i ~ sex + age + symptom_score + Group + active_treatment + (1|id/Group), list(i = as.name(x))), 
                   data = me_mem1_df, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)}
})

blups.models_nextPP = as.list(blups.models_PP)
blups.models_nextPP[sapply(blups.models_nextPP, is.null)] <- NULL

varlist1 <- lapply(blups.models_nextPP, function(f) summary(f)$call[2])
varlist1

#ICC
library(nlme)
library(multilevel)
library(psychometric)
data("bh1996")
head(bh1996)
mpp = blmer(CDCP1 ~ sex + age + symptom_score + Group + active_treatment + (1|id/Group), data=me_mem1_df)
performance::icc(mpp, ppd=FALSE)
ICC1 <- psychometric::ICC1.lme(VCAN, active_treatment, data = mpp)

ICC1 <- lapply(varlist, function(x) {
  mod = try(psychometric::ICC1.lme(substitute(i, active_treatment, list(i = as.name(x))), 
                                   data = blups.models_nextPP[i]))
})

# odds ratio
estimate_AT <- lapply(blups.models_nextPP, function(f) summary(f)$coefficients[6,1])
estimate_KOS <- lapply(blups.models_nextPP, function(f) summary(f)$coefficients[5,1])
estimate_SS <- lapply(blups.models_nextPP, function(f) summary(f)$coefficients[4,1])

#Median
med <- lapply(blups.models_nextPP, function(f) summary(f)$residuals)
med_na <- lapply(med, function(f) na.exclude(f))
med_calc <- lapply(med_na, function(f) median(f))

# SS
test_pSS <- lapply(blups.models_nextPP, function(f) parameters::p_value_wald(f)[4,])
df.pvalueSS <- as.data.frame(test_pSS)
df.pvalueSS <-t(df.pvalueSS)
df.pvalueSS = df.pvalueSS[seq(0, nrow(df.pvalueSS), 2), ]
# KOS
test_pKOS <- lapply(blups.models_nextPP, function(f) parameters::p_value_wald(f)[5,])
df.pvalueKOS <- as.data.frame(test_pKOS)
df.pvalueKOS <-t(df.pvalueKOS)
df.pvalueKOS = df.pvalueKOS[seq(0, nrow(df.pvalueKOS), 2), ]
# AT
test_pAT <- lapply(blups.models_nextPP, function(f) parameters::p_value_wald(f)[6,])
df.pvalueAT <- as.data.frame(test_pAT)
df.pvalueAT <-t(df.pvalueAT)
df.pvalueAT = df.pvalueAT[seq(0, nrow(df.pvalueAT), 2), ]

test_data = list(as.character(varlist1), med_calc, estimate_AT, estimate_KOS, estimate_SS, as.numeric(t(df.pvalueSS)), as.numeric(t(df.pvalueKOS)), as.numeric(t(df.pvalueAT)))
names(test_data) <- c('proteins', 'Median', 'Estimate_AT', 'Estimate_KOS', 'Estimate_SS', 'pvalueSS', 'pvalueKOS', 'pvalueAT')

test_final <- as.data.frame(do.call(rbind, test_data))
t(test_final)
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table')
write.csv2(t(test_final), file='protein_MEMtable_newSS.csv')



##################################
#            mRNA                #
##################################
#options(max.print = 10000)
setwd('~/Documents/Mixed_Effects_Model/')

me_memRNA = read.csv('deseq_1.csv', sep = '\t')

#me_memRNA1 = as.numeric(unlist(me_memRNA))
head(me_memRNA[1:6])
me_memRNA_df1 = as.data.frame(me_memRNA)
head(me_memRNA_df1)

#me_memRNA_df = as.data.frame(t(me_memRNA))
mRNA_list = colnames(me_memRNA_df1[-(1:6)])
head(mRNA_list)

library(brms)
m1 = brm(ENSG00000000003 ~ sex + age + symptom_score + Group + active_treatment + (1|id/Group), data=me_memRNA_df1)
plot(m1, parse=c('symptom_score', 'active_treatment'))
launch_shinystan(m1)
library('variancePartition')
library('edgeR')
library('BiocParallel')
BiocManager::install("edgeR")
BiocManager::install("variancePartition")

str(summary(m1))#$coefficients[,3]
summary(m1)$coefficients
anova(m1)
anova(m1)$'Pr(>F)'[1]
VarCorr(m1)
m1$tTable[,4]
require(broom)
glance(m1)
names(summary(m1))
anova(m1, test='Wald')

#OPTION 1
list1 = list()
for (i in 1:length(mRNA_list)) {
  variable <- mRNA_list[i]
  m <- blmer(variable ~ sex + age + symptom_score + Group + active_treatment+ (1|id/Group), data=me_memRNA_df)
  list1[[i]] <- m
}

#OPTION 2
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/')

'''me_memRNA = read.csv('deseq_ME12.csv', sep = ';', header = FALSE, row.names=1)
head(me_memRNA[1:7])
#me_memRNA_df1 = as.data.frame(me_memRNA)
me_memRNA_df = as.data.frame(t(me_memRNA))
head(me_memRNA_df)
#as.numeric(me_memRNA_df)
mRNA_list = colnames(me_memRNA_df[-(1:6)])
head(mRNA_list)'''

me_memRNA = read.csv('deseq_2_newSS.csv', sep = ';')
head(me_memRNA[1:6])
me_memRNA_df1 = as.data.frame(me_memRNA)
head(me_memRNA_df1)
mRNA_list = colnames(me_memRNA_df1[-(1:6)])
head(mRNA_list)

varlist=mRNA_list #define what vars you want
blups.models_1 <- lapply(varlist, function(x) {
  mod2 = try(blmer(substitute(i ~ sex + age + symptom_score + Group + active_treatment + (1|id/Group), list(i = as.name(x))), 
                   data = me_memRNA_df1, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)}
})

blups.models_next1 = as.list(blups.models_1)
blups.models_next1[sapply(blups.models_next1, is.null)] <- NULL
blups.models_next1[0:4]

#extract p values
'''lapply(blups.models_next1, coefficients)#$coefficients[1,4]
library(nlme)
library(lmerTest)
lapply(blups.modelslmer[1:4], function(f) summary(f)$coefficients) #before null removal
lapply(blups.models_next1[1:4], function(f) anova(f))

library(remef)
coefNames <- term2coef(fit)
coefNames
coefNames <- lapply(blups.models_next1, function(x) term2coef(x, 'active_treatment'))
lapply(blups.models_next1, function(x) coef(summary(x))[coefNames, 4])

fstat <- lapply(blups.models_next1, function(f) summary(f)$fstatistic)
library(languageR)
pvals.fnc(m1)
anova(m1, ddf="Kenward-Roger")
m1 <-blmer(ENSG00000267524 ~ sex + age + symptom_score + Group + active_treatment + (1|id/Group), data=me_memRNA_df1)
x <- summary(m1)$call
y <- x[2]'''
#m2 <- lapply(blups.models_next1, function(f) parameters::p_value_wald(f)[4:6,])

#dataframe gene p value estimates for AT
#variables
varlist1 <- lapply(blups.models_next1, function(f) summary(f)$call[2])
head(varlist1)

# odds ratio
estimate_AT <- lapply(blups.models_next1, function(f) summary(f)$coefficients[6,1])
estimate_KOS <- lapply(blups.models_next1, function(f) summary(f)$coefficients[5,1])
estimate_SS <- lapply(blups.models_next1, function(f) summary(f)$coefficients[4,1])

#Median
med <- lapply(blups.models_next1, function(f) summary(f)$residuals)
med_na <- lapply(med, function(f) na.exclude(f))
med_calc <- lapply(med_na, function(f) median(f))

# SS
test_pSS <- lapply(blups.models_next1, function(f) parameters::p_value_wald(f)[4,])
df.pvalueSS <- as.data.frame(test_pSS)
df.pvalueSS <-t(df.pvalueSS)
df.pvalueSS = df.pvalueSS[seq(0, nrow(df.pvalueSS), 2), ]
# KOS
test_pKOS <- lapply(blups.models_next1, function(f) parameters::p_value_wald(f)[5,])
df.pvalueKOS <- as.data.frame(test_pKOS)
df.pvalueKOS <-t(df.pvalueKOS)
df.pvalueKOS = df.pvalueKOS[seq(0, nrow(df.pvalueKOS), 2), ]
# AT
test_pAT <- lapply(blups.models_next1, function(f) parameters::p_value_wald(f)[6,])
df.pvalueAT <- as.data.frame(test_pAT)
df.pvalueAT <-t(df.pvalueAT)
df.pvalueAT = df.pvalueAT[seq(0, nrow(df.pvalueAT), 2), ]

test_data = list(as.character(varlist1), med_calc, estimate_AT, estimate_KOS, estimate_SS, as.numeric(t(df.pvalueSS)), as.numeric(t(df.pvalueKOS)), as.numeric(t(df.pvalueAT)))
names(test_data) <- c('genes', 'Median', 'Estimate_AT', 'Estimate_KOS', 'Estimate_SS', 'pvalueSS', 'pvalueKOS', 'pvalueAT')

test_final <- as.data.frame(do.call(rbind, test_data))
#t(test_final)
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table')
write.csv2(t(test_final), file='gene2_MEMtable_newSS.csv')



#### INPUT FOR GSEA ####
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/GSEA_ATnew/')
df <- read.csv('gene_pvalue_newALL.csv', sep=';')
x <- df[df$pvalueAT <= 0.05, ]
write.csv2(x, file='gsea_input_AT.csv')


#ENSEMBL ID to gene name
library(biomaRt)
setwd('~/Documents/Mixed_Effects_Model/')
me_memRNA = read.csv('deseq_1.csv', sep = '\t')
me_memRNA_df1 = as.data.frame(me_memRNA)
mRNA_list = colnames(me_memRNA_df1[-(1:6)])
head(mRNA_list)

mart <- biomaRt::useDataset("hsapiens_gene_ensembl", biomaRt::useMart("ensembl"))
genes <- mRNA_list
df<-me_memRNA_df1[,-1]
head(df)
G_list <- biomaRt::getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
head(G_list)
write.table(as.data.frame(G_list),file="gene1_list.csv", quote=F,sep=";",row.names=F)

df_0 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("gene", "ensembl_gene_id")
colnames(df_0) <- x
df_0 <- merge(df,G_list,by.df_0="gene",by.df_0="ensembl_gene_id") #nope
head(df_0)

###### other #####

me_mem1$fit <- predict(me_mem1.model)
me_mem1.model

newdat <- expand.grid(Group=unique(me_mem1$Group),
                      symptom_score=c(min(me_mem1$symptom_score),
                                      max(me_mem1$symptom_score)))

library(ggplot2)
p <- ggplot(me_mem1, aes(x=Group, y=symptom_score, colour='GDNF')) +
  geom_point(size=3) +
  geom_line(aes(y=predict(me_mem5.model), group=id, size="id")) +
  #geom_line(data=newdat, aes(y=predict(me_mem5.model), size="Group")) +
  scale_size_manual(name="Predictions", values=c("id"=0.5, "Group"=3)) +
  theme_bw(base_size=22) 
print(p)

#Statistical significance
#1) construct null model first
me_mem1.null = lmer(IL8 ~ sex + age + symptom_score + 
                      (1|id), data=me_mem1, REML=FALSE)
#2) Re-do full model
me_mem1.model = lmer(IL8 ~ sex + age + symptom_score + (1|Group) + 
                       (1|id), data=me_mem1, REML=FALSE)
#3) Perform likelihood ratio test
anova(me_mem1.null, me_mem1.model)

#Looks at the coefficients of the model by subject and by item
coef(me_mem1.model)

#plot and colour
(colour_plot <- ggplot(me_mem1, aes(x = id, y = CCL19, colour = sex)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none")) 

plot(me_mem1.model)
qqnorm(resid(me_mem1.model))
qqline(resid(me_mem1.model))