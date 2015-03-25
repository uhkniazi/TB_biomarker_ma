# berry_et_al_nature.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 23/04/2015
# Desc: analysis of berry et. al. dataset to select biomarkers of interest

source('tb_biomarker_ma_header.R')

## data loading
# load the data, clean and create factors
# data was loaded first time using getGEO('GSE19491', GSEMatrix = T, destdir = 'Data_external/')
# object was saved and now using that
oExp = f_LoadObject(file.choose())
oExp = oExp$GSE19491_series_matrix.txt.gz

# print samples
as.data.frame(table(oExp$source_name_ch1))

## data normalization
# normalize and log2 transform the data using lumi
oExp.lumi = lumiT(oExp, 'log2')
# remove any NA data
exprs(oExp.lumi) = na.omit(exprs(oExp.lumi))

## data quality checks 
# plot pca before and after normalization
m = exprs(oExp.lumi)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
fSamples = oExp$source_name_ch1
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
p.old = par(mfrow=c(2,2))
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2 GSE19491, not normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3 GSE19491, not normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3 GSE19491, not normalized')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components GSE19491 - berry, not normalized')
par(p.old)
# remove the outlier groups from the data
# these can be seen on the pc2 and pc3 plots
m = pr.out$x[,1:3]
m = data.frame(m, fSamples)
i = which(m$PC2 < -100 & m$PC3 > 50)
i = c(i, which(m$PC2 < 50 & m$PC3 < -100))
c = col
c[i] = 'black'
plot(pr.out$x[,c(2,3)], col=c, pch=19, xlab='Z2', ylab='Z3',
main='outlier groups')
# remove the outliers
oExp.lumi = oExp.lumi[,-i]


## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()
# pca check on quality
m = exprs(oExp)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
fSamples = as.character(oExp$source_name_ch1)
fSamples = as.factor(fSamples)
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
par(mfrow=c(2,2))
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2 GSE19491, normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3 GSE19491, normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3 GSE19491, normalized')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components GSE19491 - berry, normalized')
par(p.old)

## create groups for analysis
# plot pca on the class separation factor i.e. title
title = as.character(oExp$title)
long = grep('^\\w+?long', title, perl=T)
sep = grep('^\\w+?_sep_', title, perl=T)
still = grep('Still', x = title)
asle = grep ('ASLE', x=title)
staph = grep('Staph', x=title)
strep = grep('Strep', x=title)
psle = grep('PSLE', x=title)
hcons = grep('-H-', x=title)
ltb = grep('LTB_.+', x = title)
ptb = grep('PTB_.+', x = title)
con = grep('CON_.+', x = title)
fGroups = rep(NA, length=length(title))
fGroups[ltb] = 'LTB'
fGroups[ptb] = 'PTB'
fGroups[con] = 'CON'
fGroups[still] = 'STILL'
fGroups[asle] = 'ASLE'
fGroups[staph] = 'STAPH'
fGroups[strep] = 'STREP'
fGroups[psle] ='PSLE'
fGroups[hcons] = 'HCON'
fGroups[long] = 'LONG'
fGroups[sep] = 'SEP'
fGroups = factor(fGroups)
# create a second factor with only 2 levels
# keep ptb at 1 for downstream predictions
fGroups.2 = as.character(fGroups)
fGroups.2[fGroups.2 != 'PTB'] = 'Others'
fGroups.2 = factor(fGroups.2, levels = c('Others', 'PTB'))

# plot pca on these classes 
fSamples = fGroups
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
pch.p = 1:length(col.p)
pch = pch.p[as.numeric(fSamples)]
# plot the pca components
plot.new()
legend('center', legend = unique(fSamples), col=col.p[as.numeric(unique(fSamples))], pch=pch.p[as.numeric(unique(fSamples))],
       lwd=2)
par(mfrow=c(2,2))
plot(pr.out$x[,1:2], col=col, pch=pch, xlab='Z1', ylab='Z2', lwd=2,
     main='PCA comp 1 and 2 GSE19491, normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=pch, xlab='Z1', ylab='Z3', lwd=2,
     main='PCA comp 1 and 3 GSE19491, normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=pch, xlab='Z2', ylab='Z3', lwd=2,
     main='PCA comp 2 and 3 GSE19491, normalized')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=pch, xlab='Z1', ylab='Z2', zlab='Z3', lwd=2,
            main='Plot of first 3 components GSE19491 - berry, normalized')
par(p.old)

### analysis
## extract data
# count matrix
mDat = exprs(oExp)
# phenotypic data
dfPdata = pData(oExp)
dfPdata$fGroups = fGroups
dfPdata$fGroups.2 = fGroups.2

## data cleaning and formatting before statistical analysis
# remove any rows with NAN data
f = is.finite(rowSums(mDat))
table(f)
mDat = mDat[f,]
# scale across samples i.e. along columns
mDat = scale(mDat)
# select a subset of genes based on coefficient of variation.
mDat = t(mDat)
# calculate the coef of var for each gene
cv = apply(mDat, 2, function(x) sd(x)/abs(mean(x)))
# check cv
summary(cv)
# cut data into groups based on quantiles of cv
cut.pts = quantile(cv, probs = 0:10/10)
groups = cut(cv, breaks = cut.pts, include.lowest = T)
groups = cut(cv, breaks = cut.pts, include.lowest = T, labels = 0:9)
iMean = apply(mDat, 2, mean)
iVar = apply(mDat, 2, var)
coplot((cv) ~ iMean | groups)
coplot(iVar ~ iMean | groups)
# choose genes with small cv
# f = cv < 0.2
# choosing groups from quantile 0 to 40
mDat = mDat[,groups %in% c(0, 1, 2, 3)]

# select a subset of genes that show differential expression
p.t = apply(mDat, 2, function(x) t.test(x ~ fGroups.2)$p.value)
p.w = apply(mDat, 2, function(x) wilcox.test(x ~ fGroups.2)$p.value)
p.t.adj = p.adjust(p.t, 'BH')
p.w.adj = p.adjust(p.w, 'BH')
t = names(p.t.adj[p.t.adj < 0.1])
w = names(p.w.adj[p.w.adj < 0.1])
n = unique(c(w, t))
f1 = n %in% t
f2 = n %in% w
f = f1 & f2
n2 = n[f]
mDat.sub = mDat[, colnames(mDat) %in% n2]
# keep this variable for later usage
mDat.sub.pval = mDat.sub
### model fitting and variable selection
## use random forest on all data for variable selection
dfData = as.data.frame(mDat.sub)
dfData$fGroups.2 = dfPdata$fGroups.2
set.seed(1)
rf.fit.1 = randomForest(fGroups.2 ~., data=dfData, importance = TRUE)

# save the results to save time for next time
dir.create('Objects', showWarnings = F)
save(rf.fit.1, file='Objects/berry.rf.fit.1.rds')

# variables of importance
varImpPlot(rf.fit.1)
df = importance(rf.fit.1)
dfImportant = df[order(df[,'MeanDecreaseAccuracy'], decreasing = T),]
rm(df)
# the prediction accuracy is perfect on the full training set
p = predict(rf.fit.1, dfData)
table(p, o=dfData$fGroups.2)
# select the top few genes
### TAG 1
x = dfImportant[,'MeanDecreaseAccuracy']
summary(x)
hist(x)
# select the top genes, change this cutoff if required
f = x >= 3
cvTopGenes = names(x)[f]
ivTopGenes.Error = x[f]

# subset the data based on these selected genes
mDat.sub = mDat[,colnames(mDat) %in% cvTopGenes]
dfData = as.data.frame(mDat.sub)
dfData$fGroups = fGroups.2

### Further variable classification check
### using CV and ROC
dfData.full = dfData
iBoot = 100

## as the 2 class proportions are not equal, fit random forest multiple times on random samples
## containing equal proportions of both classes and check variable importance measures
# fit random forests multiple times
# store results 
lVarImp = vector('list', iBoot)
for (i in 1:iBoot) {
  # get indices of the particular factors in data table
  ind.o = which(dfData.full$fGroups == 'Others')
  ind.p = which(dfData.full$fGroups == 'PTB')
  # take sample of equal in size from group Others to PTB
  ind.o.s = sample(ind.o, size = length(ind.p), replace = F)
  # sample of PTB groups, 
  ind.p.s = sample(ind.p, size=length(ind.p), replace=F)
  # join the sample indices together
  ind = sample(c(ind.o.s, ind.p.s), replace=F)
  # take sample from the full dataset
  dfData = dfData.full[ind,]
  # fit model
  fit.rf = randomForest(fGroups ~., data=dfData, importance = TRUE, ntree = 5000)
  # get variables importance
  df = importance(fit.rf)
  df = df[order(df[,'MeanDecreaseAccuracy'], decreasing = T),]
  # put in list
  lVarImp[[i]] = df
} # for

## put data for each boot of each gene together in a dataframe
df = NULL
for (i in 1:iBoot) df = rbind(df, lVarImp[[i]])
# convert rownames i.e. gene names to factors
f = as.factor(rownames(df))
# calculate mean and sd for each gene
tapply(df[,'MeanDecreaseAccuracy'], f, mean)
tapply(df[,'MeanDecreaseAccuracy'], f, sd)
# annotation of genes
dfAnnot = fData(oExp)[cvTopGenes,c('Symbol', 'ID')]
df = as.data.frame(df)
df$Symbol = as.character(dfAnnot[rownames(df),'Symbol'])
# boxplots 
par(mar=c(6,4,3,2)+0.1)
boxplot(df$MeanDecreaseAccuracy ~ df$Symbol, las=2)
rm(dfData)

## stop and decide
## if we like the variable count here, then its fine, or else go back to TAG 1 and select more variables
# remove the variable of choice after testing
i = which(cvTopGenes == 'ILMN_1676448')
cvTopGenes = cvTopGenes[-i]

### cross validation with ROC
#### CV with ROC
library(ROCR)
library(MASS)
mDat.sub = mDat[,colnames(mDat) %in% cvTopGenes]
dfData = as.data.frame(mDat.sub)
dfData$fGroups = fGroups.2
### Further variable classification check
### using CV and ROC
dfData.full = dfData
set.seed(1)
lPred = vector(mode = 'list', length = 20)
lLab = vector(mode = 'list', length=20)
iCv.error = NULL
for (oo in 1:20){
  t.lPred = NULL
  t.lLab = NULL
  # select a subset of equal numbers for others and PTB
  ind.o = which(dfData.full$fGroups == 'Others')
  ind.p = which(dfData.full$fGroups == 'PTB')
  ind.o.s = sample(ind.o, size = length(ind.p), replace = F)
  ind.p.s = sample(ind.p, size=length(ind.p), replace=F)
  ind = sample(c(ind.o.s, ind.p.s), replace=F)
  dfData = dfData.full[ind,]
  for (o in 1:1){
    # perform 10 fold cross validation
    k = 10
    folds = sample(1:k, nrow(dfData), replace = T, prob = rep(1/k, times=k))
    # choose the fold to fit and test the model
    for (i in 1:k){
      # check if selected fold leads to 0 for a class
      if ((length(unique(dfData$fGroups[folds != i])) < 2) || (length(unique(dfData$fGroups[folds == i])) < 2)) next
      # check if fold too small to fit model
      if (nrow(dfData[folds != i,]) < 3) next
      # fit model on data not in fold
      fit = lda(fGroups ~ ., data=dfData[folds != i,])
      # predict on data in fold
      pred = predict(fit, newdata = dfData[folds == i,])$posterior[,'PTB']
      name = paste('pred',oo, o, i,sep='.' )
      t.lPred[[name]] = pred
      name = paste('label',oo,o, i,sep='.' )
      t.lLab[[name]] = dfData$fGroups[folds == i] == 'PTB'
      pred = predict(fit, newdata = dfData[folds == i,])$class
      iCv.error = append(iCv.error, mean(pred != dfData$fGroups[folds == i]))
    }
  }
  t.lPred = unlist(t.lPred)
  t.lLab = unlist(t.lLab)
  lPred[[oo]] = t.lPred
  lLab[[oo]] = t.lLab
}

pred = prediction(lPred, lLab)
perf = performance(pred, 'tpr', 'fpr')
auc = performance(pred, 'auc')

plot(perf, main=paste('ROC Prediction of for', 'TB'),
     spread.estimate='stddev', avg='vertical', spread.scale=2)
auc = paste('auc=', signif(mean(as.numeric(auc@y.values)), digits = 3))
cv = paste('CV Error=', signif(mean(iCv.error), 3))
legend('bottomright', legend = c(auc, cv))
abline(0, 1, lty=2)

dfAnnot = fData(oExp)[cvTopGenes,c('Symbol', 'ID', 'Definition')]
# save the results
dir.create('Data_results', showWarnings = F)
write.csv(dfAnnot, file='Data_results/berry.var.importance.csv')

### graph analysis section

