# kaforou_et_al_nejm.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 25/04/2015
# Desc: analysis of kaforou dataset to select biomarkers of interest from new england journal of medicine

source('tb_biomarker_ma_header.R')

## data loading
# load the data, clean and create factors
# data was loaded first time using gse = getGEO('GSE39941', GSEMatrix = T, destdir = 'Data_external/')
# object was saved and now using that
oExp = f_LoadObject(file.choose())
oExp = oExp$GSE39941_series_matrix.txt.gz

# print samples
temp = pData(oExp)
# as.data.frame(table(temp$source_name_ch1))
# as.data.frame(table(temp$characteristics_ch1))
# as.data.frame(table(temp$characteristics_ch1.1))
# as.data.frame(table(temp$characteristics_ch1.2))
# as.data.frame(table(temp$characteristics_ch1.3))

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
# create factors
fSamples = as.character(oExp.lumi$characteristics_ch1)
i = grep('tuberculosis', fSamples, perl = T)
fSamples[i] = 'ATB'
fSamples[-i] = 'OD'
table(fSamples)
fSamples = factor(fSamples, levels = c('OD', 'ATB'))
oExp.lumi$fSamples = fSamples
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
p.old = par(mfrow=c(2,2))
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, not normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, not normalized')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components - not normalized')
par(p.old)
# remove the outlier groups from the data
# these can be seen on the pc2 and pc3 plots
m = pr.out$x[,1:3]
m = data.frame(m, fSamples)
i = which(m$PC1 < -180 & m$PC2 < 150)
i = unique(c(i, which(m$PC1 < -180 & m$PC3 > -100)))
c = col
c[i] = 'black'
# plot(pr.out$x[,c(2,3)], col=c, pch=19, xlab='Z2', ylab='Z3',
#      main='outlier groups')
f_Plot3DPCA(pr.out$x[,1:3], c, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components - not normalized')
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
fSamples = oExp$fSamples

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
par(mfrow=c(2,2))
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, normalized')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components - normalized')
par(p.old)

oExp$fGroups.2 = oExp$fSamples

##############################################
#### setup data for analysis
mDat = exprs(oExp)
dfSamples = pData(oExp)
mDat = t(mDat)
mDat.sub = mDat[,colnames(mDat) %in% cvTopGenes]

# create test and training datasets
test = sample(1:nrow(dfSamples), size = nrow(dfSamples) * 0.2, replace = F)
dfSamples.train = dfSamples[-test,]
dfSamples.test = dfSamples[test,]
mDat.sub.train = mDat.sub[-test,]
mDat.sub.test = mDat.sub[test,]

############################################
#### CV with ROC using LDA
# choose all data together for nested 10 fold cv
par(p.old)
dfData = as.data.frame(mDat.sub.train)
dfData$fGroups = dfSamples.train$fGroups.2
### Further variable classification check
### using CV and ROC
dfData.full = dfData
set.seed(1)
lPred = vector(mode = 'list', length = 50)
lLab = vector(mode = 'list', length=50)
iCv.error = NULL
for (oo in 1:50){
  t.lPred = NULL
  t.lLab = NULL
  # select a subset of equal numbers for others and ATB
  ind.o = which(dfData.full$fGroups == 'OD')
  ind.p = which(dfData.full$fGroups == 'ATB')
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
      pred = predict(fit, newdata = dfData[folds == i,])$posterior[,'ATB']
      name = paste('pred',oo, o, i,sep='.' )
      t.lPred[[name]] = pred
      name = paste('label',oo,o, i,sep='.' )
      t.lLab[[name]] = dfData$fGroups[folds == i] == 'ATB'
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

plot(perf, main=paste('ROC Prediction of for', 'ATB'),
     spread.estimate='stddev', avg='vertical', spread.scale=2)
auc.cv = paste('auc=', signif(mean(as.numeric(auc@y.values)), digits = 3))
cv.err = paste('CV Error=', signif(mean(iCv.error), 3))
#legend('bottomright', legend = c(auc, cv))
abline(0, 1, lty=2)

## fit model and roc without cross validation, just on test and training data
dfData.train = as.data.frame(mDat.sub.train)
dfData.train$fGroups = dfSamples.train$fGroups.2

dfData.test = as.data.frame(mDat.sub.test)
dfData.test = dfData.test[,colnames(dfData.test) %in% cvTopGenes]
dfData.test$fGroups = dfSamples.test$fGroups.2

fit = lda(fGroups ~ ., data=dfData.train)
# predict on data in fold
pred = predict(fit, newdata = dfData.test)$posterior[,'ATB']
ivPred = pred
ivLab = dfData.test$fGroups == 'ATB'
pred = predict(fit, newdata = dfData.test)$class
iCv.error = mean(pred != dfData.test$fGroups)

pred = prediction(ivPred, ivLab)
perf = performance(pred, 'tpr', 'fpr')
auc = performance(pred, 'auc')

plot(perf, add=T, lty=3, lwd=2, col=2)#main=paste('ROC Prediction of for', 'SA'))
auc.t = paste('t.auc=', signif(mean(as.numeric(auc@y.values)), digits = 3))
err.t = paste('t Error=', signif(mean(iCv.error), 3))
legend('bottomright', legend = c(auc.cv, cv.err, auc.t, err.t))
abline(0, 1, lty=2)

## plot these expression values for these genes
par(mfrow=c(1,2))
x = stack(dfData.train)
x$f = dfData.train$fGroups
boxplot(values ~ f+ind, data=x, las=2, par=par(mar=c(8, 4, 2, 2)+0.1), main='Training Data')

x = stack(dfData.test)
x$f = dfData.test$fGroups
boxplot(values ~ f+ind, data=x, las=2, par=par(mar=c(8, 4, 2, 2)+0.1), main='Test Data')

# get the gene annotation
dfAnnotation = fData(oExp)
dfAnnotation = dfAnnotation[,c('ID', 'ILMN_Gene', 'Entrez_Gene_ID', 'Symbol')]

df = dfAnnotation[cvTopGenes,]

################################################################################
#### Pathway analysis
mDat = exprs(oExp)
dfSamples = pData(oExp)
# scale across samples i.e. along columns
mDat = scale(mDat)
mDat = t(mDat)
# select a subset of genes based on coefficient of variation.
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
fGroups.2 = dfSamples$fGroups.2

# select a subset of genes that show differential expression
p.t = apply(mDat, 2, function(x) t.test(x ~ fGroups.2)$p.value)
p.w = apply(mDat, 2, function(x) wilcox.test(x ~ fGroups.2)$p.value)
p.t.adj = p.adjust(p.t, 'BH')
p.w.adj = p.adjust(p.w, 'BH')
t = names(p.t.adj[p.t.adj < 0.001])
w = names(p.w.adj[p.w.adj < 0.001])
n = unique(c(w, t))
f1 = n %in% t
f2 = n %in% w
f = f1 & f2
n2 = n[f]
mDat.sub = mDat[, colnames(mDat) %in% n2]
# keep this variable for later usage
mDat.sub.pval = mDat.sub

#### pathway clustering

### libraries
library(annotate)
library(org.Hs.eg.db)
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
library(igraph)
source('../CGraphClust/CGraphClust.R')
library(reactome.db)


# expression data
mExp = mDat.sub


### processing
# assign names to the genes from illumina ids to enterez ids
n = colnames(mExp)
i = match(n, dfAnnotation$ID)
# replacement names i.e. enterez ids
n2 = as.character(dfAnnotation$Entrez_Gene_ID[i])
colnames(mExp) = n2

# remove any duplicated enterez ids
f = duplicated(n2)
mExp = mExp[,!f]

### create clusters based on various criteria
# difference between hiv positive and negative atb
f = dfSamples$fSamples
i = which(f == 'ATB')
dfSamples.atb = dfSamples[i,]
mExp.atb = mExp[i,]

## get the reactome ids
dfGraph = AnnotationDbi::select(reactome.db, keys = colnames(mExp.atb), columns = 'REACTOMEID', keytype = 'ENTREZID')
dfGraph = na.omit(dfGraph)

# remove genes from expression matrix that do not have the 
# matching enterez id and reactome id
n = unique(dfGraph$ENTREZID)
mExp.atb = mExp.atb[,n]

# create correlation matrix
mCor = cor(mExp.atb)
par(p.old)
# choose cutoff based on histogram
hist(sample(mCor, 10000))

oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6)
ig = getFinalGraph(oGr)

# plot results
plot(ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot(getCommunity(oGr), ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot.mean.expressions(oGr, t(mExp.atb), fGroups = dfSamples.atb$characteristics_ch1.1)

### other criteria, remove hiv positive
f = as.character(dfSamples$characteristics_ch1.1)
i = grep(pattern = 'negative', f, perl = T)
dfSamples.atb = dfSamples[i,]
mExp.atb = mExp[i,]

# f = dfSamples.atb$fGroups.2
# i = which(f %in% c('ATB', 'LTBI', 'Ctrl'))
# dfSamples.atb = dfSamples.atb[i,]
# mExp.atb = mExp.atb[i,]

## get the reactome ids
dfGraph = AnnotationDbi::select(reactome.db, keys = colnames(mExp.atb), columns = 'REACTOMEID', keytype = 'ENTREZID')
dfGraph = na.omit(dfGraph)

# remove genes from expression matrix that do not have the 
# matching enterez id and reactome id
n = unique(dfGraph$ENTREZID)
mExp.atb = mExp.atb[,n]

# create correlation matrix
mCor = cor(mExp.atb)

# choose cutoff based on histogram
hist(sample(mCor, 10000))

oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.7)
ig = getFinalGraph(oGr)

f = as.character(dfSamples.atb$characteristics_ch1)
f2 = rep(NA, length.out=length(f))
i = grep('active', x = f, perl = T)
f2[i] = 'ATB'
i = grep('latent', x = f, perl = T)
f2[i] = 'LTBI'
i = grep('other', x = f, perl = T)
f2[i] = 'OD'


# plot results
plot(ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot(getCommunity(oGr), ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot.mean.expressions(oGr, t(mExp.atb), fGroups = as.character(f2), legend.pos = 'topleft', main='Correlation 0.7 and above')
m = t(mExp.atb)
colnames(m) = f2
m = m[,order(f2)]
png(filename = 'Temp/hm.png')
plot.heatmap(oGr, mCounts = m)
dev.off(dev.cur())

png(filename = 'Temp/hm2.png')
plot.heatmap.means(oGr, mCounts = m)
dev.off(dev.cur())

# get the dataframe cluster labels
df = getClusterMapping(oGr)
df = df[order(df$type.2),]
write.csv(df, file='Temp/df.csv')

temp = scan(what = character())
i = which(dfAnnotation$Entrez_Gene_ID %in% temp)
cvSelectedGenes = unique(as.character(dfAnnotation$Entrez_Gene_ID[i]))

# plot these genes
m = mExp.atb[,cvSelectedGenes]
m = as.data.frame(m)

x = stack(m)
x$f = f2
boxplot(values ~ f+ind, data=x, las=2, par=par(mar=c(8, 4, 2, 2)+0.1), main='Expression in pathway GPCR ligand binding', col=2:4)

########################################################
