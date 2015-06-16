# Name: long_all_data_test_03_ATB_LTBI.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 15/06/15
# Desc: analysis of all combined tb ma data

source('tb_biomarker_ma_header.R')

## data loading
# load the data, clean and create factors
dfExp = read.csv(file.choose(), header=T, row.names=1)

# load the sample annotation
dfSamples = read.csv(file.choose(), header=T)

# sort both the samples and expression data in same order
rownames(dfSamples) = as.character(dfSamples$Sample_ID)
dfSamples = dfSamples[colnames(dfExp),]

# load annotation data
dfAnnotation = f_LoadObject(file.choose())
rownames(dfAnnotation) = as.character(dfAnnotation$ID)

fGroups = factor(dfSamples$Illness1)
# keep only the LTBI and ATB groups
i = grep('ATB|LTBI', x = fGroups)
dfSamples = dfSamples[i,]
dfExp = dfExp[,i]
# create factors
fGroups = factor(dfSamples$Illness1)
# create a second factor with only 2 levels
# keep ptb at 1 for downstream predictions
fGroups.2 = as.character(dfSamples$Illness)
fGroups.2 = factor(fGroups.2, levels = c('OD', 'ATB'))
dfSamples$fGroups.2 = fGroups.2

## data quality checks 
m = dfExp
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
fSamples = dfSamples$fGroups.2
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
par(mfrow=c(2,2))
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components')
par(p.old)
# remove the outlier groups from the data
# these can be seen on the pc2 and pc3 plots
m = pr.out$x[,1:3]
m = data.frame(m, fSamples)
i = which(m$PC1 < -150 & m$PC2 > 50)
i = unique(c(i, which(m$PC2 > 150 & m$PC3 > 50)))
i = unique(c(i, which(m$PC3 > 150)))
c = col
c[i] = 'black'
## plot outlier groups
par(mfrow=c(2,2))
plot(pr.out$x[,1:2], col=c, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=c, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=c, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], c, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components')
par(p.old)

# remove the outliers
cvOutliers = rownames(m)[i]
m = match(cvOutliers, colnames(dfExp))
dfExp = dfExp[,-m]
m = match(cvOutliers, dfSamples$Sample_ID)
dfSamples = dfSamples[-m,]
gc()


### analysis
## extract data
# count matrix
mDat = as.matrix(dfExp)
# phenotypic data
fGroups.2 = dfSamples$fGroups.2

## data cleaning and formatting before statistical analysis
# remove any rows with NAN data
f = is.finite(rowSums(mDat))
table(f)
mDat = mDat[f,]
# select a subset of genes based on coefficient of variation.
mDat = t(mDat)
# calculate the coef of var for each gene in each group
# break data into 2 matrices based on groups
mDat.od = mDat[fGroups.2 == 'OD',]
mDat.atb = mDat[fGroups.2 == 'ATB',]
cv.od = apply(mDat.od, 2, function(x) sd(x)/abs(mean(x)))
cv.atb = apply(mDat.atb, 2, function(x) sd(x)/abs(mean(x)))
cv = apply(mDat, 2, function(x) sd(x)/abs(mean(x)))
# check cv
summary(cv.od)
summary(cv.atb)
summary(cv)
# cut data into groups based on quantiles of cv
cut.pts = quantile(cv, probs = 0:10/10)
groups = cut(cv, breaks = cut.pts, include.lowest = T, labels = 0:9)

cut.pts = quantile(cv.od, probs = 0:10/10)
groups.od = cut(cv.od, breaks = cut.pts, include.lowest = T, labels = 0:9)

cut.pts = quantile(cv.atb, probs = 0:10/10)
groups.atb = cut(cv.atb, breaks = cut.pts, include.lowest = T, labels = 0:9)


iMean = apply(mDat, 2, mean)
iVar = apply(mDat, 2, var)
coplot((cv) ~ iMean | groups)
coplot(iVar ~ iMean | groups)
coplot(iVar ~ iMean | groups.atb+groups.od)


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
# keep 30% of data as test set
test = sample(1:nrow(dfSamples), size = nrow(dfSamples) * 0.3, replace = F)
dfSamples.train = dfSamples[-test,]
dfSamples.test = dfSamples[test,]
mDat.sub.train = mDat.sub[-test,]
mDat.sub.test = mDat.sub[test,]

lData = list(test=test, sample=dfSamples, expression=mDat.sub)
lData$desc = 'Longs data set for including test set vector LTBI vs ATB - no COV'
save(lData, file='Objects/long_data_set_LTBI_VS_ATB_noCOV.rds')
### model fitting and variable selection
## use random forest on training data for variable selection
dfData = as.data.frame(mDat.sub.train)
dfData$fGroups.2 = dfSamples.train$fGroups.2

rf.fit.1 = randomForest(fGroups.2 ~., data=dfData, importance = TRUE)

# save the results to save time for next time
dir.create('Objects', showWarnings = F)
save(rf.fit.1, file='Objects/long.rf.fit.1_LTBI_VS_ATB_noCOV.rds')

# variables of importance
varImpPlot(rf.fit.1)
dfImportance = as.data.frame(importance(rf.fit.1))
dfImportance = dfImportance[order(dfImportance$MeanDecreaseAccuracy, decreasing = T),]
hist(dfImportance$MeanDecreaseAccuracy)
dfImportance.ATB = dfImportance[order(dfImportance$ATB, decreasing = T),]

# select the top few genes looking at the distribution of error rates
### TAG 1
# choose the top proteins for ATB
hist(dfImportance.ATB$ATB)
f = which(dfImportance.ATB$ATB >= 1.5)
length(f)
cvTopGenes = rownames(dfImportance.ATB)[f]

# subset the data based on these selected genes from training dataset
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
dfData$fGroups = dfSamples.train$fGroups.2

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
  ind.o = which(dfData.full$fGroups == 'OD')
  ind.p = which(dfData.full$fGroups == 'ATB')
  # take sample of equal in size from group OD and ATB
  sam.size = min(length(ind.p), length(ind.o))
  ind.o.s = sample(ind.o, size = sam.size, replace = F)
  # sample of ATB groups, i.e. take everything as it is smaller group
  ind.p.s = sample(ind.p, size=sam.size, replace=F)
  # join the sample indices together
  ind = sample(c(ind.o.s, ind.p.s), replace=F)
  # take sample from the full dataset
  dfData = dfData.full[ind,]
  # fit model
  fit.rf = randomForest(fGroups ~., data=dfData, importance = TRUE, ntree = 500)
  # get variables importance
  df = importance(fit.rf)
  df = df[order(df[,'MeanDecreaseAccuracy'], decreasing = T),]
  # put in list
  lVarImp[[i]] = df
} # for

## put data for each boot of each variable together in a dataframe
df = NULL
for (i in 1:iBoot) df = rbind(df, lVarImp[[i]])
# convert rownames i.e. gene names to factors
f = as.factor(rownames(df))
# calculate mean and sd for each gene
ivMean = tapply(df[,'MeanDecreaseAccuracy'], f, mean)
ivSD = tapply(df[,'MeanDecreaseAccuracy'], f, sd)
df = as.data.frame(df)
df$Symbol = rownames(df)
dfRF.boot = df
# boxplots 
par(mar=c(6,4,3,2)+0.1)
boxplot(df$MeanDecreaseAccuracy ~ df$Symbol, las=2, cex.axis=0.4)
# calculate coefficient of variation
cv = ivSD/abs(ivMean)
# split data into groups based on cv
g = cut(cv, breaks = quantile(cv, 0:10/10), include.lowest = T)
coplot(ivSD ~ ivMean | g)
gl = cut(cv, breaks = quantile(cv, 0:10/10), include.lowest = T, labels = 0:9)
rm(dfData)
rm(dfData.full)
par(p.old)
dfRF.boot.stats = data.frame(ivMean, ivSD, cv, groups=g, group.lab=gl)
## Decide on a cutoff here
## based on coefficient of variation
cvTopGenes.step.1 = cvTopGenes
f = cv[gl %in% c(0, 1)]
cvTopGenes = names(f)

## check for the miminum sized model using test and training sets
## use variable selection method
iBoot = 1000
mTrain = matrix(NA, nrow = length(cvTopGenes), ncol = iBoot)
mTest = matrix(NA, nrow = length(cvTopGenes), ncol = iBoot)

for(o in 1:iBoot){
  dfData.train = as.data.frame(mDat.sub.train)
  dfData.train = dfData.train[,colnames(dfData.train) %in% cvTopGenes]
  dfData.train$fGroups = dfSamples.train$fGroups.2
  # create a test set on a percentage the data
  test = sample(1:nrow(dfData.train), size =nrow(dfData.train) * 0.30, replace = F)
  dfData.test = dfData.train[test,]
  dfData.train = dfData.train[-test,]
  
  # fit model
  reg = regsubsets(fGroups ~ ., data=dfData.train, nvmax = length(cvTopGenes), method='exhaustive')
  # test for validation errors in the test set
  ivCV.train = rep(NA, length=length(cvTopGenes))
  ivCV.test = rep(NA, length=length(cvTopGenes))
  for (i in 1:length(cvTopGenes)){
    # get the genes in each subset
    n = names(coef(reg, i))[-1]
    n = c(n, 'fGroups')
    dfDat.train = dfData.train[,colnames(dfData.train) %in% n]
    dfDat.test = dfData.test[,colnames(dfData.test) %in% n]
    # fit the lda model on training dataset
    fit.lda = lda(fGroups ~ ., data=dfDat.train)
    # test error rate on test dataset
    p = predict(fit.lda, newdata=dfDat.test)
    # calculate test error 
    ivCV.test[i] = mean(p$class != dfDat.test$fGroups)  
    # calculate training error
    p = predict(fit.lda, newdata=dfDat.train)
    # calculate error
    ivCV.train[i] = mean(p$class != dfDat.train$fGroups)  
  }
  mTrain[,o] = ivCV.train
  mTest[,o] = ivCV.test
  print(o)
}
# make plots
mTrain = t(mTrain)
par(mfrow=c(1,2), mar=c(4,3,1,1))
boxplot(mTrain, main='Training Set', xlab='No. of variables', ylab='Error rate', ylim=c(0.1, 0.23))
lines(1:ncol(mTrain), colMeans(mTrain), col='red', lwd=2)

mTest = t(mTest)
boxplot(mTest, main='Test Set', xlab='No. of variables', ylab='Error rate', ylim=c(0.1, 0.23))
lines(1:ncol(mTest), colMeans(mTest), col='red', lwd=2)
y = colMeans(mTest)
abline(h = min(y), lty=2)
# calculate standard error for minimum mean test error column
i = which.min(y)
x = mTest[,i]
se = sd(x)/sqrt(length(x))
abline(h = mean(x)+se+se, lty=2, col=2)
abline(h = mean(x)-se-se, lty=2, col=2)
abline(h = mean(x)+sd(x), lty=2, col='blue')
abline(h = mean(x)-sd(x), lty=2, col='blue')

par(p.old)
# plot together
tr = colMeans(mTrain)
te = colMeans(mTest)

m = cbind(tr, te)
matplot(1:nrow(m), m, type='b', pch=20, lty = 1, lwd=2, col=1:2, xaxt='n')
legend('bottomleft', legend = colnames(m), fill=1:2)
axis(1, at = 1:nrow(m), las=2)

## choose the best model after refitting on the full training data set
## choose which model is the best?
i = 2
# refit subset using i number of genes on all data
rm(list = c('dfData.train', 'dfData.test'))
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
dfData$fGroups = dfSamples.train$fGroups.2
reg = regsubsets(fGroups ~ ., data=dfData, nvmax = length(cvTopGenes), method='exhaustive')

# choose these variables
cvTopGenes.step.3 = cvTopGenes
cvTopGenes = names(coef(reg, i))[-1]
rm(list = c('dfData', 'dfDat.test', 'dfDat.train'))

######### cross validation with ROC
dfData.train = as.data.frame(mDat.sub.train)
dfData.train = dfData.train[,colnames(dfData.train) %in% cvTopGenes]

dfData.test = as.data.frame(mDat.sub.test)
dfData.test = dfData.test[,colnames(dfData.test) %in% cvTopGenes]

oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
                           train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 500)

plot.cv.performance(oCV)

## plot these expression values for these genes
par(mfrow=c(1,2))
x = stack(dfData.train)
x$f = dfSamples.train$fGroups.2
boxplot(values ~ f+ind, data=x, las=2, par=par(mar=c(8, 4, 2, 2)+0.1), main='Training Data')

x = stack(dfData.test)
x$f = dfSamples.test$fGroups.2
boxplot(values ~ f+ind, data=x, las=2, par=par(mar=c(8, 4, 2, 2)+0.1), main='Test Data')

dfAnnotation[cvTopGenes,]


################################################################################
#### Pathway comparisons
source('../CGraphClust/CGraphClust.R')
library(igraph)
library(org.Hs.eg.db)
library(reactome.db)

i = which(dfSamples.ns$Illness1 %in% c('LTBI', 'ATB'))
fGroups = as.character(dfSamples.ns$Illness1[i])
mDat.2 = mDat.ns[i, ]
fGroups = factor(fGroups)

# get names of DE genes from t and wilc tests
t = names(p.t.adj[p.t.adj < 0.001])
w = names(p.w.adj[p.w.adj < 0.001])
n = unique(c(w, t))
f1 = n %in% t
f2 = n %in% w
f = f1 & f2
n2 = n[f]
mDat.2 = mDat.2[,n2]

# replace gene names i.e. illumina ids by enterez ids
c = colnames(mDat.2)
c2 = as.character(dfAnnotation[c,'Enterez_ID'])
# which names are NA
i = which(is.na(c2))
mDat.2 = mDat.2[,-i]

c = colnames(mDat.2)
c2 = as.character(dfAnnotation[c,'Enterez_ID'])
i = !duplicated(c2)
mDat.2 = mDat.2[,i]

c = colnames(mDat.2)
c2 = as.character(dfAnnotation[c,'Enterez_ID'])
colnames(mDat.2) = c2

# create data frame for graph
dfGraph = AnnotationDbi::select(reactome.db, colnames(mDat.2), 'REACTOMEID', 'ENTREZID')
dfGraph = na.omit(dfGraph)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mDat.2[,n]
mCor = cor(mCounts)

hist(mCor, prob=T)

# create the graph
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.7)#, iCorCut = 0.7)

# order the count matrix before making heatmaps
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomright')
plot.heatmap.means(oGr, t(mCounts))
plot.heatmap(oGr, t(mCounts))

ig = getFinalGraph(oGr)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)

df = getClusterMapping(oGr)
colnames(df) = c('gene', 'cluster')
df = df[order(df$cluster),]

write.csv(df, 'Data_results/long_all_data_test_03_pathway.csv')