# Name: june_4_Genelist1.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 04/06/15
# Desc: some new data for testing by long

source('tb_biomarker_ma_header.R')

## data loading
dfAnnotation = f_LoadObject(file.choose())
rownames(dfAnnotation) = as.character(dfAnnotation$ID)

# load the data, clean and create factors
dfExp = read.csv(file.choose(), header=T)

# load the sample annotation
dfSamples = dfExp[,c(1,2,3,4)]

# remove sarcoids
i = which(dfSamples$Disease %in% c('Non-active-sarcoidosis', 'Sarcoidosis'))
dfSamples = dfSamples[-i,]
dfExp = dfExp[-i,]

# create factors
fGroups = factor(dfSamples$Disease)
# create a second factor with only 2 levels
# keep ptb at 1 for downstream predictions
fGroups.2 = as.character(dfSamples$Disease)
fGroups.2[fGroups.2 != 'ATB'] = 'OD'
fGroups.2 = factor(fGroups.2, levels = c('OD', 'ATB'))
dfSamples$fGroups.2 = fGroups.2
dfSamples$fGroups = fGroups
rownames(dfSamples) = dfSamples$Probe_ID
rownames(dfExp) = dfExp$Probe_ID
dfExp = dfExp[,5:ncol(dfExp)]
dfExp = t(dfExp)
## data quality checks 
m = dfExp
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
fSamples = dfSamples$fGroups.2
col.p = c('red', 'blue')#rainbow(length(unique(fSamples)))
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
mDat.sub = t(mDat)
# keep 20% of data as test set
#test = sample(1:nrow(dfSamples), size = nrow(dfSamples) * 0.2, replace = F)
test = which(dfSamples$X == 'Test')
dfSamples.train = dfSamples[-test,]
dfSamples.test = dfSamples[test,]
mDat.sub.train = mDat.sub[-test,]
mDat.sub.test = mDat.sub[test,]
cvTopGenes = colnames(mDat.sub)
### model fitting and variable selection
## use random forest on training data for variable selection
dfData = as.data.frame(mDat.sub.train)
dfData$fGroups.2 = dfSamples.train$fGroups.2

rf.fit.1 = randomForest(fGroups.2 ~., data=dfData, importance = TRUE)

# save the results to save time for next time
dir.create('Objects', showWarnings = F)
save(rf.fit.1, file='Objects/long.rf.fit.1_2.rds')

# variables of importance
varImpPlot(rf.fit.1)
dfImportance = as.data.frame(importance(rf.fit.1))
dfImportance = dfImportance[order(dfImportance$MeanDecreaseAccuracy, decreasing = T),]
hist(dfImportance$MeanDecreaseAccuracy)
dfImportance.ATB = dfImportance[order(dfImportance$ATB, decreasing = T),]

# select the top few genes looking at the distribution of error rates
### TAG 1
# choose the top proteins for ATB
hist(dfImportance.ATB$ATB, main='Distribution of gene importance score in predicting ATB')
f = which(dfImportance.ATB$ATB >= 2)
length(f)
cvTopGenes = rownames(dfImportance.ATB)[f]

# subset the data based on these selected genes from training dataset
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
dfData$fGroups = dfSamples.train$fGroups.2

### Further variable classification check
### using CV and ROC
dfData.full = dfData
iBoot = 30

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
  ind.o.s = sample(ind.o, size = length(ind.p), replace = F)
  # sample of ATB groups, i.e. take everything as it is smaller group
  ind.p.s = sample(ind.p, size=length(ind.p), replace=F)
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
boxplot(df$MeanDecreaseAccuracy ~ df$Symbol, las=2, main='Gene scores for predicting ATB', cex.axis=0.3)
# calculate coefficient of variation
cv = ivSD/abs(ivMean)
# split data into groups based on cv
g = cut(cv, breaks = quantile(cv, 0:10/10), include.lowest = T)
par(p.old)
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

# annotation for these genes
dfAnnotation.genes = data.frame(dfAnnotation[cvTopGenes,], dfRF.boot.stats[cvTopGenes,])


## look at the correlation of the genes to remove colinear genes

dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
#dfData$fGroups = dfSamples.train$fGroups.2
mCor = cor(dfData)
i = findCorrelation(mCor, cutoff = 0.7)
n = colnames(mCor)[i]
# remove these correlated features 
cvTopGenes.step.2 = cvTopGenes
i = which(cvTopGenes %in% n)
cvTopGenes = cvTopGenes[-i]
rm(dfData)

## check for the miminum sized model using test and training sets
## use variable selection method
iBoot = 50
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
matplot(1:nrow(m), m, type='l', lty = 1, lwd=2, col=1:2, xaxt='n')
legend('bottomleft', legend = colnames(m), fill=1:2)
axis(1, at = 1:nrow(m), las=2)

## choose the best model after refitting on the full training data set
## choose which model is the best?
i = 3
#i = which.min(ivCV.test)[1]
# refit subset using i number of genes on all data
rm(list = c('dfData', 'dfData.train', 'dfData.test'))
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
dfData$fGroups = dfSamples.train$fGroups.2
# dfData = rbind(dfData.test, dfData.train)
reg = regsubsets(fGroups ~ ., data=dfData, nvmax = length(cvTopGenes), method='exhaustive')

# choose these variables
cvTopGenes.step.3 = cvTopGenes
cvTopGenes = names(coef(reg, i))[-1]
rm(list = c('dfData'))

### cross validation with ROC
#### CV with ROC
# choose all data together for nested 10 fold cv
par(p.old)
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
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
  # select a subset of equal numbers for others and SA
  ind.o = which(dfData.full$fGroups == 'OD')
  ind.p = which(dfData.full$fGroups == 'ATB')
  ind.o.s = sample(ind.o, size = length(ind.p), replace = T)
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
dfData.train = dfData.train[,colnames(dfData.train) %in% cvTopGenes]
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

## use the class created to do this
source('../CCrossValidation/CCrossValidation.R')
obj = CCrossValidation.LDA(test.dat = as.data.frame(mDat.sub.test[,cvTopGenes]),
                           train.dat = as.data.frame(mDat.sub.train[,cvTopGenes]), 
                           test.groups = dfSamples.test$fGroups.2, 
                           train.groups = dfSamples.train$fGroups.2, 
                           level.predict = 'ATB', boot.num=1000 )
plot.cv.performance(obj, 'bottomright')
x = obj@oAuc.cv@y.values
x = as.numeric(x)
