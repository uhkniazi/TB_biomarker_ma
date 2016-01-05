# Name: variable_selection_script.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 05/01/2016
# Desc: script for variable selection

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')

# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}

## load the test data
load(file='lData.publish.rds')
# it is a list with various components
names(lData)
print(lData$desc)

## get the components from the list
# test set index
test = lData$test
# sample annotation
dfSamples = lData$sample
# expression matrix
mDat = lData$expression
# annotation data
dfAnnotation = lData$annotation

# set variables
dfSamples.train = dfSamples[-test,]
dfSamples.test = dfSamples[test,]
mDat.train = mDat[-test,]
mDat.test = mDat[test,]

## perform nested random forest on test set
## adjust boot.num as desired
dfData = as.data.frame(mDat.train)
oVar.r = CVariableSelection.RandomForest(dfData, dfSamples.train$fGroups.2, boot.num = 100)
# plot the top 20 genes based on importance scort with 95% confidence interval for standard error
plot.var.selection(oVar.r)
# get the variables
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 genes to find top combinations of genes
dfData = as.data.frame(mDat.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.sub = CVariableSelection.ReduceModel(dfData, dfSamples.train$fGroups.2, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 2:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  print(paste('Variable Count', i))
  print(dfAnnotation[cvTopGenes.sub,])
}

## 10 fold nested cross validation with various variable combinations
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 2:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = as.data.frame(mDat.train)
  dfData.train = dfData.train[,colnames(dfData.train) %in% cvTopGenes.sub]
  
  dfData.test = as.data.frame(mDat.test)
  dfData.test = dfData.test[,colnames(dfData.test) %in% cvTopGenes.sub]
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
                             train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 500)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(dfAnnotation[cvTopGenes.sub,])
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}


