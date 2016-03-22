# Name: ltbi_selection.R
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


# load the microarray data
# dfData.all = f_LoadObject(file.choose())
# # load sample annotation
# dfSamples = read.csv(file.choose(), header=T)
# # load the illumina annotation
# dfAnno = f_LoadObject(file.choose())
# 
# # create data set 
# fSamples = as.character(dfSamples$Illness1)
# # i = grepl('LTBI', fSamples, ignore.case = T)
# # fSamples[!i] = 'OD'
# # fSamples = factor(fSamples, levels = c('OD', 'LTBI'))
# 
# mData.all = t(dfData.all)
# lData = list(matrix=mData.all, groups=fSamples, annotation=dfAnno, desc='joint data sets to predict ltbi')
# save(lData, file='Objects/lData_with_annotation_for_ltbi.rds')

lData = f_LoadObject('Objects/lData_with_annotation_for_ltbi.rds')
# it is a list with various components
names(lData)
print(lData$desc)

## what do we want to predict, LTBI vs ATB
i = grep('ATB|LTBI', lData$groups, ignore.case = T)
lData$groups = lData$groups[i]
lData$matrix = lData$matrix[i,]

lData$groups = factor(lData$groups, levels = c('ATB', 'LTBI'))

## if a different factor is required
i = grepl('LTBI', lData$groups, ignore.case=T)
lData$groups[!i] = 'OD'
lData$groups = factor(lData$groups, levels = c('OD', 'LTBI'))

## remove sarcoidosis group only
i = grep('sarcoid', lData$groups, ignore.case = T)
lData$groups = lData$groups[-i]
lData$matrix = lData$matrix[-i,]

i = grepl('LTBI', lData$groups, ignore.case=T)
lData$groups[!i] = 'OD'
lData$groups = factor(lData$groups, levels = c('OD', 'LTBI'))

## make od=atb and healthy control group
## remove sarcoidosis group only
i = grep('Ctrl|ATB|LTBI', lData$groups, ignore.case = T)
lData$groups = lData$groups[i]
lData$matrix = lData$matrix[i,]

i = grepl('LTBI', lData$groups, ignore.case=T)
lData$groups[!i] = 'OD'
lData$groups = factor(lData$groups, levels = c('OD', 'LTBI'))

## get the components from the list
# test set index
test = sample(1:nrow(lData$matrix), 0.30 * nrow(lData$matrix), replace = F)
# sample annotation
fSamples = lData$groups
# expression matrix
mDat = lData$matrix
# annotation data
dfAnnotation = lData$annotation


#### choose a subset of genes of interest
cvEnt = c('3903', '780', '23556', '89872', '80729', '121260', '196527', '6095', '94039', '84246', '170960', '80321', '1260', '5108',
          '7284', '8481', '9026', '8685', '50837', '50865', '2852', '81127', '219429', '9934', '4283', '3362', '2865')
# get annotation
dfAnnotation.sel = dfAnnotation[as.character(dfAnnotation$Enterez_ID) %in% cvEnt,]
rownames(dfAnnotation.sel) = as.character(dfAnnotation.sel$ID)
il = as.character(dfAnnotation.sel$ID)
mDat = mDat[,colnames(mDat) %in% il]
dfAnnotation.sel = dfAnnotation.sel[colnames(mDat),]
## merge the duplicated probes
fDup = factor(as.character(dfAnnotation.sel$Symbol))
mDat.m = apply(mDat, 1, function(x) tapply(x, fDup, mean))
mDat = t(mDat.m)
cn = colnames(mDat)

mDat.test = mDat[test,]
mDat.train = mDat[-test,]
## perform nested random forest on test set
## adjust boot.num as desired
dfData = as.data.frame(mDat.train)
oVar.r = CVariableSelection.RandomForest(dfData, fSamples[-test], boot.num = 100)
# plot the top 20 genes based on importance scort with 95% confidence interval for standard error
plot.var.selection(oVar.r)
# get the variables
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 genes to find top combinations of genes
dfData = as.data.frame(mDat.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.sub = CVariableSelection.ReduceModel(dfData, fSamples[-test], boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 1:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
}

## 10 fold nested cross validation with various variable combinations
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:12){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = as.data.frame(mDat.train)
  dfData.train = data.frame(dfData.train[,colnames(dfData.train) %in% cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  
  dfData.test = as.data.frame(mDat.test)
  dfData.test = data.frame(dfData.test[,colnames(dfData.test) %in% cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fSamples[test],
                             train.groups = fSamples[-test], level.predict = 'LTBI', boot.num = 100)
  
  plot.cv.performance(oCV)
  
  # print variable names and 95% confidence interval for AUC
  x = getAUCVector(oCV)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}


