# Name: long_selected_genes_ATB_LTBI.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 30/06/15
# Desc: some genes of interest


source('../CCrossValidation/CCrossValidation.R')
p.old = par()
# load data from previous script
dfGenes.score = read.csv(file.choose(), header=T, row.names=1)
csFile = 'Objects/long_data_set_LTBI_VS_ATB_noCOV.rds'
load(csFile)
dfSamples = lData$sample
mDat.sub = lData$expression
test = lData$test

dfSamples.train = dfSamples[-test,]
dfSamples.test = dfSamples[test,]
mDat.sub.train = mDat.sub[-test,]
mDat.sub.test = mDat.sub[test,]

# choose the top 30 genes
cvTopGenes = rownames(dfGenes.score)[1:30]
cvTopGenes.sub = # insert gene name here

# get the genes from the training data
mDat = mDat.sub.train[,cvTopGenes]
# create correlation
mCor = cor(mDat)
as.data.frame(sort(mCor[,cvTopGenes.sub]))

i = which(abs(mCor[,cvTopGenes.sub]) < 0.6)
cvTopGenes = names(i)

# select model reduction
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]

oVar.sub = CVariableSelection.ReduceModel(dfData, dfSamples.train$fGroups.2)
plot.var.selection(oVar.sub)
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:9){
  cvTopGenes.min = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = as.data.frame(mDat.sub.train)
  dfData.train = dfData.train[,colnames(dfData.train) %in% c(cvTopGenes.sub, cvTopGenes.min)]
  
  dfData.test = as.data.frame(mDat.sub.test)
  dfData.test = dfData.test[,colnames(dfData.test) %in% c(cvTopGenes.sub, cvTopGenes.min)]
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
                             train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 100)
  
  plot.cv.performance(oCV)
}
par(p.old)

for (i in 1:9) {print(paste('size', i)); print(CVariableSelection.ReduceModel.getMinModel(oVar.sub, i))}

# plot using a tree based classifier
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:9){
  cvTopGenes.min = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = as.data.frame(mDat.sub.train)
  dfData.train = dfData.train[,colnames(dfData.train) %in% c(cvTopGenes.sub, cvTopGenes.min)]
  
  dfData.test = as.data.frame(mDat.sub.test)
  dfData.test = dfData.test[,colnames(dfData.test) %in% c(cvTopGenes.sub, cvTopGenes.min)]
  
  oCV = CCrossValidation.Tree(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
                             train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 50)
  
  plot.cv.performance(oCV)
}
par(p.old)

# testing 2 gene combinations
gene.cor = cvTopGenes
par(mfrow=c(2,2))
for (i in 1:length(gene.cor)){
  dfData.train = as.data.frame(mDat.sub.train)
  dfData.train = data.frame(dfData.train[,colnames(dfData.train) %in% c(cvTopGenes.sub, gene.cor[i])])
  dfData.test = as.data.frame(mDat.sub.test)
  dfData.test = data.frame(dfData.test[,colnames(dfData.test) %in% c(cvTopGenes.sub, gene.cor[i])])
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
                              train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 100)
  
  plot.cv.performance(oCV)  
}

as.data.frame((paste(cvTopGenes.sub, gene.cor)))






