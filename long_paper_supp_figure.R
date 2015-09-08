# Name: long_paper_supp_figure.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 15/06/15
# Desc: analysis of all combined tb ma data, ATB vs OD to make figures

source('tb_biomarker_ma_header.R')
dfAnnotation = f_LoadObject(file.choose())
rownames(dfAnnotation) = as.character(dfAnnotation$ID)
load(file='Objects/long_data_set_noCOV.rds')
test = lData$test
lData$desc
dfSamples = lData$sample
mDat.sub = lData$expression
dfSamples.train = dfSamples[-test,]
dfSamples.test = dfSamples[test,]
mDat.sub.train = mDat.sub[-test,]
mDat.sub.test = mDat.sub[test,]
load(file='Objects/long.rf.fit.1_noCOV.rds')
varImpPlot(rf.fit.1)
dfImportance = as.data.frame(importance(rf.fit.1))
dfImportance = dfImportance[order(dfImportance$MeanDecreaseAccuracy, decreasing = T),]
hist(dfImportance$MeanDecreaseAccuracy)
dfImportance.ATB = dfImportance[order(dfImportance$ATB, decreasing = T),]
hist(dfImportance.ATB$ATB)
f = which(dfImportance.ATB$ATB >= 2)
length(f)
cvTopGenes = rownames(dfImportance.ATB)[f]
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.r = CVariableSelection.RandomForest(dfData, dfSamples.train$fGroups.2)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
cvTopGenes.step.1 = cvTopGenes
cvTopGenes = rownames(dfRF)[1:30]
dfData = as.data.frame(mDat.sub.train)
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.sub = CVariableSelection.ReduceModel(dfData, dfSamples.train$fGroups.2)
plot.var.selection(oVar.sub)
for (i in 2:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  print(paste('Variable Count', i))
  print(dfAnnotation[cvTopGenes.sub,])
}

par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 2:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = as.data.frame(mDat.sub.train)
  dfData.train = dfData.train[,colnames(dfData.train) %in% cvTopGenes.sub]
  
  dfData.test = as.data.frame(mDat.sub.test)
  dfData.test = dfData.test[,colnames(dfData.test) %in% cvTopGenes.sub]
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
                             train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 500)
  
  plot.cv.performance(oCV)
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(dfAnnotation[cvTopGenes.sub,])
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}


